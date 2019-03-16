function Li(y::Real)
  return (3*y - y/5 * (6*y^2 + y^4 - 2*y^6))/(1 - y^2);
end

modified_vonmises(κ::Real, θ0::Real, θ::Real) = exp(κ*cos(4*(θ-θ0))) * sin(θ);

struct P
  κ::Real;
  θ0::Real;
  C::Real;

  function P(κ::Real, θ0::Real)
    if κ < 0; error("κ must be nonnegative"); end;
    if θ0 < 0 || θ0 > pi/2; error("θ0 must be between 0 and π/2"); end;

    (val, err) = hquadrature(θ -> modified_vonmises(κ, θ0, θ), 0, pi/2);
    if err > 1.; error("Integration failed."); end;
    return new(κ, θ0, 1.0 / val);
  end

end;

struct Q
  pq::Union{P, Q}; # distriubtion to pullback too
  λ::Real;
end;

(p::P)(θ::Real) = p.C * modified_vonmises(p.κ, p.θ0, θ);

function (q::Q)(θ::Real)
  # atan(sqrt(λ)*sin(θ) / (cos(θ) / λ))
  # map theta back to reference then ask what its probability is
  return q.pq(atan(q.λ*sqrt(q.λ)*tan(θ)));
end

function A(θ::Real, λ::Real, N::Int, ω::Number)

  r = [sin(θ) / sqrt(λ); 0; λ*cos(θ)];
  γ = sqrt( ((sin(θ))^2/λ + λ*λ*(cos(θ))^2) / N );
  rhat = r / sqrt( ((sin(θ))^2/λ + λ*λ*(cos(θ))^2) );
  rz2 = rhat[3]*rhat[3];
  Lγ = Li(γ);
  sqrtω = (ω != Complex(0.0)) ? sqrt(Complex(ω)) : 0.0;

  ret =   (ω != 0.0) ?
          real( ω*rz2 + γ*ω / Lγ * (1 - 3*rz2) + γ*Lγ + log(Lγ / (4*pi*sinh(Lγ)))
             + (1 - γ*γ)*(-ω/3 + log(2 * sqrtω / ( sqrt(pi) * erf(sqrtω) ) ) )
            ) :
          (γ*Lγ + log(Lγ / (4*pi*sinh(Lγ)))
          );

  # make sure A comes out pure real
  @assert(abs(imag(ret)) < 1e-10);
  return real(ret);

end

function find_min_λ(pq::Union{P, Q}, N::Int, ω::Number; lower::Real=0.1,
                    upper::Real=10.0, iters::Int=10000)
  results = optimize(λ -> begin;
    (val, err) = hquadrature(θ -> pq(θ)*A(θ, λ, N, ω), 0, pi/2);
    return val;
  end, lower, upper; iterations=iters);
  return Optim.minimizer(results);
end

function map_forward(θr::Real, λ::Real)
  return atan(tan(θr) / (λ*sqrt(λ)));
end

struct Ρ
  ω::Real;
  C::Real;
  λ::Real;
  α::Real;
end;

struct U
  ω::Real;
  kT::Real;
  u0::Real;

  U(ω::Real; kT::Real=1.0) = new(ω, kT, (ω > 0) ? ω : 0.0);
end

(ρ::Ρ)(ϕ::Real, θ::Real) = ρ.C * exp(Ρ.ω*(cos(θ))^2 + ρ.λ*cos(θ) + ρ.α*cos(ϕ)*sin(θ));
(ρ::Ρ)(θs::Vector) = ρ(θs[1], θs[2]);

(u::U)(ϕ::Real, θ::Real) = u.kT * u.ω * (cos(θ))^2 + _u0_;
(u::U)(θs::Vector) = u(θs[1], θs[2]);

# solvers
integrand_functions(xs::Vector) = (
  [
    (θs::Vector) -> ρ(θs, xs)*sin(θs[2]),
    (θs::Vector) -> ρ(θs, xs)*sin(θs[2])*cos(θs[2]),
    (θs::Vector) -> ρ(θs, xs)*sin(θs[2])*sin(θs[2])*cos(θs[1])
  ]
);

_LOWER_BOUNDS_ = [0.0, 0.0];
_UPPER_BOUNDS_ = [2*pi, pi];

function integrate_sys(xs::Vector)
  return map(integrand -> hcubature(integrand, _LOWER_BOUNDS_, _UPPER_BOUNDS_),
             integrand_functions(xs));
end

function gradc(xs::Vector)
  grad = zeros(3);
  grad[1], = hcubature(θs -> ρ(θs, xs)/xs[1]*sin(θs[2]), _LOWER_BOUNDS_, _UPPER_BOUNDS_);
  grad[2], = hcubature(θs -> ρ(θs, xs)/xs[1]*sin(θs[2])*cos(θs[2]), _LOWER_BOUNDS_, _UPPER_BOUNDS_);
  grad[3], = hcubature(θs -> ρ(θs, xs)/xs[1]*sin(θs[2])*sin(θs[2])*cos(θs[1]), _LOWER_BOUNDS_, _UPPER_BOUNDS_);

  return grad;
end

function gradλ(xs::Vector)
  global RELTOL_INT;
  global MAXEVALS_INT;
  
  grad = zeros(3);
  grad[1], = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*cos(θs[2]), _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  grad[2], = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*cos(θs[2])*cos(θs[2]), _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  grad[3], = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*sin(θs[2])*cos(θs[1])*cos(θs[2]), _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);

  return grad;
end

function gradα(xs::Vector)
  global RELTOL_INT;
  global MAXEVALS_INT;
  
  grad = zeros(3);
  grad[1], = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*sin(θs[2])*cos(θs[1]), _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  grad[2], = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*cos(θs[2])*sin(θs[2])*cos(θs[1]), _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  grad[3], = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*(sin(θs[2])*cos(θs[1]))^2, _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);

  return grad;
end

function hessian(xs::Vector)
  h = zeros(3, 3);

  h[:, 1] = gradc(xs);
  h[:, 2] = gradλ(xs);
  h[:, 3] = gradα(xs);

  return h;
end

function newton_raphson_iteration(xs::Vector, rs::Vector)
  h = hessian(xs);
  return xs - inv(h) * rs;
end

function newton_raphson(xs::Vector; reltol::Real=1e-10, abstol::Real=1e-8, max_iters::Int=100000)
  initial_Is = integrate_sys(xs);
  current_residuals = residuals(initial_Is);
  residual_norm = norm(current_residuals);
  init_residual_norm = residual_norm;
  iter = 1;
    
  if residual_norm < abstol
    return (xs, current_residuals, true, iter);
  end

  for iter = 2:max_iters

    lin_alg_error_flag = false;
    try
      xs = newton_raphson_iteration(xs, current_residuals);
    catch
      lin_alg_error_flag = true;
    end

    # if Hessian is singular, or Newton iteration fails, perturb guess
    if lin_alg_error_flag
      xs += rand(3);
    end

    Is = integrate_sys(xs);
    current_residuals = residuals(Is);
    residual_norm = norm(current_residuals);

    relnorm = residual_norm / init_residual_norm;

    if residual_norm < abstol || relnorm < reltol
      return (xs, current_residuals, true, iter);
    end

    @update_user(begin
      global _NUMEVALS_;
      println("Newton iterations: $_NUMEVALS_");
      println("    current residuals: ", current_residuals);
      println("    Is: ", [Is[1][1], Is[2][1], Is[3][1]]);
    end);
  end

  return (xs, current_residuals, false, iter);
end

function simulated_annealing(xs::Vector; anneal_factor::Real=1e-3,
                             init_T::Real=10.0, final_T::Real=1e-7,
                             relax_steps::Int=100,
                             init_delta::Real=0.5,
                             delta_adapt_factor::Real=1.2,
                             adapt_acc_thrshld::Real=0.7,
                             reltol::Real=1e-10, abstol::Real=1e-8, 
                             max_iters::Int=1000000)

  initial_Is = integrate_sys(xs);
  current_residuals = residuals(initial_Is);
  residual_norm = norm(current_residuals);
  init_residual_norm = residual_norm;
  iter = 1;
  Is = initial_Is;
    
  if residual_norm < abstol
    return (xs, current_residuals, true, iter);
  end

  acc = 0;
  delta = init_delta;
  T = init_T;

  for iter = 2:max_iters

    # check to see if relaxation should occur
    if iter % relax_steps == 0
      T *= (1.0 - anneal_factor);
      if T < final_T
        break;
      end

      # adapt step size to current acceptance ratio
      if acc / relax_steps > adapt_acc_thrshld
        delta *= delta_adapt_factor;
      elseif acc / relax_steps < (1 - adapt_acc_thrshld)
        delta /= delta_adapt_factor;
      end

      # reset acceptance count
      acc = 0;
    end

    # trial move
    choice = rand(1:3);
    dx = (rand(Float64)*2 - 1) * delta;
    new_xs = copy(xs);
    new_xs[choice] += dx;

    new_Is = integrate_sys(new_xs);
    new_residuals = residuals(new_Is);

    # metropolis acceptance criteria
    if rand(Float64) < exp(-(dot(new_residuals, new_residuals) - 
                             dot(current_residuals, current_residuals)) / T)

      # accept move
      acc += 1;
      copy!(xs, new_xs);
      current_residuals = new_residuals;
      Is = copy(new_Is);

      # check for convergence
      residual_norm = norm(current_residuals);
      relnorm = residual_norm / init_residual_norm;
      if residual_norm < abstol || relnorm < reltol
        return (xs, current_residuals, true, iter);
      end
      
    end
   
    @update_user(begin
      global _NUMEVALS_;
      println("Simulated annealing iterations: $_NUMEVALS_");
      println("    current residuals: ", current_residuals);
      println("    Is: ", [Is[1][1], Is[2][1], Is[3][1]]);
      println("    T: ", T);
      println("    delta: ", delta);
    end);
  end

  return (xs, current_residuals, false, iter);

end

