fω(E0::Real, Δκ::Real, β::Real) = β * E0^2 * Δκ;
fω(E0::Vector, Δκ::Real, β::Real) = ω(dot(E0, E0), Δκ, β);
fω0(E0::Real, κ2::Real, β::Real) = β * E0^2 * κ2;
fω0(E0::Vector, κ2::Real, β::Real) = ω0(dot(E0, E0), κ2, β);
fΔκ(κ1::Real, κ2::Real) = κ2 - κ1;

#! Inverse Langevin function
function finvL(x::Real)
  x2 = x*x;
  x4 = x2*x2;
  x6 = x2*x4;

  return (3*x - x / 5 * (6*x2 + x4 - 2*x6)) / (1 - x2);
end

#! Derivative of the Inverse Langevin function
function fDinvL(x::Real)
  x2 = x*x;
  x4 = x2*x2;

  return 1 + 1.0/((x-1)^2) - 3/5*x2 - 2*x4 + 1/((x+1)^2);
end

#! Free energy in small lambda limit
function fAsl(E0::Real, κ1::Real, κ2::Real, γ::Real, ψ::Real, β::Real, n::Real, 
              N0::Real)

  const Δκ = fΔκ(κ1, κ2);
  const ω  = fω(E0, Δκ, β);
  const ω0  = fω0(E0, Δκ, β);
  const expω = exp(ω);
  const sqrtω = sqrt(ω);
  const erfsω = erf(sqrtω);
  const sqrtπ = sqrt(π);

  const term1 = (2 / (2*ω - 1) * 
                 (ω0 - ω0*ω + γ^2 * ω * (1 + 2*ω) + (1 - 2*ω) * log(2*n)));
  
  const subterm1 = (2*ω - 3) * cos(2 * ψ) / (2*ω - 1);
  const subsubterm1 = (cos(ψ))^2 / (expω*sqrtπ*erfsω - 2*sqrtω);
  const subsubterm2 = ((2*sin(ψ))^2 / 
                       ((2*ω-1)*(2*sqrtω + expω*sqrtπ*erfsω*(2*ω-1))));

  return n*N0*(term1 - 3*log(π) + 2*log(n*sqrtω/erfsω) + 
               2*γ^2*ω*(subterm1 + 4*sqrtω*(subsubterm1 - subsubterm2))
              ) / (2.0 * β);
end

#! Free energy near the stretched limit
function fAkg(E0::Real, κ1::Real, κ2::Real, γ::Real, ψ::Real, β::Real, n::Real, 
              N0::Real)

  const ω = fω(E0, Δκ, β);
  const ω0 = fω0(E0, Δκ, β);
  const Λ = finvL(γ);

  return n*N0/(β*Λ)*( (Λ-2*γ)*ω*(cos(ψ))^2 + 
                     Λ*(ω0 + γ*Λ + log(Λ / (4*pi*sinh(Λ)))) + 
                     γ*ω*(sin(ψ))^2
                    );

end

#! Maximum principal stress
function fσ(E0::Real, Δκ::Real, β::Real, n::Real, N0::Real, a::Real)

  const kT = 1.0 / β;
  const ω = fω(E0, Δκ, β);
  const a2 = a*a;
  const a3 = a*a2;
  const a6 = a3*a3;
  const a9 = a6*a3;
  const λc = sqrt( (2+a3)/(a) );
  const λcr = λc / sqrt( (3*n*a) );
  const Λ = finvL(λcr);
  const dΛ = fDinvL(λcr);
  const sqrtn = sqrt(n);
  const sqrtω = sqrt(ω);
  const sqrt3 = sqrt(3);
  const term4 = (-4 + 3*a6 + a9);
  const term2 = (2 - 3*a3 + a9);
  const term1 = (1 + 7*a3 + a6);
  const logπ4 = log(π/4);

  const mult = kT*λc^(3/2) / (9 * (2+a3)^4 * Λ^2);
  sum = 3*sqrt3*term4*sqrtn*Λ^3 + 6*λc*term2*ω*dΛ;
  sum -= 3*sqrt3*(2+a3)*sqrtn*Λ*(2*term1*ω - (-2 + a3 + a6)*dΛ);
  sum += Λ^2 * (λc* (2*(-4+3*a6+a9+27*a2*a2*n)*ω + 6*term4*log(n) +
                     3*term4*logπ4 - 6*term4*log(n*sqrtω/erf(sqrtω))) +
                3*(a3-1)*(2*a3)^2*(λc - sqrtn*sqrt3*coth(Λ))*dΛ);

  return mult * sum;
end
