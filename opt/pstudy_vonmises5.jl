using ArgParse;
using Optim;
using HCubature;
using SpecialFunctions;
using PyPlot;
import Printf;

include(joinpath("inc", "std.jl"));

s = ArgParseSettings();
@add_arg_table s begin
  "--num-monomers", "-N"
    help = "number of monomers"
    arg_type = Int
    default = 100
  "--num-iterations", "-O"
    help = "number of iterations"
    arg_type = Int
    default = 10000
end

pa = parse_args(s);
const N = pa["num-monomers"];
const iters = pa["num-iterations"];

@time for ω in [-10.0; 10.0]
  @show ω
  #κs = 0.0:0.1:10.0;
  #θ0s = pi/24:pi/48:11*pi/24;
  κs = 0.0:0.25:5.0;
  θ0s = pi/12:pi/24:5*pi/12;
  nκ, nθ = length(κs), length(θ0s);
  λrelaxs = zeros(nκ, nθ);
  λfinals = zeros(nκ, nθ);
  λelecs = zeros(nκ, nθ);
  ΔAs = zeros(nκ, nθ);
  @time for i = 1:nκ, j = 1:nθ

    κ = κs[i];
    θ0 = θ0s[j];
    
    p = P(κ, θ0);
    @show λrelax = find_min_λ(p, N, 0.0; iters=iters);
    @show λfinal = find_min_λ(p, N, ω; iters=iters);
    @show λelec = λfinal / λrelax;

    λrelaxs[i, j] = λrelax;
    λfinals[i, j] = λfinal;
    λelecs[i, j] = λelec;
    ΔA = (hquadrature(θ -> p(θ)*A(θ, λrelax, N, 0), 0, pi/2)[1]
          - hquadrature(θ -> p(θ)*A(θ, λfinal, N, ω), 0, pi/2)[1]);
    ΔAs[i, j] = ΔA;
  end

  surf(θ0s, κs, λelecs);
  xlabel("\$\\theta_0\$");
  ylabel("\$\\kappa\$");
  zlabel("\$\\lambda\$");
  show();
  clf();

  surf(θ0s, κs, ΔAs);
  xlabel("\$\\theta_0\$");
  ylabel("\$\\kappa\$");
  zlabel("\$\\Delta A\$");
  show();
  clf();
end
