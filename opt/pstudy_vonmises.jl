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

figr = figure();
fige = figure();

@time for ω in [-10.0; -5.0; -1.0; -0.5; 0.0; 0.5; 1.0; 5.0; 10.0]
  @show ω
  #κs = 0.0:0.1:10.0;
  #θ0s = pi/24:pi/48:11*pi/24;
  κs = 0.0:0.25:5.0;
  θ0s = pi/24:pi/24:11*pi/24;
  nκ, nθ = length(κs), length(θ0s);
  λrelaxs = zeros(nκ, nθ);
  λelecs = zeros(nκ, nθ);
  @time for i = 1:nκ, j = 1:nθ

    κ = κs[i];
    θ0 = θ0s[j];
    
    p = P(κ, θ0);
    @show λrelax = find_min_λ(p, N, 0.0; iters=iters);
    @show q = Q(p, λrelax);
    @show find_min_λ(q, N, 0.0; iters=iters);
    @show results = optimize(λ -> begin;
            (val, err) = hquadrature(θ -> p(θ)*A(map_forward(θ, λ), λ, N, ω), 0, pi/2);
          return val;
      end, 0.1, 10.0; iterations=iters);
    @show λelec = find_min_λ(q, N, ω; iters=iters);

    λrelaxs[i, j] = λrelax;
    λelecs[i, j] = λelec;
  end

  surf(θ0s, κs, λrelaxs; label="\$\\omega = $ω\$", figure=figr);
  xlabel("\$\\theta_0\$");
  ylabel("\$\\kappa\$");
  zlabel("\$\\lambda\$");
  title("\$\\lambda_{r}\$");

  surf(θ0s, κs, λelecs; label="\$\\omega = $ω\$", figure=fige);
  xlabel("\$\\theta_0\$");
  ylabel("\$\\kappa\$");
  zlabel("\$\\lambda\$");
  title("\$\\lambda_{e}\$");
end

savefig(joinpath("figs", "lambdar.pdf"); figure=figr);
#savefig(joinpath("figs", "lambdar.png"); figure=figr);
savefig(joinpath("figs", "lambdae.pdf"); figure=fige);
#savefig(joinpath("figs", "lambdae.png"); figure=fige);
