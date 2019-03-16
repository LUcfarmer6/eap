using ArgParse;
using Optim;
using HCubature;
using SpecialFunctions;
using PyPlot;

include(joinpath("inc", "std.jl"));

s = ArgParseSettings();
@add_arg_table s begin
  "--omega", "-w"
    help = "dimensionless electric energy"
    arg_type = Float64
    default = 1.0
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
const ω = Complex(pa["omega"]); # to work with imaginary error function
const N = pa["num-monomers"];
const iters = pa["num-iterations"];

p = P(1.0, pi/4);
@show λrelax = find_min_λ(p, N, 0.0; iters=iters);
@show q = Q(p, λrelax);
@show find_min_λ(q, N, 0.0; iters=iters);
@show results = optimize(λ -> begin;
        (val, err) = hquadrature(θ -> p(θ)*A(map_forward(θ, λ), λ, N, ω), 0, pi/2);
      return val;
  end, 0.1, 10.0; iterations=iters);
@show λelec = find_min_λ(q, N, ω; iters=iters);
@show r = Q(q, λelec);

θs = 0.0:0.01:pi/2;
plot(θs, map(p, θs), "r-"; label="prestressed");
plot(θs, map(q, θs), "b--"; label="relaxed");
plot(map(θ -> map_forward(θ, λrelax), θs), map(p, θs), "g."; label="relaxed, test");
plot(θs, map(r, θs), "k-."; label="electric");
legend();
show();
