using ArgParse;

include(joinpath("inc", "algorithm_help_msg.jl"));

s = ArgParseSettings();
@add_arg_table s begin
  "--algorithm", "-A"
    help = "algorithm (ISRES | PRAXIS | SBPLX | SLSQP | LBFGS | AUGLAG | Newton | Anneal). Run with --print-algo-help for more detailed information"
    default = "try Newton then SBPLX,Anneal,PRAXIS"
  "--maxevals-opt", "-E"
    help = "maximum evaluations for numerical integration"
    arg_type = ASCIIString
    default = "try 100 then 1000,10000,500000"
  "--dk", "-D"
    help = "kappa2 - kappa1"
    arg_type = Float64
    default = 1.0
  "--kT", "-k"
    help = "dimensionless temperature"
    arg_type = Float64
    default = 1.0
  "--num-monomers", "-N"
    help = "number of monomers"
    arg_type = Int
    default = 100
  "--spring-constant", "-k"
    help = "dimensionless stretch"
    arg_type = Float64
    default = 1.0
  "--print-every-percent-compl"
    help = "Print diagnostic info every percentage complete of maximum iterations (should be between 0 and 1)"
    arg_type = Float64
    default = 0.1
  "--verbose", "-v"
    help = "print diagnostic information"
    action = :store_true
  "--reltol-int", "-r"
    help = "relative tolerance for numerical integration"
    arg_type = Float64
    default = 1e-6
  "--maxevals-int", "-e"
    help = "maximum evaluations for numerical integration"
    arg_type = Int
    default = 2500
  "--evals", "-E"
    help = "iterations to explore residual surface"
    arg_type = Int
    default = 1000000
  "--init-x", "-X"
    help = "initial x"
    arg_type = ASCIIString
  "--temperature", "-t"
    help = "simulation temperature"
    arg_type = Float64
    default = 10.0
  "--delta", "-d"
    help = "maximum step size"
    arg_type = Float64
    default = 5.0
  "--output-steps", "-S"
    help = "simulation temperature"
    arg_type = Int
    default = 10
  "--outfile", "-F"
    help = "output file"
    default = "energy_surface.csv"
end

pa = parse_args(s);
if pa["print-algo-help"]
  println(_ALGO_HELP_MSG_);
  exit();
end
const dk = pa["dk"];
const kT = pa["kT"];
const N = pa["num-monomers"];

const verbose = pa["verbose"];
const print_every_percent_compl = pa["print-every-percent-compl"];
const RELTOL_INT = pa["reltol-int"];
const MAXEVALS_INT = pa["maxevals-int"];
const T = pa["temperature"];
const delta = pa["delta"];
const stepout = pa["output-steps"];
const outfile = pa["outfile"];
const iters = pa["evals"];

include(joinpath("inc", "solvers.jl"));
include(joinpath("inc", "free_energy.jl"));
include(joinpath("inc", "output.jl"));

_RHS_ = _RHS_GCOORD_;
ω(ϕ::Real, θ::Real) = _ω_gcoord_(ϕ, θ);

function calculate_energy(xs::Vector)
  main_algos = parse_algostr(pa["algorithm"]);
  main_max_iters = parse_max_iters(pa["maxevals-opt"]);
  (xs, rs, converged, iterations) = run_solver(main_algos, main_max_iters;
                                               x0s = DEFAULT_INIT_GUESS,
                                               reltol = reltol_opt,
                                               abstol = abstol_opt);
end

for iter = 1:iters

  # trial move
  choice = rand(1:3);
  dx = (rand(Float64)*2 - 1) * delta;
  new_xs = copy(xs);
  new_xs[choice] += dx;

  new_Is = integrate_sys(new_xs);
  new_residuals = residuals(new_Is);

  # metropolis acceptance criteria
  if rand(Float64) < exp(-(norm(new_residuals, 2) - 
                           residual_norm) / T)

    # accept move
    copy!(xs, new_xs);
    current_residuals = new_residuals;
    residual_norm = norm(current_residuals, 2);
    Is = copy(new_Is);

    # check for convergence
    residual_norm = norm(current_residuals);
  end
  
  # check to see if relaxation should occur
  if iter % stepout == 0
    write(w, "$(join(xs, ',')),$residual_norm\n");
  end
 
  @update_user(begin
    println("Residual surface exploration iteration: $_NUMEVALS_");
    println("    current residuals: ", current_residuals);
    println("    Is: ", [Is[1][1], Is[2][1], Is[3][1]]);
  end);
end

close(w);
