include(joinpath("inc", "algorithm_help_msg.jl"));

using ArgParse;

s = ArgParseSettings();
@add_arg_table s begin
  "--Ex", "-x"
    help = "dimensionless electric field in x-direction"
    arg_type = Float64
    default = 0.0
  "--Ez", "-z"
    help = "dimensionless electric field in z-direction"
    arg_type = Float64
    default = 0.0
  "--dk", "-D"
    help = "kappa2 - kappa1"
    arg_type = Float64
    default = 1.0
  "--kT", "-k"
    help = "dimensionless temperature"
    arg_type = Float64
    default = 1.0
  "--gamma", "-g"
    help = "dimensionless stretch"
    arg_type = Float64
    default = 0.0
  "--num-monomers", "-N"
    help = "number of monomers"
    arg_type = Int
    default = 100
  "--algorithm", "-A"
    help = "algorithm (ISRES | PRAXIS | SBPLX | SLSQP | LBFGS | AUGLAG | Newton | Anneal). Run with --print-algo-help for more detailed information"
    default = "try Newton then SBPLX,Anneal,PRAXIS"
  "--print-algo-help"
    help = "Get detailed information regarding algorithms themselves, and how to run to solvers in stages"
    action = :store_true
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
  "--reltol-opt", "-R"
    help = "relative tolerance for optimization"
    arg_type = Float64
    default = 0.0
  "--abstol-opt", "-B"
    help = "absolute tolerance for optimization"
    arg_type = Float64
    default = 1e-8
  "--maxevals-opt", "-E"
    help = "maximum evaluations for numerical integration"
    arg_type = ASCIIString
    default = "try 100 then 10000,100000,5000"
  "--init-guess", "-X"
    help = "initial guess"
    arg_type = ASCIIString
  "--anneal-factor", "-m"
    help = "annealing factor, gamma i.e. T_new = T * (1 - gamma)"
    arg_type = Float64
    default = 0.005
  "--anneal-init-T", "-t"
    help = "initial temperature"
    arg_type = Float64
    default = 25.0
  "--anneal-final-T", "-T"
    help = "final temperature"
    arg_type = Float64
    default = 0.0
  "--anneal-relax-steps", "-S"
    help = "number of steps before annealing and step size adaptation"
    arg_type = Int
    default = 100
  "--anneal-init-delta", "-s"
    help = "initial step size for annealing"
    arg_type = Float64
    default = 0.5
  "--anneal-adapt-factor", "-P"
    help = "step size adapation factor (should be greater than 1.0)"
    arg_type = Float64
    default = 1.2
  "--anneal-acc-thrshld", "-a"
    help = "step size adapation acceptance threshold (should be in [0, 1))"
    arg_type = Float64
    default = 0.7
end

pa = parse_args(s);
if pa["print-algo-help"]
  println(_ALGO_HELP_MSG_);
  exit();
end
const dk = pa["dk"];
const Ex = pa["Ex"];
const Ez = pa["Ez"];
const ωx = dk * Ex*Ex; 
const ωz = dk * Ez*Ez;
const ωxz = dk * Ex*Ez;
const _u0_ = (dk > 0) ? -dk*(Ex^2 + Ez^2) : 0.0; 
const γ = pa["gamma"];
const γx = NaN;
const γz = NaN;
const kT = pa["kT"];
const N = pa["num-monomers"];
const reltol_opt = pa["reltol-opt"];
const abstol_opt = pa["abstol-opt"];

const verbose = pa["verbose"];
const print_every_percent_compl = pa["print-every-percent-compl"];
const RELTOL_INT = pa["reltol-int"];
const MAXEVALS_INT = pa["maxevals-int"];
const ANNEAL_FACTOR = pa["anneal-factor"];
const ANNEAL_INIT_T = pa["anneal-init-T"];
const ANNEAL_FINAL_T = pa["anneal-final-T"];
const ANNEAL_RELAX_STEPS = pa["anneal-relax-steps"];
const ANNEAL_INIT_DELTA = pa["anneal-init-delta"];
const ANNEAL_ADAPT_FACTOR = pa["anneal-adapt-factor"];
const ANNEAL_ACC_THRSHLD = pa["anneal-acc-thrshld"];

include(joinpath("inc", "solvers.jl"));
include(joinpath("inc", "free_energy.jl"));
include(joinpath("inc", "dipole.jl"));
include(joinpath("inc", "output.jl"));

_RHS_ = _RHS_GCOORD_;
ω(ϕ::Real, θ::Real) = _ω_gcoord_(ϕ, θ);

main_algos = parse_algostr(pa["algorithm"]);
main_max_iters = parse_max_iters(pa["maxevals-opt"]);
const DEFAULT_INIT_GUESS = (
  if pa["init-guess"] != nothing
    eval(parse(pa["init-guess"]));
  elseif ωx == 0 && ωz == 0
    initial_guess(γ);
  elseif abs(ωx) + abs(ωz) < 1 && ωx != 0 && ωz != 0 
    ((1-abs(ωx)+abs(ωz))*initial_guess(γ) + 
     (abs(ωx)+abs(ωz))*initial_guess((ωx, ωz), γ));
  else
    initial_guess((ωx, ωz), γ);
  #=else
    initial_guess(γ);=#
  end
);
(xs, rs, converged, iterations) = run_solver(main_algos, main_max_iters;
                                             x0s = DEFAULT_INIT_GUESS,
                                             reltol = reltol_opt,
                                             abstol = abstol_opt);
if verbose
  println((converged) ? "Solution converged." : "Solution did NOT converge.");
  println("Iterations: ", iterations);
end
residuals_output(rs);
standard_solution_output(xs);
