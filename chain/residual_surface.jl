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
    default = 1000
  "--outfile", "-O"
    help = "output file"
    default = "residual_surface.csv"
end

pa = parse_args(s);
const dk = pa["dk"];
const Ex = pa["Ex"];
const Ez = pa["Ez"];
const ωx = dk * Ex*Ex; 
const ωz = dk * Ez*Ez;
const ωxz = dk * Ex*Ez;
const γ = pa["gamma"];
const γx = NaN;
const γz = NaN;
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

const DEFAULT_INIT_X = ((pa["init-x"] != nothing) ?
                            eval(parse(pa["init-x"])) : 
                            initial_guess(γ));

_NUMEVALS_ = 0;
_ITERS_PER_PRINT_ = convert(Int, round(print_every_percent_compl * iters));

xs = DEFAULT_INIT_X;
initial_Is = integrate_sys(xs);
current_residuals = residuals(initial_Is);
residual_norm = norm(current_residuals, 2);
w = open(outfile, "a");

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
  end);
end

close(w);
