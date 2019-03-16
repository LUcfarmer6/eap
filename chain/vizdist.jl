using ArgParse;
using Plots;
using Cubature;

@everywhere ρ(ϕ::Real, θ::Real, c::Real, ω::Real, λ::Real, α::Real, ) = 
  c * exp(-ω*(cos(θ))^2 + λ*cos(θ) + α*cos(ϕ)*sin(θ));

function viz_file(fname::String, xls::Vector, yls::Vector, zls::Vector,
                  header::Bool=false; ext::String=".eps")
  headers = [];
  if header
    data, headers = readdlm(fname, ','; header=header);
  else
    data, headers = readdlm(fname, ','), [];
  end

  pieces = split(split(fname, "/")[end], ".");
  rootname = pieces[end-1];
  subpieces = split(rootname, "_");
  E0, dk = NaN, NaN;
  foundE0Flag, foundDkFlag = false, false;
  for subpiece in subpieces
    subpiece, startswith(subpiece, "E0"), startswith(subpiece, "dk")
    if startswith(subpiece, "E0")
      E0 = parse(subpiece[4:end]) / 100.0;
      foundE0Flag = true;
      if foundDkFlag; break; end
    elseif startswith(subpiece, "dk")
      dk = (subpiece == "dk-pos") ? 1.0 : -1.0;
      foundDkFlag = true;
      if foundE0Flag; break; end
    end
  end
  if !foundE0Flag || !foundDkFlag
    return false;
  end

  ω = E0*E0*dk;
  nrows = size(data, 1);
  θs = linspace(0, pi, 100);
  ϕs = linspace(0, 2*pi, 100);
  calls = [];
  for i=1:nrows
    push!(calls, @spawn begin
      @show gamma, A, λ, α, c = data[i, 1:5];

      plotdata2d = map(θ -> begin
        I = Cubature.pquadrature(x -> ρ(x, θ, c, ω, λ, α), 0, 2*pi);
        return I[1];
      end, ϕs);
      xs = map( (ts) -> ts[1]*cos(ts[2]), zip(plotdata2d, ϕs) );
      ys = map( (ts) -> ts[1]*sin(ts[2]), zip(plotdata2d, ϕs) );
      Plots.plot(xs, ys,
                 title="\$\\gamma = $(@sprintf("%.3lf", gamma))\$");
      gammastr = @sprintf("%05d", gamma*10000);
      Plots.xlims!(xls[1], xls[2]);
      Plots.ylims!(yls[1], yls[2]);
      Plots.savefig(rootname * "_2d_gamma-$gammastr"*ext);

      nθ, nϕ = length(θs), length(ϕs);
      x3d = zeros(nθ*nϕ);
      y3d = zeros(nθ*nϕ);
      z3d = zeros(nθ*nϕ);
      for i=1:nθ, j=1:nϕ
        r3d = ρ(ϕs[j], θs[i], c, ω, λ, α);
        z3d[(i-1)*nϕ+j] = r3d * cos(θs[i]);
        x3d[(i-1)*nϕ+j] = r3d * cos(ϕs[j])*sin(θs[i]);
        y3d[(i-1)*nϕ+j] = r3d * sin(ϕs[j])*sin(θs[i]);
      end
      Plots.surface(x3d, y3d, z3d, 
                    title="\$\\gamma = $(@sprintf("%.3lf", gamma))\$");
      Plots.xlims!(xls[1], xls[2]);
      Plots.ylims!(yls[1], yls[2]);
      Plots.zlims!(zls[1], zls[2]);
      Plots.savefig(rootname * "_3d_gamma-$gammastr"*ext);
    end);
  end

  for call in calls; fetch(call); end

  return true;
end

s = ArgParseSettings();
@add_arg_table s begin
  "--pattern", "-P"
    help = "is the filename argument a regex pattern"
    action = :store_true
  "--dir", "-D"
    help = "directory that the file(s) is stored in"
    default = "."
  "--verbose", "-v"
    help = "print diagnostic information"
    action = :store_true
  "--header", "-H"
    help = "the delimited file has a header"
    action = :store_true
  "--figure-ext", "-E"
    help = "file extension for figures"
    default = ".eps"
  "--backend", "-B"
    help = "plotting backend (pyplot|gr|plotly|plotlyjs)"
    default = "pyplot"
  "--fontscale", "-S"
    help = "font scale for figures"
    default = 1.0
  "--xlims", "-X"
    help = "y limits"
    default = "[-100.0; 100.0]"
  "--ylims", "-Y"
    help = "y limits"
    default = "[-100.0; 100.0]"
  "--zlims", "-Z"
    help = "z limits"
    default = "[-100.0; 100.0]"
  "filename"
    help = "Data filename (or regex pattern)"
    required = true
end
pa = parse_args(s);

const patternflag = pa["pattern"];
const dir = pa["dir"];
const fname = pa["filename"];
const verbose = pa["verbose"];
const header = pa["header"];
const figext = pa["figure-ext"];
const xlims = eval(parse(pa["xlims"]));
const ylims = eval(parse(pa["ylims"]));
const zlims = eval(parse(pa["zlims"]));
for i=1:nprocs()
  remotecall_fetch(() -> begin
    if pa["backend"] == "pyplot"
      Plots.pyplot();
    elseif pa["backend"] == "gr"
      Plots.gr();
    elseif pa["backend"] == "plotly"
      Plots.plotly();
    elseif pa["backend"] == "plotlyjs"
      Plots.plotlyjs();
    end
    Plots.scalefontsizes(pa["fontscale"]);
  end, i);
end

if patternflag
  files = filter(Regex(fname), readdir(dir));
  if verbose; println("matched ", length(files), " files to viz"); end;
  for file in files
    if verbose; println("vizzing file: ", file); end;
    if !viz_file(joinpath(dir, file), xlims, ylims, zlims, header; ext=figext);
      warn("Could not visualize the data in file: $file");
    end
  end
else
  if verbose; println("vizzing file: ", fname); end;
  if !viz_file(joinpath(dir, fname), xlims, ylims, zlims, header; ext=figext);
    warn("Could not visualize the data in file: $fname");
  end
end
