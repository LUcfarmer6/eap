using NLopt;
using ArgParse;

include(joinpath("inc", "kinematics.jl"));

@everywhere function cost(x::Vector, grad::Vector, data::Matrix)
  if length(grad) > 0
    error("Gradient not implemented");
  end
end

@everywhere function fit_files(fnames::Vector{String}, header::Bool=false)
  headers = [];
  if header
    data, headers = readdlm(fname[1], ','; header=header);
  else
    data, headers = readdlm(fname[1], ','), [];
  end

  for i=2:length(fnames)
    if header
      idata, _headers = readdlm(fname[1], ','; header=header);
      if size(data, 2) == size(idata, 2)
        data = vcat(data, idata);
      else
        return false;
      end
    else
      idata, _headers = readdlm(fname[1], ','), [];
      if size(data, 2) == size(idata, 2)
        data = vcat(data, idata);
      else
        return false;
      end
    end
  end
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

if patternflag
  files = filter(Regex(fname), readdir(dir));
  if verbose; println("matched ", length(files), " files to fit"); end;
  if verbose; println("fitting file: ", files); end;
  results, success = fit_files(map(file -> joinpath(dir, file), files), header);
  if !success;
    warn("Could not fit the data in files: $files");
  end
else
  if verbose; println("fitting file: ", fname); end;
  results, success = fit_files([joinpath(dir, fname)], header);
  if !success;
    warn("Could not fit the data in file: $fname");
  end
end
