using NLopt;
using ArgParse;



@everywhere function fit_file(fname::String, header::Bool=false)
  headers = [];
  if header
    data, headers = readdlm(fname, ','; header=header);
  else
    data, headers = readdlm(fname, ','), [];
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
  if verbose; println("matched ", length(files), " files to sanitized"); end;
  for file in files
    if verbose; println("fitting file: ", file); end;
    if !fit_file(joinpath(dir, file), header);
      warn("Could not fit the data in file: $file");
    end
  end
else
  if verbose; println("fitting file: ", fname); end;
  if !viz_file(joinpath(dir, fname), header);
    warn("Could not fit the data in file: $fname");
  end
end
