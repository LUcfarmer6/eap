using DelimitedFiles;

mathematica_to_julia(x::Number) = x;

function mathematica_to_julia(numstr::AbstractString)
  if numstr == "Indeterminate" || (findnext(numstr, "`", 1) != nothing);
    return NaN;
  end

  newstr = replace(numstr, "*I" => "im");
  newstr = replace(newstr, "*^" => "e");
  newstr = replace(newstr, "ComplexInfinity" => "Inf");
  newstr = replace(newstr, "Infinity" => "Inf");
  newstr = replace(newstr, ".+" => "+");
  newstr = replace(newstr, ".-" => "-");

  pieces = split(newstr, "e");
  if length(pieces) > 1
    exponent = parse(Float64, split(pieces[2],"*I")[1]);
    if exponent > 100
      return NaN;
    end
  end

  x = parse(Complex{Float64}, newstr);
  return real(x);
end

function sanitize_file(fname; sanitize_function::Function = mathematica_to_julia,
                       header=false)
  unsanitized_data = zeros(1, 1);
  headers = [];
  if header
    unsanitized_data, headers = readdlm(fname, ','; header=header);
  else
    unsanitized_data, headers = readdlm(fname, ','), [];
  end

  sanitized_data = map(x -> sanitize_function(x), unsanitized_data);

  if header
    sanitized_data = vcat(headers, sanitized_data);
  end

  global overwriteflag;
  if overwriteflag
    writedlm(fname, sanitized_data, ',');
  else
    writedlm("sanitized_" * fname, sanitized_data, ',');
  end
end

using ArgParse;

s = ArgParseSettings();
@add_arg_table s begin
  "--do-not-overwrite", "-W"
    help = "do not overwrite existing files"
    action = :store_true
  "--pattern", "-P"
    help = "is the filename argument a regex pattern"
    action = :store_true
  "--dir", "-D"
    help = "directory that the file(s) is stored in"
    default = "."
  "--bound", "-b"
    help = "artifical bound on maximum abs value allowed"
    arg_type = Float64
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

const overwriteflag = !(pa["do-not-overwrite"]);
const patternflag = pa["pattern"];
const dir = pa["dir"];
const fname = pa["filename"];
const bound = pa["bound"];
const verbose = pa["verbose"];
const header = pa["header"];

if patternflag
  regex = Regex(fname);
  files = filter(f -> (match(regex, f) != nothing), readdir(dir));
  if verbose; println("matched ", length(files), " files to sanitized"); end;
  for file in files
    if verbose; println("sanitizing file: ", file); end;
    sanitize_file(joinpath(dir, file); 
      sanitize_function=(bound != nothing) ? (x -> begin
        new_x = mathematica_to_julia(x);
        return (abs(new_x) > abs(bound)) ? NaN : new_x;
    end) : mathematica_to_julia, header=header);
  end
else
  if verbose; println("sanitizing file: ", fname); end;
  sanitize_file(joinpath(dir, fname); 
      sanitize_function=(bound != nothing) ? (x -> begin
        new_x = mathematica_to_julia(x);
        return (abs(new_x) > abs(bound)) ? NaN : new_x;
    end) : mathematica_to_julia, header=header);
end
