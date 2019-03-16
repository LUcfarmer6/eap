@everywhere Exs = [0.0; 0.05; 0.25; 0.5; sqrt(0.75); 1.0; 3.0; 5.0];
@everywhere Ezs = copy(Exs);
@everywhere dks = [-1; 1];

if !isdir(joinpath("data", "small-omega_G-coord"))
  mkdir(joinpath("data", "small-omega_G-coord"))
end

@show pmap(Ex -> begin
  for Ez in Ezs, dk in dks
    if abs(Ex) < 1e-12 && abs(Ez) < 1e-12 && dk == -1 # we require ωx and ωz not be opposite signs
      continue;
    end

    dk_str = (dk == 1) ? "pos" : "neg";

    println("Creating graph for Ex = $Ex, Ez = $Ez, dk = $dk");
    outfile = joinpath("data", "small-omega_G-coord", @sprintf("opt_Ex-%05d_Ez-%05d_dk-%s.csv", round(Ex*100), round(Ez*100), dk_str));

    if isfile(outfile)
      println(outfile, " already exists. Skipping...");
      continue;
    end

    run(`julia create_graph.jl --Ex $Ex --Ez $Ez --dk $dk --num-data-points=25 -O $outfile --verbose`);
  end
  return 0;
end, Exs);
