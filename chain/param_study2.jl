#E0s = [0.0; 0.5; 1.0; 3.0; 5.0; 6.0];
@everywhere E0s = [0.0; 0.5; 1.0; 3.0; 5.0; 6.0];
@everywhere dks = [-1; 1];
@everywhere fracs = Tuple{Int,Int}[(1, 9999); (1, 12); (1, 6); (1, 4); (1, 3); (5, 12); (1, 2)];

if !isdir(joinpath("data", "small-lambda_E-coord"))
  mkdir(joinpath("data", "small-lambda_E-coord"))
end

@show pmap(frac -> begin
  for dk in dks, E0 in E0s
    if abs(E0) < 1e-12 && (dk == -1 || frac != (1, 9999))
      continue;
    end

    dk_str = (dk == 1) ? "pos" : "neg";
    theta_str = ((frac[2] == 9999) ? "0" : 
                 ((frac[1] == 1) ? "pid$(frac[2])" : "$(frac[1])pid$(frac[2])"));
    theta = (frac[2] == 9999) ? 0 : frac[1]*pi / frac[2];

    println("Creating graph for E0 = $E0, dk = $dk, theta = $theta_str");
    outfile = joinpath("data", "small-lambda_E-coord", @sprintf("opt_E0-%05d_dk-%s_theta-%s.csv", round(E0*100), dk_str, theta_str));

    if isfile(outfile)
      println(outfile, " already exists. Skipping...");
      continue;
    end

    run(`julia create_graph_ecoord.jl --E0 $E0 --dk $dk --theta $theta --num-data-points=25 -O $outfile --verbose`);
  end
  return 0;
end, fracs);
