using HCubature;
using Optim;

struct Ρ
  C::Real;
  λ::Real;
  α::Real;
end;

(ρ::Ρ)(ω::Real, ϕ::Real, θ::Real) = ρ.C * exp(ω*(cos(θ))^2 + ρ.λ*cos(θ) + ρ.α*cos(ϕ)*sin(θ));

ωrange = -20.0:0.1:20.0;



for ω in ωrange
  ρ = Ρ(); 
end
