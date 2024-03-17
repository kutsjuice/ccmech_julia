using Gridap
using Gridap.Arrays
using Gridap.Fields
using LinearAlgebra

tcf(filename::String) = (@__DIR__) * "\\" * filename;

const E = 2.1e11
const ν = 0.3

const λ = E*ν/(1+ν)/(1-2*ν)
const μ = E/2/(1+ν)


# Model
flag = 1;
k1 = ν/(1-ν);
k2 = 0.5*(1-2ν)/(1-ν);
k3 = (1-ν)/2;

plain_stress = true;
mE = Matrix{Float64}(undef, 3, 3);
if plain_stress
    mE[:,:] = [
        1  k1  0;
        k1  1  0;
        0   0 k2
    ]*E*(1-ν)/((1+ν)*(1-2ν));   
else
    mE[:,:] = [
        1  ν  0;
        ν  1  0;
        0  0 k3;
    ]*E/(1-ν^2);
end

function voigtstrain(∇u)
    ε_x = ∇u[1,1];
    ε_y = ∇u[2,2];
    γ_xy = ∇u[1,2] + ∇u[2,1];

    return VectorValue(ε_x, ε_y, γ_xy);
end

function voigtstress(ε)
    σ_x = mE[1,1] * ε[1] + mE[1,2] * ε[2] + mE[1,3] * ε[3];
    σ_y = mE[2,1] * ε[1] + mE[2,2] * ε[2] + mE[2,3] * ε[3];
    τ_xy = mE[3,1] * ε[1] + mE[3,2] * ε[2] + mE[3,3] * ε[3];
    return VectorValue(σ_x, σ_y, τ_xy);
end

ϵ(u) = Operation(voigtstrain)(∇(u))
σ(u) = Operation(voigtstress)(ϵ(u))

# Model

domain = (0,1,0,1)
partition = (40,40)
model = CartesianDiscreteModel(domain,partition)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"diri_0",[1,3,7])
add_tag_from_tags!(labels,"diri_1",[2,4,8])

degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

reffe = ReferenceFE(lagrangian,VectorValue{2,Float64},1)
V = TestFESpace(model,reffe,conformity=:H1,dirichlet_tags = ["diri_0", "diri_1"])

g0 = VectorValue(0.0,  0.0)
g1 = VectorValue(0.75, 0.0)

U = TrialFESpace(V, [g0,g1])

# weak form 

# f = VectorValue(0.0, - gravity * ρ)
A(u,v) = ∫( ϵ(v)⋅σ(u))*dΩ
L(v) = 0; #∫(dot(f, v))*dΩ

op = AffineFEOperator(A, L, U,V)

x0 = zeros(Float64, num_free_dofs(V))
uh_lin = FEFunction(U,x0)
ls = LUSolver()
solver = LinearFESolver(ls)

uh_lin,_ = solve!(uh_lin,solver,op)

res_file = "results_new" |> tcf;


writevtk(Ω, res_file ,  cellfields=["uh"=>uh_lin,"sigma"=>σ(uh_lin)])
