using Gridap
using Gridap.Arrays
using Gridap.Fields
using LinearAlgebra

tcf(filename::String) = (@__DIR__) * "\\" * filename;


const E = 70.0e9
const ν = 0.33
const λ = (E*ν)/((1+ν)*(1-2*ν))
const μ = E/(2*(1+ν))
const ρ = 7.9e3
σ(ε) = λ*tr(ε)*one(ε) + 2*μ*ε

domain = (0,0.10,0,0.01,0,0.01)
partition = (100,10,10)
model = CartesianDiscreteModel(domain,partition)


labels = get_face_labeling(model)
add_tag_from_tags!(labels,"diri_0",[25])
# add_tag_from_tags!(labels,"diri_1",[2,4,8])

degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

model_file = "model" |> tcf;
writevtk(model, model_file)


reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},1)
V = TestFESpace(model,reffe,conformity=:H1,dirichlet_tags = ["tag_25"])

g0 = VectorValue(0.0,  0.0, 0.0)
# g1 = VectorValue(0.0, 0.0)

U = TrialFESpace(V, [g0])

# weak form 

# f = VectorValue(0.0, - gravity * ρ)
bfK(u,v) = ∫( ε(v) ⊙ (σ∘ε(u)) )*dΩ
L(v) = 0; #∫(dot(f, v))*dΩ


# K*x = f

# Mẍ + Kx = 0 

op_K = AffineFEOperator(bfK, L, U, V)
mK = get_matrix(op_K)

bfM(u,v) = ρ*∫(v⋅u)dΩ
op_M = AffineFEOperator(bfM, L, U, V)

mM = get_matrix(op_M)

using Arpack

res = eigs(mM, mK; nev = 10, which=:LM)
res[2]
vec = real.(res[2])


uh_lin = FEFunction(U,vec[:,3])


res_file = "eigen_freq" |> tcf;

writevtk(Ω, res_file ,  cellfields=["uh"=>uh_lin])
