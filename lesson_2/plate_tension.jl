using Gridap
using LinearAlgebra

tcf(filename::String) = (@__DIR__) * "\\" * filename;

const λ = 100.0
const μ = 1.0

# const E = 2.1e11
# const ν = 0.3


# function σ(eps)
#     # λ = E*ν/(1+ν)/(1-2*ν)
#     # μ = E/2/(1+ν)
#     # return lambda_ * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)
#     return λ*tr(eps)*I(2) + 2*μ*eps
# end

σ(e) = λ*tr(e)*one(e) + 2*μ*e
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
a(u,v) = ∫( ε(v) ⊙ (σ∘ε(u)) )*dΩ
l(v) = 0; #∫(dot(f, v))*dΩ

op = AffineFEOperator(a,l,U,V)

x0 = zeros(Float64, num_free_dofs(V))
uh_lin = FEFunction(U,x0)
ls = LUSolver()
solver = LinearFESolver(ls)

uh_lin,_ = solve!(uh_lin,solver,op)

res_file = "results_new" |> tcf;


writevtk(Ω, res_file ,  cellfields=["uh"=>uh_lin,"sigma"=>σ∘ε(uh_lin)])

mat = get_matrix(op)
