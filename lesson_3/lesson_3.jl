tcf(filename::String) = (@__DIR__) * "\\" * filename;
using Gmsh

# ------------------------------------------------------------------------------
#
#  Gmsh Julia tutorial 20
#
#  STEP import and manipulation, geometry partitioning
#
# ------------------------------------------------------------------------------

# The OpenCASCADE CAD kernel allows to import STEP files and to modify them. In
# this tutorial we will load a STEP geometry and partition it into slices.


gmsh.initialize()
gmsh.option.setNumber("Mesh.MeshSizeMax", 3.0) #Максимальный размер элементов внутри всей модели - 3мм
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 50)#колличество элементов на 2pi радиан
gmsh.option.setNumber("Mesh.ElementOrder", 2)#колличество элементов на 2pi радиан

gmsh.model.add("model")

# Load a STEP file (using `importShapes' instead of `merge' allows to directly
# retrieve the tags of the highest dimensional imported entities):
path = "lug.stp" |> tcf
v = gmsh.model.occ.importShapes(path)

# If we had specified
#
# gmsh.option.setString('Geometry.OCCTargetUnit', 'M')
#
# before merging the STEP file, OpenCASCADE would have converted the units to
# meters (instead of the default, which is millimeters).

# Get the bounding box of the volume:
#xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(
 #   v[1][1], v[1][2])




gmsh.model.occ.synchronize()
for plane in gmsh.model.getEntities(2)
    println(gmsh.model.getColor(plane...))
    gmsh.model.addPhysicalGroup(2, [plane[2]], -1, "$(plane[2])");

end

#colors4 = [gmsh.model.getColor(plane...) for plane in gmsh.model.getEntities(2)]
#Поиск граней с заданным цветом
faces = Dict()
for plane in gmsh.model.getEntities(2)
    faces[gmsh.model.getColor(plane...)] = plane[2];
    # gmsh.model.addPhysicalGroup(2, [plane[2]], -1, "$(plane[2])");
end
faces[(100,0,0,255)] #<- returns face_id

d3 = gmsh.model.getEntities(3)
gmsh.model.addPhysicalGroup(3, [d3[1][2]], -1, "domain")

d2 = gmsh.model.getEntities(2)
gmsh.model.addPhysicalGroup(2, [faces[(0,100,0,255)]], -1, "loaded")
gmsh.model.addPhysicalGroup(2, [faces[(100,0,0,255)]], -1, "fixed")

gmsh.model.occ.synchronize()

# xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(
#     d2[13][1], d2[13][2])
# Δ = 0.1
# gmsh.model.mesh.field.add("Box", 17)
# gmsh.model.mesh.field.setNumber(17, "VIn", 0.5)
# gmsh.model.mesh.field.setNumber(17, "VOut", 3.0)
# gmsh.model.mesh.field.setNumber(17, "XMin", xmin-Δ)
# gmsh.model.mesh.field.setNumber(17, "XMax", xmax+Δ)
# gmsh.model.mesh.field.setNumber(17, "YMin", ymin-Δ)
# gmsh.model.mesh.field.setNumber(17, "YMax", ymax+Δ)
# gmsh.model.mesh.field.setNumber(17, "ZMin", zmin-Δ)
# gmsh.model.mesh.field.setNumber(17, "ZMax", zmax+Δ)
# gmsh.model.occ.synchronize()
# gmsh.option.setNumber("Mesh.Algorithm", 5)
# gmsh.option.setNumber("Mesh.Algorithm3D", 5)



gmsh.model.mesh.generate(3)
msh_file = "lug.msh" |> tcf;
gmsh.write(msh_file)

# gmsh.fltk.run()

gmsh.finalize()


using Gridap
using GridapGmsh

model = GmshDiscreteModel(msh_file);

vtk_file = "model" |> tcf;
writevtk(model, vtk_file);


##





model = GmshDiscreteModel(msh_file);

reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},1) #Как и в предыдущем уроке, мы строим непрерывную интерполяцию Лагранжа первого порядка...
#...Векторно-значная интерполяция выбирается с помощью "VectorValue".

V = TestFESpace(model,reffe,conformity=:H1,dirichlet_tags = ["fixed"])

g(x) = VectorValue(0.0,  0.0, 0.0) # Граничное условие на левой грани – перемещение 0.

U = TrialFESpace(V, [g]) # Создаем пространство тестовых функций.
EL_ORDER = 1
degree = EL_ORDER*2
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

neumanntags = ["loaded"]
Γ = BoundaryTriangulation(model,tags=neumanntags)
dΓ = Measure(Γ,degree)

f(x) = VectorValue(0.0, 0.0, 0.10);

const E = 2.1e5
const ν = 0.3
const λ = (E*ν)/((1+ν)*(1-2*ν))
const μ = E/(2*(1+ν))

σ(ε) = λ*tr(ε) + 2*μ*ε

a(u,v) = ∫(ε(v) ⊙ (σ∘ε(u)) )*dΩ
b(v) = ∫( dot(v, f))*dΓ 


op = AffineFEOperator(a,b,U,V)
x0 = zeros(Float64, num_free_dofs(V))
uh_lin = FEFunction(U,x0)
ls = BackslashSolver()
solver = LinearFESolver(ls)

uh_lin, _ = solve!(uh_lin,solver,op)

res_file = "results_new" |> tcf;

function mises(s)
    return 0.5 * sqrt((s[1,1] - s[2,2])^2 + (s[2,2] - s[3,3])^2 + (s[3,3] - s[1,1])^2 + 6 * (s[2,3]^2 + s[3,1]^2 + s[1,2]^2))
end

writevtk(Ω, res_file ,  cellfields=["uh"=>uh_lin,"sigma"=>σ∘ε(uh_lin)])