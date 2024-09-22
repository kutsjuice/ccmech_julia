# # Построение геометрии


# некоторые вспомагательные функции, которые понадобятся в процессе

findin(a, b) = a[findall(in(b), a)]
tcf(filename::String) = (@__DIR__) * "\\" * filename;


# Подключаем gmsh и инциализируем его, после чего моздаем модель и задаем базовый размер репрезентативной ячейки:
using Gmsh

gmsh.initialize()

factory = gmsh.model.occ
gmsh.model.add("model")

cell_size = 1;

# Исследуемая геометрия чем-то напомнимает сыр и представляет собой прямоугольник с крыглыми отверстиями. 
# Для этого создадим прямоугольник размером `cell_size`, затем некоторое количество окружностей, которые в дальнейшем будут вырезаны из прямоугольной области. 
# Следует заметить, что булевы операции возможны лишь для 2D или 3D объектов. 
# Тут необходимо учесть некоторую особонность построения геометрии в **gmsh** - при создании прямоугольника, будет созданы необходимые отрезки, замкнутый контур и участок поверхности, ограниченный данным контуром. 
# В случае же окружностей будут созданы только лишь замкнутые дуги, т.е. 1D объекты. 
# Поэтому прежде чем произвести операцию вырезки, необходимо создать контур для каждой окружности (состоящий из одной окружности) и "натянуть" на созданные контура поверхности при помощи функции `addPlaneSurface`.

rect = factory.addRectangle(-0.5 * cell_size, -0.5 * cell_size, 0, cell_size, cell_size);

#src # factory.getEntities(2)
#src # full_domain = factory.addPlaneSurface([rect])

holes = [];

for i in 1:20
    xy = rand(2) .- 0.5;;
    r = 0.05;
    circle = factory.addCircle(xy[1], xy[2], 0, r);
    inner_loop = factory.addCurveLoop([circle]);
    hole = factory.addPlaneSurface([inner_loop]);
    push!(holes, (2, hole));
end

# Теперь можно удалить созданные кругляши из прямоугольника, воспользовавшись функцией `cut`

diff = factory.cut([(2, rect)], holes, -1, true, true)

factory.synchronize()

# В некоторых случаях может оставаться избыточная геометрия, которая будет мешать создать модель. 
# Поэтому удалим все 1D объекты, которые не входят в границу нашей ячейки

entities = factory.getEntities(2)

surface_boundary = gmsh.model.getBoundary(entities)

edges_for_leave = [abs(tag[2]) for tag in surface_boundary]

all_edges = [abs(tag[2]) for tag in factory.getEntities(1)]

edges_for_remove = []
for edge in all_edges
    if  !(edge in edges_for_leave)
        push!(edges_for_remove, edge)
    end
end

dimtag_for_remove = [(1, tag) for tag in edges_for_remove]
factory.remove(dimtag_for_remove)

factory.synchronize()

# Далее необходимо создаnm физические группы для созданной геометрией. Группа для всей геометрии:
for dim in 0:2
    entities = gmsh.model.occ.getEntities(dim)
    tags = [entity[2] for entity in entities]
    gmsh.model.addPhysicalGroup(dim, tags, -1, "domain")
end

# Далее создаем физические группы для границ.
# Границы можно выделить при помощи инструмента `boundingBox`, что напоминает выделение мышкой на экране прямоугольной области. 
# Для начала создадим 4 области, в которые входят боковые границы ячейки: 

eps = cell_size * 1e-4;
left_boundary_bbox = [
    (- cell_size/2 - eps),
    (- cell_size/2 - eps),
    (- eps),
    (- cell_size/2 + eps),
    (+ cell_size/2 + eps),
    (+ eps)
]
right_boundary_bbox = [
    (+ cell_size/2 - eps),
    (- cell_size/2 - eps),
    (- eps),
    (+ cell_size/2 + eps),
    (+ cell_size/2 + eps),
    (+ eps)
]
bottom_boundary_bbox = [
    (- cell_size/2 - eps),
    (- cell_size/2 - eps),
    (- eps),
    (+ cell_size/2 + eps),
    (- cell_size/2 + eps),
    (+ eps)
]
top_boundary_bbox = [
    (- cell_size/2 - eps),
    (+ cell_size/2 - eps),
    (- eps),
    (+ cell_size/2 + eps),
    (+ cell_size/2 + eps),
    (+ eps)
]

bounds = [left_boundary_bbox, bottom_boundary_bbox, right_boundary_bbox, top_boundary_bbox]


#src # # Далее сгруппируем эти области в 2 словаря: первый используется, для выделения боковых граней
#src # edge_bboxes = Dict(
#src #     "left" => left_boundary_bbox, 
#src #     "right" => right_boundary_bbox,
#src #     "bottom" => bottom_boundary_bbox,
#src #     "top" => top_boundary_bbox
#src # )

#src # node_bboxes = Dict(
#src #     "top_right" => (top_boundary_bbox, right_boundary_bbox),
#src #     "top_left" => (top_boundary_bbox, left_boundary_bbox), 
#src #     "bottom_right" => (bottom_boundary_bbox, right_boundary_bbox), 
#src #     "bottom_left" => (bottom_boundary_bbox, left_boundary_bbox)
#src # )

#src # for (boundary_name, bbox) in edge_bboxes
#src #     dim = 1
#src #     entities = gmsh.model.occ.getEntitiesInBoundingBox(bbox..., dim)   
#src #     tags = [entity[2] for entity in entities]
#src #     gmsh.model.addPhysicalGroup(dim, tags, -1, boundary_name)
#src # end

#src # for (boundary_name, bbox) in node_bboxes
#src #     dim = 0
#src #     entities1 = gmsh.model.occ.getEntitiesInBoundingBox(bbox[1]..., dim)   
#src #     entities2 = gmsh.model.occ.getEntitiesInBoundingBox(bbox[2]..., dim)   
#src #     entities = findin(entities1, entities2)
#src #     tags = [entity[2] for entity in entities]
#src #     gmsh.model.addPhysicalGroup(dim, tags, -1, boundary_name)
#src # end

# Далее создадим две физические группы: с узлами (1D объекты) и с ребрами (2D объекты)

boundary_name = "boundary";
for dim in 0:1
    entities = []
    for bbox in bounds
        push!(entities, factory.getEntitiesInBoundingBox(bbox..., dim)...)
    end
    tags = [entity[2] for entity in entities]
    gmsh.model.addPhysicalGroup(dim, tags, -1, boundary_name)
end
gmsh.model.occ.synchronize()

# Теперь осталось задать параметры сетки и запустить генерацию сетки, после чего созданная геометрия может быть сохранена в файл


gmsh.option.setNumber("Mesh.MeshSizeMin", 0.01)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.05)
#src # gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.Algorithm", 8)
#src # gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3) # or 3
#src # gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 20)

gmsh.model.mesh.generate(2)
#src # gmsh.model.mesh.recombine()
#src # gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
#src # gmsh.model.mesh.refine()
#src
#src # if !("-nopopup" in ARGS)
#src #      gmsh.fltk.run()
#src # end

msh_file = "porous_media.msh" |> tcf;

gmsh.write(msh_file)

gmsh.finalize()
#src
# Теперь необходимо загрузить нашу геометрию при помощи Gridap и 

using Gridap
using GridapGmsh

model = GmshDiscreteModel(msh_file);

vtk_file = "porous_model" |> tcf;
#src # writevtk(model, vtk_file);


# Зададим основные механические параметры для материала и определим закон для напряжений (закон Гука):
const E = 70.0e9
const ν = 0.33
const λ = (E*ν)/((1+ν)*(1-2*ν))
const μ = E/(2*(1+ν))
const ρ = 2.8e3
σ(ε) = λ*tr(ε)*one(ε) + 2*μ*ε

# Задаем пространство тестовых функций. В качестве границы Дирихле укажем границу с тегом `boundary` которую мы создали в **gmsh**
degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

reffe = ReferenceFE(lagrangian, VectorValue{2,Float64}, 1)
V = TestFESpace(model, reffe, 
    conformity=:H1,
    dirichlet_tags = [boundary_name]
)

# Теперь нам необходимо решить 3 краевых задач, каждой из которых соответствуют свои граничные условия.

g1(x) = VectorValue(x[1],  0.0)
g2(x) = VectorValue(0.0,  x[2])
g3(x) = VectorValue(x[2]/(2cell_size),  x[1]/(2cell_size))
    
U1 = TrialFESpace(V, [g1])
U2 = TrialFESpace(V, [g2])
U3 = TrialFESpace(V, [g3])

U = [U1, U2, U3]

# Слабая постановка для всех трех задач будет одинаковая

A(u,v) = ∫( ε(v) ⊙ (σ∘ε(u)) )*dΩ
L(v) = 0; 

# Теперь решаем все три задачи и сохраняем результат в файл
res_file = "results_new" |> tcf;
result = []
for i in 1:3 
    op = AffineFEOperator(A, L, U[i], V)

    x0 = zeros(Float64, num_free_dofs(V))
    uh_lin = FEFunction(U[i], x0)
    ls = LUSolver()
    solver = LinearFESolver(ls)
    uh_lin ,_ = solve!(uh_lin,solver,op)
    push!(result, uh_lin);

end
writevtk(Ω, res_file ,  cellfields=[
    "x_tension_displ" => result[1],
    "y_tension_displ" => result[2],
    "xy_shear_displ" => result[3],
    "x_tension_stress"=> σ∘ε(result[1]),
    "y_tension_stress"=> σ∘ε(result[2]),
    "xy_shear_stress"=> σ∘ε(result[3]),
    ]
)