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

# Исследуемая геометрия чем-то напомнимает сыр и представляет собой прямоугольник с крыглыми отверстиями. Для этого создадим прямоугольник размером `cell_size`, 
# затем некоторое количество окружностей, которые в дальнейшем будут вырезаны из прямоугольной области. Следует заметить, что булевы операции возможны лишь для 2D или 3D объектов,
# однако прямоугольник и окружность это по своей сути 1D объект. Поэтому прежде чем произвести операцию вырезки, необходимо "натянуть" на созданные контура прямоугольника и окружности 
# поверхности при помощи функции `addPlaneSurface`.

rect = factory.addRectangle(-0.5 * cell_size, -0.5 * cell_size, 0, cell_size, cell_size);

factory.getEntities(2)
# full_domain = factory.addPlaneSurface([rect])

holes = [];
for i in 1:20
    xy = rand(2) .- 0.5;;
    r = 0.05;
    circle = factory.addCircle(xy[1], xy[2], 0, r);
    inner_loop = factory.addCurveLoop([circle]);
    hole = factory.addPlaneSurface([inner_loop]);
    push!(holes, (2, hole));
end


diff = factory.cut([rect], holes, -1)

factory.synchronize()

#№ После соз

#№ entities = gmsh.model.occ.getEntities(2)
#№ gmsh.model.occ.remove(entities[1:1])


#№ factory.synchronize()


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

dim = 2
entities = gmsh.model.occ.getEntities(dim)
tags = [entity[2] for entity in entities]
gmsh.model.addPhysicalGroup(dim, tags, -1, "domain")

# Создаем физические группы для границ

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

edge_bboxes = Dict(
    "left" => left_boundary_bbox, 
    "right" => right_boundary_bbox,
    "bottom" => bottom_boundary_bbox,
    "top" => top_boundary_bbox
)

node_bboxes = Dict(
    "top_right" => (top_boundary_bbox, right_boundary_bbox),
    "top_left" => (top_boundary_bbox, left_boundary_bbox), 
    "bottom_right" => (bottom_boundary_bbox, right_boundary_bbox), 
    "bottom_left" => (bottom_boundary_bbox, left_boundary_bbox)
)

for (boundary_name, bbox) in edge_bboxes
    dim = 1
    entities = gmsh.model.occ.getEntitiesInBoundingBox(bbox..., dim)   
    tags = [entity[2] for entity in entities]
    gmsh.model.addPhysicalGroup(dim, tags, -1, boundary_name)
end

for (boundary_name, bbox) in node_bboxes
    dim = 0
    entities1 = gmsh.model.occ.getEntitiesInBoundingBox(bbox[1]..., dim)   
    entities2 = gmsh.model.occ.getEntitiesInBoundingBox(bbox[2]..., dim)   
    entities = findin(entities1, entities2)
    tags = [entity[2] for entity in entities]
    gmsh.model.addPhysicalGroup(dim, tags, -1, boundary_name)
end


gmsh.model.occ.synchronize()

gmsh.option.setNumber("Mesh.MeshSizeMin", 0.01)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.05)
## gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.Algorithm", 8)
## gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3) # or 3
## gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 10)

gmsh.model.mesh.generate(2)
## gmsh.model.mesh.recombine()
## gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
## gmsh.model.mesh.refine()
## gmsh.model.mesh.refine()



## if !("-nopopup" in ARGS)
##     gmsh.fltk.run()
## end

msh_file = "porous_media.msh" |> tcf;
print(msh_file)
gmsh.write(msh_file)

gmsh.finalize()

#
using Gridap
using GridapGmsh

model = GmshDiscreteModel(msh_file);

vtk_file = "porous_model" |> tcf;
writevtk(model, vtk_file);
