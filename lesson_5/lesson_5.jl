
using Gmsh

gmsh.initialize()
tcf(filename::String) = (@__DIR__) * "\\" * filename;


lc = 5e-2

factory = gmsh.model.occ
gmsh.model.add("model")


# p1 = factory.addPoint(0,0,0, lc);
# p2 = factory.addPoint(1,0,0, lc);
# p3 = factory.addPoint(1,1,0, lc);
# p4 = factory.addPoint(0,1,0, lc);

# l1 = factory.addLine(p1, p2);
# l2 = factory.addLine(p2, p3);
# l3 = factory.addLine(p3, p4);
# l4 = factory.addLine(p4, p1);


rect = factory.addRectangle(-0.5, -0.5, 0, 1, 1);
# outter_loop = factory.addCurveLoop([c1]);
full_domain = factory.addPlaneSurface([rect])

holes = [];
for i in 1:20
    xy = rand(2) .- 0.5;;
    r = 0.05;
    circle = factory.addCircle(xy[1], xy[2], 0, r);
    inner_loop = factory.addCurveLoop([circle]);
    hole = factory.addPlaneSurface([inner_loop]);
    push!(holes, (2, hole));
end


diff = factory.cut([(2,full_domain)], holes, -1, false, true)

# println(diff)
# frag = factory.fragment(diff[0], [(2, Surface2)])    
# print(frag1)

# filled_entities = diff[2][end]
# main_entities = copy(diff[1])
# for entity in filled_entities
#   if entity in main_entities
    # main_entities.remove(entity)
#   end
# end
# print("  ", main_entities, "\n")
factory.synchronize()

entities = gmsh.model.occ.getEntities(2)
gmsh.model.occ.remove(entities[1:2])

factory.synchronize()

for dim in 2:2
    entities = gmsh.model.occ.getEntities(dim)
    tags = [entity[2] for entity in entities]
    gmsh.model.addPhysicalGroup(dim, tags, -1, "domain")
end

gmsh.model.occ.synchronize()

gmsh.option.setNumber("Mesh.MeshSizeMin", 0.01)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.05)
# gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.Algorithm", 8)
# gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3) # or 3
# gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 30)

gmsh.model.mesh.generate(2)
# gmsh.model.mesh.recombine()
# gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
# gmsh.model.mesh.refine()
# gmsh.model.mesh.refine()


# RecombineMesh
if !("-nopopup" in ARGS)
    gmsh.fltk.run()
end

msh_file = "porous_media.msh" |> tcf;
print(msh_file)
gmsh.write(msh_file)

gmsh.finalize()

##
using Gridap
using GridapGmsh

model = GmshDiscreteModel(msh_file);

vtk_file = "model" |> tcf;
writevtk(model, vtk_file);
