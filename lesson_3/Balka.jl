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

gmsh.model.add("model")

# Load a STEP file (using `importShapes' instead of `merge' allows to directly
# retrieve the tags of the highest dimensional imported entities):
path = "Balka.stp" |> tcf
v = gmsh.model.occ.importShapes(path)

# If we had specified
#
# gmsh.option.setString('Geometry.OCCTargetUnit', 'M')
#
# before merging the STEP file, OpenCASCADE would have converted the units to
# meters (instead of the default, which is millimeters).

# Get the bounding box of the volume:
xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(
    v[1][1], v[2][2])

    zmin=-13.0
gmsh.model.occ.synchronize()
for plane in gmsh.model.getEntities(2)
    println(gmsh.model.getColor(plane...))
    gmsh.model.addPhysicalGroup(2, [plane[2]], -1, "$(plane[2])");

end



colors4 = [gmsh.model.getColor(plane...) for plane in gmsh.model.getEntities(2)]



d3 = gmsh.model.getEntities(2)
gmsh.model.addPhysicalGroup(2, [d3[1][2]], -1, "down1")
gmsh.model.addPhysicalGroup(2, [d3[3][2]], -1, "down2")
gmsh.model.occ.synchronize()


gmsh.model.mesh.generate(3)
msh_file = "model.msh" |> tcf;
gmsh.write(msh_file)

# gmsh.fltk.run()

gmsh.finalize()


using Gridap
using GridapGmsh

model = GmshDiscreteModel(msh_file);

vtk_file = "model" |> tcf;
writevtk(model, vtk_file);