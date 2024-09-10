using Gmsh

gmsh.initialize()


lc = 5e-2

factory = gmsh.model.occ


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
for i in 1:10
    xy = rand(2) .- 0.5;;
    r = 0.1;
    circle = factory.addCircle(xy[1], xy[2], 0, r);
    inner_loop = factory.addCurveLoop([circle]);
    hole = factory.addPlaneSurface([inner_loop]);
    push!(holes, (2, hole));
end


diff = factory.cut([(2,full_domain)], holes, -1, true, true)

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

gmsh.model.occ.synchronize()

gmsh.model.mesh.generate(2)

if !("-nopopup" in ARGS)
    gmsh.fltk.run()
end

gmsh.finalize()