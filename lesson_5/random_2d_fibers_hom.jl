using Gmsh

gmsh.initialize()

##
gmsh.clear()

function addFiber(cx, cy, th, L, θ, tag = -1)
    x = cx - th/2;
    y = cy - L/2;
    dx = th;
    dy = L;

    v_corner = [cx, cy, 0];
    v_direction = [0, 0, 1];

    rect_tag = gmsh.model.occ.addRectangle(x, y, 0, dx, dy, tag);
    gmsh.model.occ.synchronize();
    gmsh.model.occ.rotate([(2, rect_tag)], v_corner..., v_direction..., θ);
    return rect_tag;
end

fibers = []
for i in 1:10
    cx, cy = rand(2)
    th = 0.1;
    θ = randn();
    L = 2;

    push!(fibers, addFiber(cx, cy, th, L, θ))
end
first = [(2, tag) for tag in fibers[1:5]]
second = [(2, tag) for tag in fibers[6:10]]

tags = gmsh.model.occ.fuse(first, second)
gmsh.fltk.run()
##
gmsh.finalize()