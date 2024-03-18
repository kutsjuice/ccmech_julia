using WriteVTK
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using PyPlot

remove(a, b) = findall(!in(b), a)
findin(a, b) = findall(in(b), a)


tcf(filename::String) = (@__DIR__) * "\\" * filename;

H = 1;
W = 1;
num = 50;

x = LinRange(0, W, num);
y = LinRange(0, H, num);

XX = reshape(repeat(x, inner = (1,num) )', num^2);
YY = reshape(repeat(y, inner = (1,num)), num^2);

points = [XX'; YY']

connections = Matrix{Int64}(undef, 3, (num-1)^2 * 2)
k = 1;
for i in 1:num-1
    for j in 1:num-1
        p1 = (i-1) * num + j;
        p2 = p1 + 1;
        p3 = i*num +j +1;
        p4 = p3-1;
        connections[:, k] .= [p1,p3,p2];
        connections[:, k+1] .= [p1,p4,p3];
        k+=2; 
    end
end

cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, connection) for connection in eachcol(connections)]

filename = "geo" |> tcf
vtk_grid(filename, points, cells) do vtk
end

const E = 2.1e11
const ν = 0.3

# Model
flag = 1;
k1 = ν/(1-ν);
k2 = 0.5*(1-2ν)/(1-ν);
k3 = (1-ν)/2;

plain_stress = true;

mE = E/(1-ν^2)* [
    1  ν  0;
    ν  1  0;
    0  0 k3;
];


function calcMKLocal(pts::Matrix{Float64})::Matrix{Float64}
    mB = Matrix{Float64}(undef, 3,6);
    ijk = [1,2,3];
    for _ in eachindex(pts)
        i,j,k = ijk;
        circshift!(ijk,1)
        # a = pts[1,j]*pts[2,k] - pts[1,k]*pts[2,j];
        b = pts[2,j] - pts[2,k];
        c = pts[1,k] - pts[1,j];
        mB[:, 2i-1:2i] .= [
            b 0; 
            0 c; 
            c b
            ];
    end
    Δ2 = det([
        1 pts[1,1] pts[2, 1];
        1 pts[1,2] pts[2, 2];
        1 pts[1,3] pts[2, 3];
        ]); 
    mB /= Δ2;
    mKloc = mB' * mE * mB * Δ2/2;
    return mKloc;
end

function stress(pts::AbstractMatrix{<:Real}, disps::AbstractMatrix{<:Real})::Vector{Float64}
    mB = Matrix{Float64}(undef, 3,6);
    ijk = [1,2,3];
    for _ in eachindex(pts)
        i,j,k = ijk;
        circshift!(ijk,1)
        # a = pts[1,j]*pts[2,k] - pts[1,k]*pts[2,j];
        b = pts[2,j] - pts[2,k];
        c = pts[1,k] - pts[1,j];
        mB[:, 2i-1:2i] .= [
            b 0; 
            0 c; 
            c b
            ];
    end
    Δ2 = det([
        1 pts[1,1] pts[2, 1];
        1 pts[1,2] pts[2, 2];
        1 pts[1,3] pts[2, 3];
        ]); 
    mB /= Δ2;
    σ = (mE * mB * reshape(disps, 6, 1))[:];
    return σ
end

mK = spzeros(num^2 * 2, num^2 * 2);
vF = spzeros(num^2 * 2);

for el in eachcol(connections)
    mKloc = calcMKLocal(points[:, el]);
    inds = [2*el[1]-1, 2*el[1], 2*el[2]-1, 2*el[2], 2*el[3]-1, 2*el[3] ]
    mK[inds, inds] += mKloc;
end
rank(mK)


fixed = findall(points[1,:] .== 0)
loaded = findall(points[1,:] .== 1)

function applydirichelet!(mK::AbstractMatrix{<:Number}, vF::AbstractVector{<:Number},indices::AbstractVector{<:Integer}, dirsmask::BitVector, func::Function, points::AbstractMatrix{<:Number})
    dims = size(dirsmask,1);
    dirs = findall(dirsmask);
    for ind in indices
        
        p = points[:, ind];
        dofsvalue = func(p);

        for dir in dirs
            dof = dims*(ind-1) + dir;            
            mK[dof, dof] = 1;
            vF[dof] = dofsvalue[dir];
            
            for k in eachindex(vF)
                if dof != k
                    if mK[dof, k] != 0
                        vF[k] -= dofsvalue[dir] * mK[dof, k];
                        mK[dof,k] = mK[k, dof] = 0;
                    end
                end
            end
        end 
    end

end
# function applydirichelet!(mK, vF, dofs, value)
    
#     for dof in dofs
#         for 
#     end
# end



f(p) = [0,0,0];
applydirichelet!(mK, vF, fixed, BitVector([1,1]), f, points)
rank(mK)


function applyMPC(  
            mK::AbstractMatrix{<:Number}, 
            vF::AbstractVector{<:Number}, 
            primary_dofs::AbstractVector{<:Integer}, 
            secondary_dofs::AbstractVector{<:Integer}, 
            koefs::AbstractVector{<:Number}, 
            offset::AbstractVector{<:Number}
        )

    @assert size(primary_dofs, 1) == size(secondary_dofs,1)
    @assert size(primary_dofs, 1) == size(koefs,  1)
    @assert size(primary_dofs, 1) == size(offset, 1)
    
    n = size(vF, 1);
    m = n - size(secondary_dofs, 1);
    uncommitted = remove(1:n, [primary_dofs;secondary_dofs]);
    uncommitted_and_primary = remove(1:n, secondary_dofs);

    new_primary = findall(in(primary_dofs), uncommitted_and_primary)
    new_uncommitted = findall(in(uncommitted), uncommitted_and_primary);

    mT = spzeros(n,m);
    vG = spzeros(n);
    for (j,i) in enumerate(uncommitted_and_primary)
        mT[i,j] = 1;
    end
    for el in eachindex(primary_dofs)
        i = secondary_dofs[el];
        j = new_primary[el];
        mT[i,j] = koefs[el];
        vG[i] = offset[el];
    end
    vFnew = mT' * (vF - mK * vG) 
    mKnew = mT' * mK * mT;

    return (mKnew,vFnew, uncommitted_and_primary, new_primary, mT)
end

const EPS = 1e-6

isupper(p) = abs(p[2] - H) < EPS;
islower(p) = abs(p[2]) < EPS;


function finddofs(points::AbstractMatrix{<:Number}, condition::Function, number_of_directions::Integer)
    indices = findall(condition.(eachcol(points)));
    dofs = Vector{Int64}(undef, size(indices, 1) * number_of_directions);
    k = 0;
    for ind in indices
        for direction in 1:number_of_directions
            k += 1; 
            dofs[k] = number_of_directions*(ind-1) + direction;
        end
    end
    return dofs
end

upper_dofs = finddofs(points, isupper, 2) 
lower_dofs = finddofs(points, islower, 2)

upper_dofs
findall(points[2, :] .== 1)

for p in loaded
    vF[2p] = 100;
end

mKn,vFn, uncommitted_and_primary, new_primary, mT = applyMPC(mK, vF, upper_dofs, lower_dofs, ones(size(lower_dofs)), zeros(size(lower_dofs)))

# uncommitted_and_primary
# fixed


# function applydirichelet!(mK, vF, old_dofs, uncommitted_and_primary, primary, secondary, value) 
#     new_dofs = findall(in(uncommitted_and_primary), old_dofs)
#     for dof in new_dofs
        
#     end 
# end

# fixed_dofs = Vector{Int64}(undef, size(fixed,1)*2)
# for i in eachindex(fixed)
#     fixed_dofs[2i-1] = fixed[i]*2 - 1; 
#     fixed_dofs[2i] = fixed[i]*2 
# end 
# fixed_dofs




# res = qr(mKn)

δn = cg(mKn,vFn)

δ = Vector{Float64}(undef, size(vF))
δ[uncommitted_and_primary] = δn;
δ[lower_dofs] = δn[new_primary]

stresses = Matrix{Float64}(undef, 3, length(cells))
for i in 1:size(connections,2)
    stresses[:, i] = stress(points[:, connections[:,i]], reshape(δ, 2, num^2)[:, connections[:,i]]);
end

vtk_grid(filename*"2", points, cells) do vtk

    disps = zeros(3, num^2);
    disps[1:2,:] = reshape(δ, 2, num^2)
    vtk["disps"] = disps;
    vtk["stress"] = stresses;
end





# all_dofs = collect(1:10)

# pr_dofs = [3,4]
# sn_dofs = [1,2]

# un_dofs = remove(all_dofs, [pr_dofs;sn_dofs])

# un_pr_dofs = remove(all_dofs, sn_dofs)
# pr_new_dofs = findall(in(pr_dofs), un_pr_dofs)


# un_new_dofs = findall(in(un_dofs), un_pr_dofs)
