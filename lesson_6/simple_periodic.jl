using WriteVTK
using LinearAlgebra
using SparseArrays

tcf(filename::String) = (@__DIR__) * "\\" * filename;

H = 1;
W = 1;
num = 100;

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
# rank(mK)


fixed = findall(points[1,:] .== 0)
loaded = findall(points[1,:] .== 1)


function applydirichelet!(mK::AbstractMatrix{<:Number}, vF::AbstractVector{<:Number},ind::AbstractVector{<:Integer}, dirs::BitVector, func::Function, points::AbstractMatrix{<:Numbers})
    m = size(dofs,1);
    for (i,s) in dirs
        if s 
            for j in ind
                p = points[:, j];
                mK[m*j - 3 + i, m*j - 3 + i] = 1;

                vF[m*j - 3 + i] = 2*func(p)[i];
                
                for k in eachindex(vF)
                    if mK[m*j - 3 + i, k] != 0
                        vF[k] -= func(p)[i] * mK[m*j - 3 + i, k];
                    end
                end 
            end
        end
    end

end

for p in fixed
    mK[2p, 2p] += 1e16;
    mK[2p-1, 2p-1] += 1e16;
end

# rank(mK)

for p in loaded
    vF[2p] = 100;
end
using IterativeSolvers
using LinearAlgebra

δ = cg(mK,vF)

stresses = Matrix{Float64}(undef, 3, length(cells))
for i in 1:size(connections,2)
    stresses[:, i] = stress(points[:, connections[:,i]], reshape(δ, 2, num^2)[:, connections[:,i]]);
end

vtk_grid(filename, points, cells) do vtk

    disps = zeros(3, num^2);
    disps[1:2,:] = reshape(δ, 2, num^2)
    vtk["disps"] = disps;
    vtk["stress"] = stresses;
end
