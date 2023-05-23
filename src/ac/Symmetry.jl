#!/usr/bin/env julia

function populate_by_symmetry!(Data,Symmetry::String)
    nx = size(Data,2)
    ny = size(Data,3)
    nx2 = Int64(nx/2)
    ny2 = Int64(ny/2)
    if Symmetry == "square"
        for x in 2:nx2+1
            for y in 1:x-1
                Data[:,y,x] = Data[:,x,y]
            end
        end
    end
    Data[:,1:nx2 + 1,2+ny2:ny] = reverse(Data[:,1:nx2 + 1,2:ny2],dims=3)
    Data[:,nx2 + 2:nx,1:ny] = reverse(Data[:,2:nx2,1:ny],dims=2)
    return nothing
end

function get_symmetry_vecs(nx,ny,symmetry)
    xlow = 1
    ylow = 1    
    if symmetry == "none"
        xhigh = nx
        yhigh = ny
        pairlen = nx * ny
    elseif symmetry == "square" || symmetry == "rectangle"        
        xhigh = Int64(nx/2) + 1
        yhigh = Int64(ny/2) + 1
        if symmetry == "rectangle"
            pairlen = Int64(0.25*nx * ny + 0.5*(nx + ny) + 1)
        else 
            pairlen = Int64(0.125 * nx * nx + 0.75 * nx + 1)
        end
    end
    x_vec = zeros(Int64,(pairlen,))
    y_vec = zeros(Int64,(pairlen,))
    pos = 1
    for x in xlow:xhigh
        if (symmetry == "square") 
            yhigh = x
        end
        for y in ylow:yhigh
            x_vec[pos] = x
            y_vec[pos] = y
            pos += 1
            
        end
    end
    sym_dict = Dict{String,Any}(
        "xvec" => x_vec,
        "yvec" => y_vec,
        "pairlen" => pairlen
    )
    return sym_dict
end