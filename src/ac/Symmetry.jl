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