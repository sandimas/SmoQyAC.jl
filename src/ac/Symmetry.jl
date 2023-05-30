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
            pairlen =   Int64(0.125 * nx * nx + 0.75 * nx + 1)
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


function get_symmetric_point(data,x,y,symmetry)
    Nx = size(data,2)
    Ny = size(data,3)
    Nx2 = Int64(Nx/2)+1
    Ny2 = Int64(Ny/2)+1
    
    # use symmetry if possible, bozo!
    if symmetry == "none"
        return 0,0
    end
    #  Γ, M are singly degenerate
    if ((x==1)&&(y==1)) || ((x==Nx2)&&(y==Ny2))
        return 0,0
    end


    if symmetry == "rectangle"
        
        if ((x>Nx2)&&(y<=Ny2)) # Top right, done
            return 0,0
        end
        # doubly degenerate points
        if (y==1)  
            if x != Nx2
                return 2*Nx2-x,y
            else
                return 0,0
            end
        end
        if (x==1) 
            if (y<Ny2)
                return x,2*Ny2-y
            else 
                return 0,0
            end
        end
        if (x==Nx2) 
            if (y<Ny2)
                return Nx2,2*Ny2-y
            else
                return 0,0
            end
        end
        if (y==Ny2)
            if (x<Nx2) 
                return 2*Nx2-x,Ny2
            else
                return 0,0
            end
        end
        # reflect from top left to bottom left
        if (x < Nx2) && (y< Ny2)
            return x, 2*Ny2 -y
        end
        if (x < Nx2) && (y > Ny2)
            return 2*Nx2-x,y
        end
        if (x > Nx2) && (y > Ny2)
            return x, 2*Ny2-y
        end
        return 0,0

    end

    if symmetry == "square"
        # Deal with x==1 and rotations
        if (y==1) 
            if (x<=Nx2)
                return y,x 
            else
                return 0,0
            end
        end
        if (x==1) 
            if (y<Ny2)
                return x,2*Ny2-y

            elseif (y==Ny2)
                return 0,0
            else
                return y,x
            end
        end
        # deal with x==Nx2 and rotations
        if (x==Nx2)
            return y,x
        end
        if (y==Ny2)
            if (x==1) || (x>Nx2)
                return 0,0
            else
                return Nx2,Ny2+Nx2-x
            end
        end
        # diagonals
        if (x==y)
            return x,2*Ny2-y
        end
        if (Nx2-x==y-Ny2)
            if (x<=Nx2)
                return 2*Nx2-x,y
            else 
                return 0,0
            end
        end
        
        # Pinwheel portion
        # top left
        if (x<Nx2) && (y<Ny2) 
            # First slice
            if (x>y)
                return y,x
            # second slice
            else
                return x,2*Ny2-y
            end

        end
        # bottom left
        if (x<Nx2) && (y>Ny2)
            if (x-1<=Ny-y)
                return 2*Nx2-y,Ny+2-x
            else
                return 2*Nx2-x,y
            end
        end
        # bottom right
        if (y>Ny2) 
            if (x<y)
                return y,x
            else
                return x, 2*Ny2-y
            end
        end
        # top right
        if (x-Nx2>Ny2-y)
            return Nx+2-y,2*Nx2-x
        else
            return 0,0
        end
        

    end
    println("case not found!!!")
    return 0,0

  
end

