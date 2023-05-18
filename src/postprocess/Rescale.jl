using Interpolations
using FourierTools

function Upsample_Interpolations(G::AbstractArray{ComplexF64},GreensType::String,ωs::AbstractArray{Float64},new_nx::Int64,new_ny::Int64;Ω_0::Float64=1.0)
    nω = size(ωs,1)        
    Δx = (size(G,2)-1)/new_nx
    Δy = (size(G,3)-1)/new_ny
        
    if GreensType != "phonon"
        
        
        ϵ_k_old = build_ϵ_k(size(G,2),size(G,3))
        ϵ_k_new = build_ϵ_k(new_nx,new_ny)
        Σ_old = fermion_get_Σ(G,ϵ_k_old,ωs)
        Σ_new = zeros(ComplexF64,(nω,Nx,Ny))
           
        itp = interpolate(Σ_old[:,:,:],BSpline(Quadratic(Reflect(OnCell()))))
    
      
        for ω in 1:ωs
            for kx in 1:new_nx
                for ky in 1:new_ny
                    Σ_new[ω,kx,ky] = itp(Float64(ω),Δx*kx+1.0,Δy*ky+1.0)
                end
            end
        end
        
        G_new = Σ_to_G(Σ_new,ϵ_k_new,ωs)

        A = G_to_A(G_new)
        return G_new, A
    else
     
        Π_old = boson_get_Π(G,ωs,Ω_0)
        Π_new = zeros(ComplexF64,(nω,new_nx,new_ny))
    
        itp = interpolate(Π_old[:,:,:],BSpline(Quadratic(Reflect(OnCell()))))
        
        for ω in 1:nω
            for kx in 1:new_nx
                for ky in 1:new_ny
                    Π_new[ω,kx,ky] = itp(Float64(ω),Δx*kx+1.0,Δy*ky+1.0)
                end
            end

        end
    
        D_new = Π_to_D(Π_new,ωs,Ω_0)
        B = D_to_B(D_new,ωs,Ω_0)
        return D_new, B
    end
    
    
end



function build_ϵ_k(nx::Int64,ny::Int64)
    ϵ_k = zeros(Float64,(nx,ny))
    Δx = 2 * π / (nx)
    Δy = 2 * π / (ny)
    for x in 1:nx
        calc_x = cos((x-1)*Δx)
        for y in 1:ny
            calc_y = cos((y-1)*Δy)
            ϵ_k[x,y] = - 2 *  (calc_x+calc_y)
        end
    end
    return ϵ_k
end


function fermion_get_Σ(G::AbstractArray{ComplexF64},ϵ_k::AbstractArray{Float64},ωs::AbstractArray{Float64})
    Σ = zeros(ComplexF64,size(G))
    nω = size(ωs,1)
    Lx = size(ϵ_k,1)
    Ly = size(ϵ_k,2)
    for ω in 1:nω
        for x in 1:Lx
            for y in 1:Ly
                Σ[:,ω,x,y] .= ωs[ω]-ϵ_k[x,y] .- 1/G[:,ω,x,y]
            end
        end
    end

    return Σ
end

function Σ_to_G(Σ::AbstractArray{ComplexF64},ϵ_k::AbstractArray{Float64},ωs::AbstractArray{Float64})
    G = zeros(ComplexF64,size(Σ))
    nω = size(Σ,1)
    nkx = size(Σ,2)
    nky = size(Σ,3)

    for kx in 1:nkx
        for ky in 1:nky
            for ω in 1:nω
                G[ω,kx,ky] = 1/(ωs[ω]-ϵ_k[kx,ky]-Σ[ω,kx,ky])
            end
        end
    end
    
    return G
end

function G_to_A(G::AbstractArray{ComplexF64})
    A = -(1/π) * imag.(G)
    return A
end

function boson_get_Π(D::AbstractArray{ComplexF64},ωs::AbstractArray{Float64},Ω_0::Float64)
    Π = zeros(ComplexF64,size(D))
    nω = size(ωs,1)
    Lx = size(D,2)
    Ly = size(D,3)
    for ω in 1:nω
        for x in 1:Lx
            for y in 1:Ly
                
                Π[ω,x,y] = (ωs[ω]^2 - (Ω_0)^2) /(2*(Ω_0)) - 1/D[ω,x,y]
            end
        end
    end
    
    return Π
end

function Π_to_D(Π::AbstractArray{ComplexF64},ωs::AbstractArray{Float64},Ω_0::Float64)
    D = zeros(ComplexF64,size(Π))
    nω = size(Π,1)
    nkx = size(Π,2)
    nky = size(Π,3)
    
    for kx in 1:nkx
        for ky in 1:nky
            for ω in 1:nω
                D[ω,kx,ky] = 2*Ω_0/(ωs[ω]^2 - (Ω_0)^2 - 2* Ω_0 * Π[ω,kx,ky])
            end
        end
    end
    return D
end

function D_to_B(D::AbstractArray{ComplexF64},ωs::AbstractArray{Float64},Ω_0::Float64)
    B = zeros(Float64,size(D))
    for ω in 1:size(D,1)
        for kx in 1:size(D,2)
            for ky in 1:size(D,3)
                B[ω,kx,ky] = -(1/π)*ωs[ω]*imag(D[ω,kx,ky])
            end
        end
    end
    return B
end

