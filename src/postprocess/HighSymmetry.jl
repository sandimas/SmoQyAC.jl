#!/usr/bin/env julia
function get_high_symmetry_2D(data_dict::Dict{String,Any}) 
    A = data_dict["A"]
    nω = size(A,1)
    L = size(A,2)
    L2 = Int64(L/2)
    points = 3*L2+1
    
    data=zeros(Float64,(points,nω))
    for i in 1:L2+1
        data[i,:] = A[:,i,1]
    end
    for i in 1:L2
        data[i+L2+1,:] =A[:,L2+1,i+1]
    end
    for i in 1:L2
        data[i+2*L2+1,:] = A[:,L2+1-i,L2+1-i] 
    end
    xtick_info = Dict{String,Any}(
        "x_tick_position" => [0.001,0.333,0.666,0.999],
        "x_tick_labels" => ["Γ","X","M","Γ"],
        "x_positions" => collect(range(0.0,1.0,size(data,1))),
    )
    return data, xtick_info
end


function get_high_symmetry_1D(A::AbstractArray{T}) where {T}
    L = size(A,1)
    L2 = Int64(L/2)
    points = 3*L2+1
    
    data=zeros(T,(points,))
    for i in 1:L2+1
        data[i] = A[i,1]
    end
    for i in 1:L2
        data[i+L2+1] =A[L2+1,i+1]
    end
    for i in 1:L2
        data[i+2*L2+1] = A[L2+1-i,L2+1-i] 
    end
    xtick_info = Dict{String,Any}(
        "x_tick_position" => [0.001,0.333,0.666,0.999],
        "x_tick_labels" => ["Γ","X","M","Γ"],
        "x_positions" => collect(range(0.0,1.0,size(data,1))),
    )

    return data, xtick_info
end