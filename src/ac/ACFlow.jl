#!/usr/bin/env julia
using ACFlow
using CSV
using DataFrames
using Base.Threads

function run_MEM_AC(SimulationFolder::String,
                    Correlation::String,
                    inputData::AbstractArray,
                    inputError::AbstractArray,
                    β::AbstractFloat;
                    ω_max=20.,nω=100, symmetry::String="none")

    fermion = (Correlation == "greens_up" || Correlation == "greens_dn" || Correlation ==  "greens") 
    if fermion
        ω_low = -ω_max 
        ω_high = ω_max
        nω_adj = 2*nω
    else
        ω_low = 0.0
        ω_high = ω_max
        nω_adj = nω
    end

    nτ = size(inputData,1) 
    # nx = size(inputData,2)
    # ny = size(inputData,3)
    # grid = collect(range(0.0,β,nτ))
    # meshω = collect(range(ω_low,ω_high,nω_adj))
    kernel = (fermion) ? "fermi" : "bsymm"
    grid_type =  (fermion) ? "ftime" : "btime"

    B_dict = Dict{String,Any}(
        "solver" => "MaxEnt",
        "ktype"  => kernel,
        "mtype"  => "flat",
        "mesh"   => "linear",
        "grid"   => grid_type,
        "nmesh"  => nω_adj,
        "ngrid"  => nτ,
        "wmax"   => ω_high,
        "wmin"   => ω_low,
        "beta"   => β,
    )

    S_dict = Dict{String,Any}( )
    data_out = run_ACFlow(SimulationFolder,Correlation,inputData,inputError,B_dict,S_dict,symmetry=symmetry)
    return data_out
end

function run_SK_AC(SimulationFolder::String,
                    Correlation::String,
                    inputData::AbstractArray,
                    inputError::AbstractArray,
                    β::AbstractFloat;
                    ω_max=20.,nω=100, symmetry::String="none")

    fermion = (Correlation == "greens_up" || Correlation == "greens_dn" || Correlation ==  "greens") 
    if fermion
    ω_low = -ω_max 
    ω_high = ω_max
    nω_adj = 2*nω
    else
    ω_low = 0.0
    ω_high = ω_max
    nω_adj = nω
    end

    nτ = size(inputData,1) 
    # nx = size(inputData,2)
    # ny = size(inputData,3)
    # grid = collect(range(0.0,β,nτ))
    # meshω = collect(range(ω_low,ω_high,nω_adj))
    kernel = (fermion) ? "fermi" : "bsymm"
    grid_type =  (fermion) ? "ftime" : "btime"

    B_dict = Dict{String,Any}(
    "solver" => "StochSK",
    "ktype"  => kernel,
    "mtype"  => "flat",
    "mesh"   => "linear",
    "grid"   => grid_type,
    "nmesh"  => nω_adj,
    "ngrid"  => nτ,
    "wmax"   => ω_high,
    "wmin"   => ω_low,
    "beta"   => β,
    )

    S_dict = Dict{String,Any}( )
    data_out = run_ACFlow(SimulationFolder,Correlation,inputData,inputError,B_dict,S_dict,symmetry=symmetry)
    return data_out
end

function run_OM_AC(SimulationFolder::String,
    Correlation::String,
    inputData::AbstractArray,
    inputError::AbstractArray,
    β::AbstractFloat;
    ω_max=20.,nω=100, symmetry::String="none")

    fermion = (Correlation == "greens_up" || Correlation == "greens_dn" || Correlation ==  "greens") 
    if fermion
        ω_low = -ω_max 
        ω_high = ω_max
        nω_adj = 2*nω
    else
        ω_low = 0.0
        ω_high = ω_max
        nω_adj = nω
    end

    nτ = size(inputData,1) 
    # nx = size(inputData,2)
    # ny = size(inputData,3)
    # grid = collect(range(0.0,β,nτ))
    # meshω = collect(range(ω_low,ω_high,nω_adj))
    kernel = (fermion) ? "fermi" : "bsymm"
    grid_type =  (fermion) ? "ftime" : "btime"

    B_dict = Dict{String,Any}(
        "solver" => "StochOM",
        "ktype"  => kernel,
        "mtype"  => "flat",
        "mesh"   => "linear",
        "grid"   => grid_type,
        "nmesh"  => nω_adj,
        "ngrid"  => nτ,
        "wmax"   => ω_high,
        "wmin"   => ω_low,
        "beta"   => β,
    )

    S_dict = Dict{String,Any}( )
    data_out = run_ACFlow(SimulationFolder,Correlation,inputData,inputError,B_dict,S_dict,symmetry=symmetry)
    return data_out
end

function run_ACFlow(SimulationFolder::String,
                    Correlation::String,
                    inputData::AbstractArray,
                    inputError::AbstractArray,
                    B_dict::Dict{String,Any},
                    S_dict::Dict{String,Any};
                    symmetry::String="none")
    
    output_dir = SimulationFolder *"/AC_out/"
    try
        mkdir(output_dir)
    catch
    end
    try
        output_dir *= Correlation * "/"
        mkdir(output_dir)
    catch
    end
    try
        output_dir *= B_dict["solver"]*"/"
        mkdir(output_dir)
    catch
    end
    
    nτ = B_dict["ngrid"]
    β = B_dict["beta"]
    nx = size(inputData,2)
    ny = size(inputData,3)
    grid = collect(range(0.0,β,nτ))
    meshω = collect(range(B_dict["wmin"],B_dict["wmax"],B_dict["nmesh"]))
    nω_adj = B_dict["nmesh"]
    ACFlow.setup_param(B_dict,S_dict)
    cur_dir = pwd()
    cd(mktempdir(cur_dir))
    A_out = zeros(Float64,(nω_adj,nx,ny))
    G_out = zeros(ComplexF64,(nω_adj,nx,ny))

    sym = get_symmetry_vecs(nx,ny,symmetry)
    
    Threads.@threads for pair in 1:sym["pairlen"]
        x = sym["xvec"][pair]
        y = sym["yvec"][pair]
        try
            _, A_out[:,x,y], G_out[:,x,y] = ACFlow.solve(grid,inputData[:,x,y],inputError[:,x,y])
        catch
            fail = true
            println("fail ",x," ",y," ",B_dict["solver"])
            xnew,ynew = get_symmetric_point(inputData,x,y)
            while (xnew != 0 && fail)
                try
                    _, A_out[:,x,y], G_out[:,x,y] = ACFlow.solve(grid,inputData[:,xnew,ynew],inputError[:,xnew,ynew])
                    println(xnew," ",ynew," worked for ",x," ",y)
                    fail = false
                catch
                    println("fail ",xnew," ",ynew," ",B_dict["solver"])
                end
                xnew,ynew = get_symmetric_point(inputData,xnew,ynew)
            end
            ### TODO, if fail try using symmetry
        end            
            
#
#        outfile = open("../"*output_dir*string(x) * "_" * string(y)*".csv","w")
#        write(outfile,"OMEGA A GREENS_R GREENS_I\n")
#        for w in 1:nω_adj
#            write(outfile,string(meshω[w]) * " " * string(A_out[w,x,y]) * " " * string(real(G_out[w,x,y])) * " " * string(imag(G_out[w,x,y])) * "\n")
#        end
#        close(outfile)
    end
    cd("../")
    
    ### Symmetry reflections
    if symmetry != "none"
        populate_by_symmetry!(A_out,symmetry)
        populate_by_symmetry!(G_out,symmetry)
    end
    data_out = Dict{String,Any}(
        "A" => A_out,
        "ωs" => meshω,
        "G" => G_out
    )
    return data_out
end

# function load_from_MEM(SimulationFolder::String,Correlation::String,nx::Int64,ny::Int64)

#     ωfile = SimulationFolder * "/AC_out/"*Correlation*"/MaxEnt/1_1.csv"
#     df = CSV.read(ωfile,DataFrame)
#     freqlist = df[:,"OMEGA"]
#     nω = size(freqlist,1)
#     A = zeros(Float64,(nω,nx,ny))
#     G = zeros(ComplexF64,(nω,nx,ny))
#     for x in 1:nx
#         for y in 1:ny
#             file = SimulationFolder * "/AC_out/"*Correlation*"/MaxEnt/" * string(x) * "_" * string(y) * ".csv"
            
#             df = CSV.read(file,DataFrame)
#             A[:,x,y] = df[:,"A"]
#             G[:,x,y] = df[:,"GREENS_R"] + 1im * df[:,"GREENS_I"]
            
            
#         end
#     end
#     data_out = Dict{String,Any}(
#         "A" => A,
#         "ωs" => freqlist,
#         "G" => G
#     )
#     return data_out
# end


