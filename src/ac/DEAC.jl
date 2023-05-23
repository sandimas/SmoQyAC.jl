using Statistics


function deac_2D_generate_input_files(SimulationFolder::String,Correlation::String,inputData::AbstractArray,inputError::AbstractArray,β::AbstractFloat)

    nτ = size(inputData,1) 
    nx = size(inputData,2)
    ny = size(inputData,3)
    Δτ = β / (nτ - 1)
    τs = collect(range(0.0,(Float64(nτ)-1)*Δτ,nτ))
    
    # fermion = (Correlation == "greens_up" || Correlation == "greens_dn") 
    nτ2 = Int64((nτ-1)/2)+1 # (fermion) ? nτ : Int64((nτ-1)/2)+1
    deac_dir = SimulationFolder * "/deac_inputs/" 
    try
        mkdir(deac_dir)
        
    catch
    end
    try
        mkdir(deac_dir * Correlation * "/")
    catch
    end
    for x in 1:nx
        for y in 1:ny
            out_arr = zeros(Float64,(nτ2,3))
            out_arr[:,1] = τs[1:nτ2]
            out_arr[:,2] = inputData[1:nτ2,x,y]
            out_arr[:,3] = inputError[1:nτ2,x,y]
            fname = deac_dir * Correlation *"/"*string(x) * "_" * string(y) * ".bin"
            open(fname,"w") do file
                write(file,out_arr)
            end
            # read!(f,freqlist)
    
        end
    end
    return nx, ny
end

function run_deac_AC_2D(SimulationFolder::String,Correlation::String,nx,ny,β;nStatistics=1000,baseSeed=1000,ω_max=20,nω=100, symmetry::String="none",deac_exe_dir::String="")
    fermion = (Correlation == "greens_up" || Correlation == "greens_dn") 
    executable = (fermion) ? "deac.f" : "deac.b"
    if deac_exe_dir != ""
        executable = deac_exe_dir * "/" * executable
    end
    save_dir = SimulationFolder * "/AC_out/" 
    in_dir = SimulationFolder * "/deac_inputs/"*Correlation*"/"
    temperature = 1/β
    try
        mkdir(save_dir)
    catch
    end
    try
        save_dir *= Correlation
        mkdir(save_dir)
    catch
    end
    try
        save_dir *= "/DEAC/"
        mkdir(save_dir)
    catch
    end
    default_flags = "--number_of_generations 1600000 --temperature "*string(temperature)* 
                    " --population_size 8 --genome_size "*string(nω)*
                    " --normalize --omega_max " *string(ω_max)* " --stop_minimum_fitness 1.0"

    ### Symmetry stuff
    sym = get_symmetry_vecs(nx,ny,symmetry)


    Threads.@threads for thd in 1:sym["pairlen"]*nStatistics
        seed_offset = ceil(Int,thd/nStatistics)
        pair = thd % sym["pairlen"] +1
        x = sym["xvec"][pair]
        y = sym["yvec"][pair]
        fname = string(x) * "_" * string(y)
        try
            mkdir(save_dir * fname)
        catch
        end
        seed = baseSeed + seed_offset
        # for seed = 1+baseSeed:nStatistics+baseSeed
        flags = split(default_flags * " --save_directory "*save_dir *fname *" --seed "* string(seed)*" "*in_dir*fname*".bin")
        Base.run(`$executable $flags`)
        # end
    end
    
end

function load_from_deac(SimulationFolder::String,Correlation::String,nx::Int64,ny::Int64;symmetry="none")
    
    dir_f = SimulationFolder * "/AC_out/"*Correlation*"/DEAC/1_1/"
    flist_f = filter(x->contains(x,"_frequency_"),readdir(dir_f))
    nω = Int64(filesize(dir_f*flist_f[1])/sizeof(Float64))
    freqlist = Vector{Float64}(undef,nω)
    ff = open(dir_f*flist_f[1])
    read!(ff,freqlist)
    close(ff)
    out_dir = SimulationFolder * "/AC_out/"*Correlation*"/DEAC/"
    data = zeros(Float64,(nω,nx,ny))
    N_data = zeros(Int64,(nx,ny))
    
    sym = get_symmetry_vecs(nx,ny,symmetry)

    Threads.@threads for pair in 1:sym["pairlen"]
        x = sym["xvec"][pair]
        y = sym["yvec"][pair]

        dir = SimulationFolder * "/AC_out/"*Correlation*"/DEAC/" * string(x) * "_" * string(y) * "/"
        
        flist = filter(x->contains(x,"_dsf_"),readdir(dir))
        nfile = size(flist,1)
        datums = zeros(Float64,(nω,nfile))
        err = zeros(Float64,(nω,))
        outdat = Vector{Float64}(undef,nω)
        N_data[x,y] = nfile
        for fn in 1:nfile
            f = open(dir*flist[fn])
            indat = Vector{Float64}(undef,nω)
            read!(f,indat)
            datums[:,fn]= indat
            close(f)
        end

        if nfile > 1
            outdat, err = JackKnife(datums)
        else
            outdat = datums[:,1]
        end
        data[:,x,y] = outdat
        
        
        
        outfile = open(out_dir*string(x) * "_" * string(y)*".csv","w")
        write(outfile,"OMEGA MEAN_I STD\n")
        for w in 1:nω
            write(outfile,string(freqlist[w]) * " " * string(outdat[w]) * " " * string(err[w]) * "\n")
        end
        close(outfile)
        
    end
    ### Symmetry reflections
    if symmetry != "none"
        populate_by_symmetry!(data,symmetry)
    end
    data_dict = Dict{String,Any}(
        "A" => data,
        "ωs" => freqlist,
        "N" => N_data
    )
    return data_dict
end

function merge_DEAC_outputs(data_1::Dict{String,Any},data_2::Dict{String,Any})
    n_total = data_1["N"] +  data_2["N"]
    A1 = data_1["A"]
    A2 = data_2["A"]
    nω = size(data_1["ωs"],1)
    merged_data = zeros(Float64,size(A1))
    for ω in 1:nω
        merged_data[ω,:,:] = (data_1["N"] .* A1[ω,:,:] .+ data_2["N"] .* A2[ω,:,:] ) ./ n_total
    end
    data_out = Dict{String,Any}(
        "A" => merged_data,
        "ωs" => data_1["ωs"],
        "N" => n_total
    )
    return data_out
end
