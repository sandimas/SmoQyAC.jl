using Statistics
using FileIO

function deac_2D_generate_input_files(SimulationFolder::String,Correlation::String,inputData::AbstractArray,inputError::AbstractArray,β::AbstractFloat)

    nτ = size(inputData,1) 
    nx = size(inputData,2)
    ny = size(inputData,3)
    Δτ = β / (nτ - 1)
    τs = collect(range(0.0,(Float64(nτ)-1)*Δτ,nτ))
    
    fermion = (Correlation == "greens_up" || Correlation == "greens_dn" || Correlation ==  "greens") 
    nτ2 =  (fermion) ? nτ : Int64((nτ-1)/2)+1
    # nτ2 = Int64((nτ-1)/2)+1
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

function run_deac_AC_2D(SimulationFolder::String,Correlation::String,nx,ny,β;nStatistics=1000,baseSeed=1000,ω_max=20,nω=100, symmetry::String="none",deac_exe_dir::String="",bin_size=100)
    fermion = (Correlation == "greens_up" || Correlation == "greens_dn" || Correlation ==  "greens") 
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
    for stat in 1:nStatistics
        Threads.@threads for pair in 1:sym["pairlen"]
            
            x = sym["xvec"][pair]
            y = sym["yvec"][pair]
            
            fname = string(x) * "_" * string(y)
            try
                mkdir(save_dir * fname)
            catch
            end
            seed = baseSeed + stat
            # for seed = 1+baseSeed:nStatistics+baseSeed
            flags = split(default_flags * " --save_directory "*save_dir *fname *" --seed "* string(seed)*" "*in_dir*fname*".bin")
            Base.run(`$executable $flags`)
            if stat % bin_size == 0
                bin_DEAC_data(SimulationFolder,Correlation,x,y,bin_size)
            end
        end
        
    end
    
end

function load_from_deac(SimulationFolder::String,Correlation::String,nx::Int64,ny::Int64;symmetry="none",delete_files=false)
    out_dir = SimulationFolder * "/AC_out/"*Correlation*"/DEAC/"
    dir_f = out_dir * "1_1/"
    ωs = load(dir_f*"omega.jld2")["ω"]
    nω = size(ωs,1)
    data = zeros(Float64,(nω,nx,ny))
    N_data = zeros(Int64,(nx,ny))
    
    sym = get_symmetry_vecs(nx,ny,symmetry)

    # Threads.@threads 
    for pair in 1:sym["pairlen"]
        x = sym["xvec"][pair]
        y = sym["yvec"][pair]

        dir = SimulationFolder * "/AC_out/"*Correlation*"/DEAC/" * string(x) * "_" * string(y) * "/"
        
        flist = filter(x->contains(x,"bin_"),readdir(dir))
        nfile = size(flist,1)
        datums = zeros(Float64,(nω,nfile))
        err = zeros(Float64,(nω,))
        outdat = Vector{Float64}(undef,nω)
        n_total = 0
        for fn in 1:nfile
            bindict = load(dir*flist[fn])
            
            datums[:,fn]= bindict["bin"]
            n_total += bindict["N"]
        end
        N_data[x,y] = n_total
        if nfile > 1
            outdat, err = JackKnife(datums)
        else
            outdat = datums[:,1]
        end
        data[:,x,y] = outdat
        
        if delete_files 
            for file in flist
                rm(dir*file)
            end
        end
    end
    ### Symmetry reflections
    if symmetry != "none"
        populate_by_symmetry!(data,symmetry)
    end
    data_dict = Dict{String,Any}(
        "A" => data,
        "ωs" => ωs,
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


function bin_DEAC_data(SimulationFolder::String,Correlation::String,kx::Int64,ky::Int64,bin_size::Int64)
    
    folder = SimulationFolder * "/AC_out/"*Correlation*"/DEAC/"*string(kx)*"_"*string(ky)*"/"
    bin_freq_files = filter(x->contains(x,"omega.jld2"),readdir(folder))
    flist_f = filter(x->contains(x,"_frequency_"),readdir(folder))
    nω = Int64(filesize(folder*flist_f[1])/sizeof(Float64))
    
    if size(bin_freq_files,1) == 0
        freqlist = Vector{Float64}(undef,nω)
        ff = open(folder*flist_f[1])
        read!(ff,freqlist)
        close(ff)
        freqdict = Dict{String,Any}( "ω" => freqlist)
        save(folder*"omega.jld2",freqdict)
    end
    # remove freq files
    for file in flist_f
        rm(folder*file)
    end
    # remove logs
    flist_f = filter(x->contains(x,"_log_"),readdir(folder))
    for file in flist_f
        rm(folder*file)
    end

    # bin data
    flist_f = filter(x->contains(x,"_dsf_"),readdir(folder))
    if size(flist_f,1) < bin_size
        println("Error binning in "*folder)
        println("Only "*string(size(flist_f,1))*" of "*string(bin_size)*" dsf files found")
        exit()
    end
    data = zeros(Float64,(nω,))
    for file in flist_f
        tmp = Array{Float64}(undef,(nω,))
        ff = open(folder*file)
        read!(ff,tmp)
        close(ff)
        data += tmp
    end
    data = data ./ bin_size
    bins = filter(x->contains(x,"bin_"),readdir(folder))
    binnum = size(bins,1) + 1
    datadict = Dict{String,Any}( "bin" => data, "N" => bin_size)
    save(folder*"bin_"*string(binnum)*".jld2",datadict)
    for file in flist_f
        rm(folder*file)
    end
end
