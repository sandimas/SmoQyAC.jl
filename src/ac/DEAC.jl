using Statistics


function deac_2D_generate_input_files(SimulationFolder::String,Correlation::String,inputData::AbstractArray,inputError::AbstractArray,β::AbstractFloat)

    nτ = size(inputData,1) 
    nx = size(inputData,2)
    ny = size(inputData,3)
    Δτ = β / (nτ - 1)
    τs = collect(range(0.0,(Float64(nτ)-1)*Δτ,nτ))
    nτ2 = Int64((nτ-1)/2)+1
    deac_dir = SimulationFolder * "/deac_inputs/" 
    try
        mkdir(deac_dir)
        mkdir(deac_dir * Correlation * "/")
    catch
        print("")
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
            
        end
    end
    return nx, ny
end

function run_deac_AC_2D(SimulationFolder::String,Correlation::String,nx,ny,β;nStatistics=1000,baseSeed=1000,ω_max=20,nω=100)
    fermion = (Correlation == "greens_up" || Correlation == "greens_dn") 
    executable = (fermion) ? "deac.f" : "deac.e"
    save_dir = SimulationFolder * "/AC_out/" 
    in_dir = SimulationFolder * "/deac_inputs/"
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

    for x in 1:nx
        for y in 1:ny
            fname = string(x) * "_" * string(y)
            try
                mkdir(save_dir * fname)
            catch
            end
            for seed = 1+baseSeed:nStatistics+baseSeed
                flags = split(default_flags * " --save_directory "*save_dir *fname *" --seed "* string(seed)*" "*in_dir*fname*".bin")
                Base.run(`$executable $flags`)
            end
        end
    end
end

function load_from_deac(SimulationFolder::String,Correlation::String,nx::Int64,ny::Int64)
    
    dir = SimulationFolder * "/AC_out/"*Correlation*"/DEAC/1_1/"
    flist = filter(x->contains(x,"deac_freq"),readdir(dir))
    nω = Int64(filesize(dir*flist[1])/sizeof(Float64))
    out_dir = SimulationFolder * "/AC_out/"*Correlation*"/DEAC/"
    data = zeros(Float64,(nω,nx,ny))
    freqlist = Vector{Float64}(undef,nω)
    for x in 1:nx
        for y in 1:ny
            dir = SimulationFolder * "/AC_out/"*Correlation*"/DEAC/" * string(x) * "_" * string(y) * "/"
            
            flist = filter(x->contains(x,"deac_dsf"),readdir(dir))
            nfile = size(flist,1)
            datums = zeros(Float64,(nω,nfile))
            err = zeros(Float64,(nω,))
            outdat = Vector{Float64}(undef,nω)
            
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
            flist = filter(x->contains(x,"deac_frequency"),readdir(dir))
            f = open(dir*flist[1])
            read!(f,freqlist)
            close(f)
            outfile = open(out_dir*string(x) * "_" * string(y)*".csv","w")
            write(outfile,"OMEGA MEAN_I STD\n")
            for w in 1:nω
                write(outfile,string(freqlist[w]) * " " * string(outdat[w]) * " " * string(err[w]) * "\n")
            end
            close(outfile)
        end
    end
    return data, freqlist
end