
function get_2D_time_displaced(SimFolder::String,GreensType::String;orbital=1,space::String="momentum")
    _,real,_,std, β = load_from_SmoQyDQMC(simulationfolder=SimFolder,
                                       correlation=GreensType,
                                       space=space,
                                       type="time-displaced")
    real_2D = real[orbital,orbital,:,:,:,1]
    err_2D = std[orbital,orbital,:,:,:,1]

    return real_2D, err_2D, β
end

function get_2D_integrated(SimFolder::String,GreensType::String;orbital=1)
    _,real,_,std, β = load_from_SmoQyDQMC(simulationfolder=SimFolder,
                                       correlation=GreensType,
                                       space="momentum",
                                       type="integrated")
    real_2D = real[orbital,orbital,1,:,:,1]
    err_2D = std[orbital,orbital,1,:,:,1]

    return real_2D, err_2D, β
end

function get_DOS_AC_input(SimFolder::String,GreensType::String;orbital=1)
    _,real,_,std, β = load_from_SmoQyDQMC(simulationfolder=SimFolder,
                                        correlation=GreensType,
                                        space="position",
                                        type="time_displaced")
    real_2D = real[orbital,orbital,:,1,1,1]
    err_2D = std[orbital,orbital,:,1,1,1]

    return real_2D, err_2D, β
end

function get_2D_for_renormalized_Ω(SimFolder::String,Ω_0::Float64;orbital=1)
    _,real,_,std, β = load_from_SmoQyDQMC(simulationfolder=SimFolder,
                                        correlation="phonon_greens",
                                        space="momentum",
                                        type="integrated")
    real_2D = real[orbital,orbital,1,:,:,1]
    err_2D = std[orbital,orbital,1,:,:,1]
    return real_2D, err_2D, β
end


