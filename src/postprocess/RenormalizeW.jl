#!/usr/bin/env julia

function renormalize_Ω(SimulationFolder::String,Ω_0::Float64;orbital=1)
    data, _, _ = get_2D_for_renormalized_Ω(SimulationFolder,Ω_0;orbital=orbital)
    # data = -2 .* data .* Ω_0
    Ω_new = sqrt.(1 ./ data)
    return Ω_new
end

