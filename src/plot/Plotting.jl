#!/usr/bin/env julia

using CairoMakie
using ColorSchemes

function Plot_cut(data::AbstractArray,ωs::AbstractArray,outfile_name::String;
                  ω_low::Float64=0.0,ω_high::Float64=2.0,
                  title::String="",x_label::String="",y_label::String="ω",
                  xtick_info::Dict{String,Any}=(),overlay_1D_data::AbstractArray=[],
                  clip_vals::Vector{Float64}=[])

    figure_info = (; resolution=(1000,600))
    axis = (; xlabel=x_label, ylabel=y_label,title=title)
    ω_low_index = findfirst(ωs .>= ω_low)
    ω_high_index = findlast(ωs .<= ω_high)
    x_positions = collect(range(0.0,1.0,size(data,1)))
    

    if !isempty(xtick_info)
        x_positions = xtick_info["x_positions"]
        x_string = xtick_info["x_tick_labels"]
        x_tick_position = xtick_info["x_tick_position"]
    end
    
    if isempty(clip_vals)
        fig, ax, pltobj = heatmap(x_positions,ωs[ω_low_index:ω_high_index],data[:,ω_low_index:ω_high_index]; axis=axis, figure=figure_info)
    else
        
        fig, ax, pltobj = heatmap(x_positions,ωs[ω_low_index:ω_high_index],data[:,ω_low_index:ω_high_index]; axis=axis, figure=figure_info, colorrange=clip_vals, highclip=:red, lowclip=:black)
    end
    if !isempty(xtick_info)
        ax.xticks = (x_tick_position,x_string)
    end
    if !isempty(overlay_1D_data)
        scatter!(ax,x_positions,overlay_1D_data; color=:white, marker=:star4, markersize=20)
    end



    Colorbar(fig[1,2], pltobj)
    current_figure()
    save(outfile_name,fig)

end