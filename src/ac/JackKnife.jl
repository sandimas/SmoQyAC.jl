function JackKnife(datums::AbstractArray{T}) where {T}
    nbin = size(datums,2)
    avg = mean(datums, dims=2)
    
    datums = datums .- avg
    err = sum(datums .* datums,dims = 2) ./ (nbin * ( nbin -1 ))
    err = sqrt.(err)
    return avg, err
end