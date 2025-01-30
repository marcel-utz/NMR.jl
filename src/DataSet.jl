## (c)2025 Marcel Utz

abstract type CoordMap{T} <: AbstractVector{T} end

struct SpectData{T,N} <: AbstractArray{T,N}
    dat::AbstractArray{T,N}
    coord::NTuple{N,AbstractVector}
end

import Base.size
import Base.getindex
import Base.setindex!
import Base.IndexStyle

size(S::SpectData) = size(S.dat)
getindex(S::SpectData, k::Integer) = getindex(S.dat,k)
setindex!(S::SpectData, v, k::Integer) = setindex!(S.dat,v,k)
IndexStyle(S::SpectData) = IndexStyle(S.dat)

coords(S::SpectData) = S.coord
coords(S::SpectData,k::Integer) = S.coord[k]
