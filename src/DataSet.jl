#module DataSet

# using AbstractFFTs
using FFTW
using Statistics

export Data1D,ind2pos,pos2ind,val,ind,cut
export FourierTransform,PhaseCorrect,BaseLineCorrect
export integral,derivative,set!
export plot,plot!
export length,shift,hard_shift
export interp

ENV["GKS_ENCODING"]="utf-8"
import Plots

import Base.iterate

#import NMR.SimplePlot


# Data structures for storing digitised data sets
# Definitions:
#
# - Index: the horizontal index of a data point
# - Position: the numerical position (integer) of the data point in
#           the internal storage array



mutable struct Data1D{Tdata,Tindex}
   dat::Array{Tdata,1}
   istart::Tindex
   istop::Tindex
end

Base.iterate(d::NMR.Data1D,state=1) =
    state > length(d) ? nothing :
      ((d.istart+state*(d.istop-d.istart)/(length(d)-1),d.dat[state]),state+1)

function ind2pos(d::Data1D,ind)
    return round(Int64,(ind-d.istart)*(length(d.dat)-1)/(d.istop-d.istart))+1
end

function pos2ind(d::Data1D,pos::Integer)
   return ind(d)[pos]
end

function ind(d::Data1D)
	return LinRange(d.istart,d.istop,length(d.dat))
end

function val(d::Data1D)
    return d.dat
end


function val(d::Data1D,ind)
    i1=(ind-d.istart)*length(d.dat)/(d.istop-d.istart)
    p1=floor(Int64,i1)+1
    ic=i1-p1+1
    return ((1-ic)*d.dat[p1]+ic*d.dat[p1+1])
end

function cut(d::Data1D,i1,i2)
    p1=ind2pos(d,i1)
    p2=ind2pos(d,i2)
    return Data1D(d.dat[p1:p2],pos2ind(d,p1),pos2ind(d,p2))
end



"""
`FourierTransform(d::Data1D)` performs an FFT assuming that the
input data set `d` is a free induction decay, and returns
a `DataSet` object with the resulting complex data.
"""
function FourierTransform(fid::Data1D;
        CTR=0.0,
        LB=0.0,
        SI=length(fid.dat),
        PPM=1,FFTSHIFT=true)

    apod=[exp(-t*LB) for t=ind(fid)];
    s=zeros(Complex{Float64},SI);
    s[1:length(apod)]=apod.*fid.dat;
    s[1]*=0.5;
    SW=length(fid.dat)/(fid.istop-fid.istart);
    spect=fft(s);
    if(FFTSHIFT)
        spect=fftshift(spect)
    end;
    return Data1D(spect,CTR-SW/2/PPM,CTR+SW/2/PPM)
end

"""
`PhaseCorrect(spect::Data1D;Ph0=0.0,Ph1=0.0,Pivot=0.0)`
returns a phase corrected `DataSet`, using the values indicated.
The parameters are given in rad, k, and rad/k, where k are the units
of the spectral domain (typically, Hz or ppm).
"""
function PhaseCorrect(spect::Data1D;Ph0=0.0,Ph1=0.0,Pivot=0.0)
    out=deepcopy(spect);
    out.dat = spect.dat .* [exp(im*(Ph0+(k-Pivot)*Ph1)) for k=ind(spect)]
    return out;
end

polyfit(x::Vector, y::Vector, deg::Int) = pinv(collect(v ^ p for v in x, p in 0:deg)) * y

function horner(x::Vector, c::Vector)
	r=ones(length(x))*c[end]
    for k=(length(c)-1):-1:1
        r.*=x
        r.+=c[k]
    end
    return r
end

"""
`BaseLineCorrect(spect::Data1D;regions=128,kfactor=5)` corrects the
base line of `spect` using the algorithm
by S. Golotvin, A. Williams, J. Magn. Reson. 146 (2000) 122-125.

`BaseLineCorrect(spect::Data1D,pts::Array{Float64,1})`
uses the points in `pts` to determine the baseline. This should be a list
of positions in the spectrum that are known not to contain signal.

Note that if the input data is complex, the imaginary part is ignored,
and the result is always real.
"""
function BaseLineCorrect(spect::Data1D;regions=128,kfactor=5,wdw=32,order=5,verbose=false)

    out=Data1D(real.(spect.dat),spect.istart,spect.istop);
    n=length(out.dat)
    lreg=floor(Int64,n/regions)

    # find smallest standard deviation in regions
    mindev=minimum(abs.([std(out.dat[k:k+lreg]) for k=1:lreg:(n-lreg)]))

    # find base line points
    select=[abs(std(out.dat[(k-wdw):(k+wdw)])) < kfactor*mindev for k=(wdw+1):(n-wdw)]

    pts=((wdw+1):(n-wdw))[select];
    verbose && println("Fitting $(length(pts)) Points.")

    # polynomial fitting
    coeffs=polyfit(pts/n,out.dat[pts],order)
    verbose && println("Coefficients: $coeffs")

    res=sum(abs.(out.dat[pts]-horner(pts/n,coeffs)))
    verbose && println("Residual: $res")

    out.dat .-= horner(collect(1:n)/n,coeffs)

    return (out)

end

function BaseLineCorrect(spect::Data1D,ind::Array{Float64,1};order=5)
    out=Data1D(real.(spect.dat),spect.istart,spect.istop);
    n=length(out.dat)
    pts=[ind2pos(out,k) for k in ind]
    coeffs=polyfit(pts/n,out.dat[pts],order)
    out.dat .-= horner(collect(1:n)/n,coeffs)
    return(out)
end


"""
`integral(spect::Data1D)` computes the numerical integral of `spect` by
the trapezoid rule.
"""
function integral(spect::Data1D)
    s=sum(spect.dat[2:(end-1)])+0.5*spect.dat[1]+0.5*spect.dat[end]
    return s*(spect.istop-spect.istart)/length(spect.dat)
end


"""
`derivative(spect::Data1D)` computes the derivative of `spect` and returns
it as another Data1D object
"""
function derivative(spect::Data1D)
    s=spect.dat
    inc=(spect.istop-spect.istart)/length(spect.dat)
#    d[1:(end-1)]=(spect.dat[2:end]-spect.dat[1:(end-1)] )/inc
#    d[end]=d[end-1]
    d = 1.0/12*(8*[s[2:end];0]-8*[0;s[1:(end-1)]] + [s[3:end];0;0] - [0;0;s[1:(end-2)]] )/inc
    return Data1D(d,spect.istart,spect.istop)
end


function Plot(s::Data1D;opts...)
    return Plot(ind(s),real(val(s)),style=Dict(["stroke-width"=>"1"]),Reverse=[true,false];opts...)
end


function set!(spect::Data1D,start,stop,value=0)
    rge=ind2pos(spect,start):ind2pos(spect,stop)
    spect.dat[rge]=value
end

import Base.*, Base./, Base.+,Base.-,Base.abs

function *(d::Data1D,n::Number)
    c=deepcopy(d);
    c.dat .*= n;
    return c
end

*(n::Number,d::Data1D)=d*n

function *(d1::Data1D,d2::Data1D)
	ind(d1) == ind(d2) || error("Data1D objects incompatible for multiplication")
	c=deepcopy(d1);
	c.dat = d1.dat .* d2.dat
	return c
end

function /(d::Data1D,n::Number)
    c=deepcopy(d);
    c.dat ./= n;
    return c
end

function /(d1::Data1D,d2::Data1D)
	ind(d1) == ind(d2) || error("Data1D objects incompatible for division")
	c=deepcopy(d1);
	c.dat = d1.dat ./ d2.dat
	return c
end

function +(d::Data1D,n::Number)
    c=deepcopy(d);
    c.dat .+= n;
    return c
end

+(n::Number,d::Data1D)=d+n

function +(d1::Data1D,d2::Data1D)
	ind(d1) == ind(d2) || error("Data1D objects incompatible for addition")
	c=deepcopy(d1);
	c.dat = d1.dat .+ d2.dat
	return c
end

function -(d::Data1D,n::Number)
    c=deepcopy(d);
    c.dat .-= n;
    return c
end

function -(d1::Data1D,d2::Data1D)
	ind(d1) == ind(d2) || error("Data1D objects incompatible for subtraction")
	c=deepcopy(d1);
	c.dat = d1.dat .- d2.dat
	return c
end

function abs(d::Data1D)
    c=Data1D(abs.(d.dat),d.istart,d.istop)
    return c
end

function plot(d::Data1D;opts...)
	return Plots.plot(real.(ind(d)),real.(val(d));opts...)
end

import Plots.RecipesBase.plot!

function plot!(d::Data1D;opts...)
	return Plots.plot!(real.(ind(d)),real.(val(d));opts...)
end

import Base.length

function length(d::Data1D)
	return length(d.dat)
end

function shift(d::Data1D,δ::Number)
	c=deepcopy(d);
	c.istart += δ
	c.istop += δ
	return c
end

"""
`function hard_shift(d::Data1D,δ::Number)`

shifts the data in `d` by `δ`. The data is re-interpolated as needed; the 
span of the index of the new data set is reduced by |`δ`|.
"""

function hard_shift(d::Data1D,δ::Number)
    # compute integer offset and fractional contribution from left
    inc=(d.istop-d.istart)/length(d.dat)
    n = -convert(Integer,div(δ,inc))
    f = mod(δ,inc)/inc
    if n>=0
        ndata=(1-f)*d.dat[(1+n):(end-1)]+f*d.dat[(2+n):end]
        nstart=d.istart
        nstop=d.istop+δ
    else
        ndata=(1-f)*d.dat[1:(end+n-1)]+f*d.dat[2:(end+n)]
        nstart=d.istart+δ
        nstop=d.istop
    end
    return Data1D(ndata,nstart,nstop)
end

#end
