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
export medianBaseline

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


"""
    Data1D(f<:AbstractArray,start<:Real,stop<:Real) 

returns a data structure to hold a 1D NMR dataset either in the time or the frequency domain. `Data1D` objects
store an *index* (range of time points or frequency points) along with the corresponding y-axis *data*.
`Data1D` objects can be added and subtracted from one another, as long as their  indices match.
"""
mutable struct Data1D{Tdata,Tindex}
   dat::Tdata
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
    i1=(ind-d.istart)*(length(d.dat))/(d.istop-d.istart)
    p1=floor(Int64,i1+1)
    ic=i1-p1+1
    p2 = p1+1
    if p2>length(d.dat) # this ugly hack is needed in case p1 happens to be the
                        # last point in the vector 
        p2=p1
    end
    return ((1.0-ic)*d.dat[p1]+ic*d.dat[p2])
end

function cut(d::Data1D,i1,i2)
    p1=ind2pos(d,i1)
    p2=ind2pos(d,i2)
    return Data1D(d.dat[p1:p2],pos2ind(d,p1),pos2ind(d,p2))
end

resolution(d::Data1D) = (d.istop-d.istart)/length(d.dat)


"""
    FourierTransform(d::Data1D) 
    
performs an FFT assuming that the
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

@doc raw"""
## Automatic Baseline Correction
`BaseLineCorrect(spect::Data1D;regions=128,kfactor=5,wdw=32)` corrects the
base line of `spect` using the algorithm
by S. Golotvin, A. Williams, J. Magn. Reson. 146 (2000) 122-125.



### Optional Parameters:
- `regions`: the number of regions the spectrum is divided into to estimate
   the noise variance. This should be chosen large enough so that at least one
   of these regions is free from signals, and shows a flat baseline. The 
   variance ``\\sigma^2_{min}`` in this region is then used to distinguish signals and noise.

- `wdw`: the window size around each data point that is used to decide whether the
  data point is part of a signal, or should be counted as base line. 

- `kfactor`: a point will be counted as part of the baseline if and only if the local
  variance in its window is less than `kfactor`*``\\sigma^2_{min}``.

- `output`: defaults to `:spectrum`. `output=:points` will return an array of points,
  which indicate baseline positions (locations free from peaks) as determined
  by the algorithm. This is useful in situations where a large number of similar
  spectra need to be processed.

### Notes
- if the input data is complex, the imaginary part is ignored,
  and the result is always real.

- the first and last point of `spect` are always used as part of the baseline. 
  The first and last points of the output are therefore guaranteed to be zero. This
  helps to keep the interpolation numerically stable.


## Baseline Correction With Known Locations

`BaseLineCorrect(spect::Data1D,pts::Array{Float64,1})`
uses the points in `pts` to determine the baseline. This should be a list
of positions in the spectrum that are known not to contain signal; it can be obtained
by using the option `output=:points` in the first method of `BaseLineCorrect()`. 

"""
function BaseLineCorrect(spect::Data1D;regions=128,kfactor=5,lw=0.0,wdw=32,order=5,verbose=false,output=:spectrum)

    out=Data1D(real.(spect.dat),spect.istart,spect.istop);
    n=length(out.dat)
    lreg=floor(Int64,n/regions)
    
    # compute wdw from lw parameter if set
    if lw>0.0
        wdw = round(Int32,4*lw/resolution(spect))
    end

    # find smallest standard deviation in regions
    mindev=minimum(abs.([std(out.dat[k:k+lreg]) for k=1:lreg:(n-lreg)]))

    # find base line points
    select=[abs(std(out.dat[(k-wdw):(k+wdw)])) < kfactor*mindev for k=(wdw+1):(n-wdw)]

    pts=((wdw+1):(n-wdw))[select];
    verbose && println("Fitting $(length(pts)) Points.")
    push!(pts,1);
    push!(pts,n);
    
    # polynomial fitting
    coeffs=polyfit(pts/n,out.dat[pts],order)
    verbose && println("Coefficients: $coeffs")

    res=sum(abs.(out.dat[pts]-horner(pts/n,coeffs)))
    verbose && println("Residual: $res")

    out.dat .-= horner(collect(1:n)/n,coeffs)
    
    if output==:points
        return([pos2ind(spect,p) for p in pts])
    else
        return (out)
    end
end

function BaseLineCorrect(spect::Data1D,ind::Array{Float64,1};order=5)
    out=Data1D(real.(spect.dat),spect.istart,spect.istop);
    n=length(out.dat)
    pts=[ind2pos(out,k) for k in ind]
    coeffs=polyfit(pts/n,out.dat[pts],order)
    out.dat .-= horner(collect(1:n)/n,coeffs)
    return(out)
end


wrap(n,l)=[mod(k,l) for k in n]

function extrema(X::Array{Float64,1})
    Y=Array{Float64,1}()
    for k=2:(length(X)-1)
        if (X[k]>X[k-1] && X[k]>X[k+1]) || X[k]<X[k-1] && X[k]<X[k+1]
            push!(Y,X[k])
        end
    end
    return Y
end

@doc raw"""
    function medianBaseline(s::NMR.Data1D;wdw=length(s)>>6)

Compute baseline for the real part of `s` by the algorithm of M. S. Friedrichs,
*Journal of Biomolecular NMR*,  **5** (1995) 147  153.
The window is automatically calculated as ``\Delta\omega/2^6``, unless
overridden by giving a value to the optional key `wdw`.

"""
function medianBaseline(s::NMR.Data1D{Tdata,Tindex};wdw=length(s)>>6) where { Tdata<:Real, Tindex}
    r=s.dat
    L=length(r)
    b=zeros(L)
    c=zeros(L)
    for k=1:L
        b[k]=median(extrema(r[ wrap(k.+(-wdw:wdw),1:L)]))
    end
    gauss=exp.(-((wdw:wdw)./wdw).^2)
    gauss=gauss/sum(gauss)
    for k=1:length(r)
        c[k]=sum(gauss.*b[wrap(k.+(-wdw:wdw),1:L)])/(2*wdw+1)
    end
    return NMR.Data1D(c,s.istart,s.istop)
end

function medianBaseline(s;opts...)
    return medianBaseline(real(s))+im*medianBaseline(imag(s);opts...)
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
    c = n * d.dat;
    return Data1D(c,d.istart,d.istop)
end

*(n::Number,d::Data1D)=d*n

function *(d1::Data1D,d2::Data1D)
	ind(d1) == ind(d2) || error("Data1D objects incompatible for multiplication")
	return Data1D(d1.dat.*d2.dat,d1.istart,d1.istop)
end

function /(d::Data1D,n::Number)
    return Data1D(d.dat./n,d.istart,d.istop)
end

function /(d1::Data1D,d2::Data1D)
	ind(d1) == ind(d2) || error("Data1D objects incompatible for division")
	return Data1D(d1.dat./d2.dat,d1.istart,d1.istop)
end

function +(d::Data1D,n::Number)
    return Data1D(n.+d.dat,d.istart,d.istop)
end

+(n::Number,d::Data1D)=d+n

function +(d1::Data1D,d2::Data1D)
	ind(d1) == ind(d2) || error("Data1D objects incompatible for addition")
	return Data1D(d1.dat.+d2.dat,d1.istart,d1.istop)
end

function -(d::Data1D,n::Number)
    return Data1D(-n.+d.dat,d.istart,d.istop)
end

function -(d1::Data1D,d2::Data1D)
	ind(d1) == ind(d2) || error("Data1D objects incompatible for subtraction")
	return Data1D(d1.dat.-d2.dat,d1.istart,d1.istop)
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

shifts the data in `d` by `δ`. The data is re-interpolated as needed. The
resulting `Data1D` object is guaranteed to have the same index as `d`, such that
they are compatible for the purposes of addition etc. Data shifted out of the 
window is lost; data points shifted into view are set to zero.
"""
function hard_shift(d::Data1D,δ::Number)
    # compute integer offset and fractional contribution from left
    inc=(d.istop-d.istart)/(length(d.dat)-1)
    n = convert(Integer,div(δ,inc))
    f = mod(δ,inc)/inc
    ndata=zero(d.dat)
    nstart=d.istart
    nstop=d.istop
    
    if n<0
        ndata[1:(end+n-1)]=(1-f)*d.dat[(1-n):(end-1)]+f*d.dat[(2-n):end]
    else
        ndata[n+1:end-1]=(1-f)*d.dat[1:end-n-1]+f*d.dat[2:end-n]
    end
    
    return Data1D(ndata,nstart,nstop)
end

import Base.real, Base.imag, Base.abs, Base.angle

function Base.real(s::NMR.Data1D)
    return(NMR.Data1D(real.(s.dat),s.istart,s.istop))
end

# function Base.abs(s::NMR.Data1D)
#     return(NMR.Data1D(abs.(s.dat),s.istart,s.istop))
# end

function Base.imag(s::NMR.Data1D)
    return(NMR.Data1D(imag.(s.dat),s.istart,s.istop))
end

function Base.angle(s::NMR.Data1D)
    return(NMR.Data1D(angle.(s.dat),s.istart,s.istop))
end

#end
