module DataSet

export Data1D,ind2pos,pos2ind,val,ind,cut
export FourierTransform,PhaseCorrect,BaseLineCorrect

# Data structures for storing digitised data sets
# Definitions:
#
# - Index: the horizontal index of a data point
# - Position: the numerical position (integer) of the data point in
#           the internal storage array

type Data1D{Tdata,Tindex}
   dat::Array{Tdata,1}
   istart::Tindex
   istop::Tindex
end

function ind2pos(d::Data1D,ind)
    return round(Int64,(ind-d.istart)*length(d.dat)/(d.istop-d.istart))+1
end

function pos2ind(d::Data1D,pos::Integer)
   return d.istart+(pos-1)/length(d.dat)*(d.istop-d.istart)
end

function ind(d::Data1D)
   return d.istart+((1:length(d.dat))-1)/(length(d.dat)-1)*(d.istop-d.istart)
end

function val(d::Data1D)
    return d.dat
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
        PPM=1)

    apod=[exp(-t*LB) for t=ind(fid)];
    s=zeros(Complex{Float64},SI);
    s[1:length(apod)]=apod.*fid.dat;
    s[1]*=0.5;
    SW=length(fid.dat)/(fid.istop-fid.istart);
    spect=fftshift(fft(s));
    return Data1D(spect,CTR-SW/2/PPM,CTR+SW/2/PPM)
end

"""
`PhaseCorrect(spect::Data1D;Ph0=0.0,Ph1=0.0,Pivot=0.0)`
returns a phase corrected `DataSet`, using the values indicated.
The parameters are given in rad, k, and rad/k, where k are the units
of the spectral domain (typically, Hz or ppm).
"""
function PhaseCorrect(spect::Data1D;Ph0=0.0,Ph1=0.0,Pivot=0.0)
    out=spect;
    out.dat = spect.dat .* [exp(im*(Ph0+(k-Pivot)*Ph1)) for k=ind(spect)]
    return out;
end

polyfit(x::Vector, y::Vector, deg::Int) = collect(v ^ p for v in x, p in 0:deg) \ y

function horner(v::Vector, c::Vector)
    r=zeros(v)
    for k=length(c):-1:2
      r += c[k]
      r .*= v
    end
    r += c[1]
    return r
end

"""
`BaseLineCorrect(spect::Data1D;regions=128,kfactor=5)` corrects the
base line of `spect` using the algorithm
by S. Golotvin, A. Williams, J. Magn. Reson. 146 (2000) 122-125.
"""
function BaseLineCorrect(spect::Data1D;regions=128,kfactor=5,wdw=32,order=5)
    out=spect;
    n=length(out.dat)
    lreg=floor(Int64,n/regions)

    # find smallest standard deviation in regions
    mindev=minimum(abs.([std(out.dat[k:k+lreg]) for k=1:lreg:(n-lreg)]))

    # find base line points
    select=[abs(std(out.dat[k:(k+wdw)])) < kfactor*mindev for k=1:(n-wdw)]

    pts=(1:(n-wdw))[select];

    # polynomial fitting
    coeffs=polyfit(pts,out.dat[pts],order)

    out.dat -= horner(collect(1:n),coeffs)

    return out

end

end
