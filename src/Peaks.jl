export peaks, PeakStruct

mutable struct PeakStruct
    positions::Array{Float64,1}
    heights::Array{Float64,1}
    widths::Array{Float64,1}
    intensities::Array{Float64,1}
    deconvolution::Data1D
end

"""
`smooth(d::Data1D,n=10)`
returns a smoothed version of the input data set `d`, using
convolution with a Gaussian peak.
"""
function smooth(d::Data1D,n=10)
      l=length(d.dat);
      kern=[exp(-5*j^2/n^2) for j=-n:n]
      kern /= sum(kern)
      return(Data1D( [sum([d.dat[k+j]*kern[j+n+1] for j=-n:n]) for k=(n+1):(l-n)],pos2ind(d,n+1),pos2ind(d,l-n-1)));
end

lorentzian(x0::Float64,σ::Float64,x::Float64) = 1.0/π*sqrt(σ)/(1.0+σ*(x-x0)^2)

"""
`peaks(dinput::Data1D;threshold=1,athresh=1,regions=128)`
identifies the peaks in `dinput`. Peak positions are located roughly
by searching for the downward zero crossings of the first derivative.
The data points surrounding each peak are then used for a quadratic regression,
from which the position, curvature, and height of the maximum are interpolated.
The signal `dinput` is then decomposed into a linear superposition of Lorentzian
peaks at the determined positions and with the determined curvatures.

`peaks()` returns a data structure
```
struct PeakStruct
    positions::Array{Float64,1}
    heights::Array{Float64,1}
    widths::Array{Float64,1}
    intensities::Array{Float64,1}
    deconvolution::Data1D
  end
```

where `positions` and `heights` are the interpolated peak positions and heights,
respectively; `widths` are the peak curvatures, which correspond
to 1/FWHM line widths of the Lorentzian peaks. `intensities` are the
peak intensities resulting from deconvolution. These values are
proportional to the underlying concentrations in an NMR spectrum.
Finally, `deconvolution` contains the deconvoluted approximation to the
original input signal.
"""
function peaks(dinput::Data1D;threshold=1,athresh=1,regions=128)
    d=smooth(dinput)
    dr=derivative(d);
    n=length(d.dat)
    δ=(d.istop-d.istart)/n
    lreg=floor(Int64,n/regions)

    # find smallest standard deviation in regions
    mindev=minimum(abs.([std(dinput.dat[k:(k+lreg)]) for k=1:lreg:(n-lreg)]))
    mindevd=minimum(abs.([std(dr.dat[k:(k+lreg)]) for k=1:lreg:(n-lreg)]))

    # locate indices of peak positions roughly on the basis of the derivative
    # and original data set
    ind=10:1:(length(dr.dat)-10);
    peaks=filter(k->dr.dat[k-1]>0 && dr.dat[k]<0 && dr.dat[k-1]-dr.dat[k] > athresh*mindevd && dinput.dat[k]>mindev*threshold, ind);

    # For each peak, perform a second order polynomial fit to the
    # peak tip. More accurate locations and equivalent peak width
    # information is extracted from the resulting coefficients

    dpeaks=[ind2pos(dinput,pos2ind(dr,k)) for k in peaks] # a hack to eliminate aligment differences
    coeffs = [NMR.polyfit(δ*collect(-9:9),dinput.dat[(p-9):(p+9)],2)
                    for p in dpeaks]

    # compute updated peak positions, intensities, and lorentzian equivalent
    # line widths from the expansion coefficients

    σ=[abs(C[3]) for C in coeffs] ;    # inverse line widths
    cpks=[pos2ind(dr,k) for k in peaks]+[-0.5*C[2]/C[3] for C in coeffs]; # corrected peak positions
    hpks=[C[1]-0.25*C[2]*C[2]/C[3] for C in coeffs]; # corrected peak heights

    # devonvolute the spectrum with Lorentzian peaks of the determined
    # positions and widths

    B=collect(lorentzian(cpks[k],0.66*σ[k],x) for x in NMR.ind(dinput), k in 1:length(cpks))  ;
    c=pinv(B)* NMR.val(dinput);

    return(PeakStruct(cpks,hpks,σ,c,Data1D(B*c,dinput.istart,dinput.istop)))
end
