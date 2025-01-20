#module AutoPhase

#using NMR.DataSet
import Optim

export entropy, AutoPhaseCorrect, AutoPhaseCorrectChen

ent(x) = x*log(x)

"""
`entropy(spect::Data1D)` computes the entropy of the first derivative
in the real part of an
NMR spectrum as defined by Chen et al. in
*Journal of Magnetic Resonance* **158** (2002) 164–168.
This quantity can be optimised with respect to zero- and first-order
phase correction for automatic (unsupervised) phase correction.
"""
function entropy(spect::Data1D)
    h=abs.(real.(val(derivative(spect))))
    h/=sum(h)
    return -sum(ent.(h))
end

# penalty(x) computes the sum of squares of all negative points in x
function penalty(x)
    x /= sum(abs.(x))
    return sum([k<0.0 ? k*k : 0.0 for k in x])
end

# this is the minimisation target for automatic phase correction
function goalfun(x,spect,γ=0.0e-5)
    piv=0.5*(spect.istop+spect.istart)
    c=PhaseCorrect(spect,Ph0=x[1]*10,Ph1=x[2],Pivot=piv);
    return entropy(c)+γ*penalty(real.(val(c)));
end

"""
`AutoPhaseCorrectChen(spect::Data1D;verbose=false)`
performs zero- and first-order phase correction of `spect`
using a minimum entropy algorithm  by Chen et al. in
*Journal of Magnetic Resonance* **158** (2002) 164–168.
It uses the Nelder-Mead algorithm, as implemented in the `Optim.jl`
package.
"""
function AutoPhaseCorrectChen(spect::Data1D;verbose=false,γ=0.0)
  piv=0.5*(spect.istop+spect.istart)
  result=Optim.optimize(x->goalfun(x,spect,γ),[0.1,0.2],Optim.NelderMead(),Optim.Options(show_trace=false,g_tol=1.0e-9));
  if verbose print(result) end;
  pc=Optim.minimizer(result);
  scorr = PhaseCorrect(spect,Ph0=pc[1]*10,Ph1=pc[2],Pivot=piv);
  if (real(integral(scorr))<0.0) scorr = -1.0*scorr end;
  return scorr ;
end


function unwrap(a;p=π,tol=0.25)
    l=length(a)
    w=zeros(Float64,l)
    offset=0
    for k=1:(l-1)
        if( (p-a[k])<p*tol && (a[k+1]+p)<p*tol)
            offset+=2*p;
        elseif( (a[k]+p)<p*tol && (p-a[k+1])<p*tol)
            offset-=2*p;
        end
        w[k+1]=offset;
    end

    return a+w
end

"""
    AutoPhaseCorrect(spect::Data1D; exclude=false,threshold=10)

performs zero- and first-order phase correction of `spect`
using a the algorithm by van der Waals and Geerens,
*Journal of Magnetic Resonance* **86** (1990) 127-154.
It works by recognising the peaks in the magnitude mode spectrum
and then representing the phase at each peak location
through a linear regression.

Passing a pair of indices to `exclude` disregards any peaks 
between these two boundaries. This is useful in spectra that
contain a solvent artefact that may not be phased correctly.

`threshold` defines a threshold for peak detection (cf. `NMR.peaks()`).
"""
function AutoPhaseCorrect(s::NMR.Data1D;exclude=false,threshold=10)
    # --- Step 1: basic phase adjustment
    phase0=angle(sum(s.dat))
    s1=NMR.PhaseCorrect(s,Ph0=-phase0)

    # --- Step 2: recognise peaks in the absolute mode spectrum
    s2=abs(s1);
    pks=NMR.peaks(s2,threshold=threshold);

    # --- Step 3: obtain list of peak phases
    ppos=[NMR.ind2pos(s1,x) for x in pks.positions]

    if exclude!=false
        i=NMR.ind2pos(s,exclude[1]);
        j=NMR.ind2pos(s,exclude[2]);
        ppos=filter(x->(x<i||x>j),ppos)
    end
    pind = [NMR.pos2ind(s,k) for k in ppos]
    pphase=[angle(s1.dat[p]) for p in ppos]
    α=NMR.polyfit(pind,unwrap(pphase),1)
    # print(α,"\n")
    s3=NMR.PhaseCorrect(s1,Ph0=-α[1],Ph1=-α[2])
    return s3
end



#end
