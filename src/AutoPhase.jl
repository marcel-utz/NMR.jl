#module AutoPhase

#using NMR.DataSet
import Optim

export entropy, AutoPhaseCorrect

ent(x) = x*log(x)

"""
`entropy(spect::Data1D)` computes the entropy in the real part of an
NMR spectrum as defined by Chen et al. in
*Journal of Magnetic Resonance* **158** (2002) 164–168.
This quantity can be optimised with respect to zero- and first-order
phase correction for automatic (unsupervised) phase correction.
"""
function entropy(spect::Data1D)
    h=abs.(real.(val(spect)))
    h/=sum(h)
    return -sum(ent.(h))
end

# penalty(x) computes the sum of squares of all negative points in x
function penalty(x)
    x /= sum(abs.(x))
    return sum([k<0.0 ? k*k : 0.0 for k in x])
end

# this is the minimisation target for automatic phase correction
function goalfun(x,spect,γ=1.e-5)
    piv=(spect.istart+spect.istop)/2
    c=PhaseCorrect(spect,Ph0=x[1],Ph1=x[2],Pivot=piv);
    return entropy(derivative(c))+γ*penalty(real.(val(c)));
end

"""
`AutoPhaseCorrect(spect::Data1D;verbose=false)`
performs zero- and first-order phase correction of `spect`
using a minimum entropy algorithm  by Chen et al. in
*Journal of Magnetic Resonance* **158** (2002) 164–168.
It uses the Nelder-Mead algorithm, as implemented in the `Optim.jl`
package.
"""
function AutoPhaseCorrect(spect::Data1D;verbose=false,γ=1.e-5)
  piv=(spect.istart+spect.istop)/2
  result=Optim.optimize(x->goalfun(x,spect,γ),[π/2,-0.1],Optim.NelderMead());
  if verbose print(result) end;
  pc=Optim.minimizer(result);
  return PhaseCorrect(spect,Ph0=pc[1],Ph1=pc[2],Pivot=piv);
end



#end
