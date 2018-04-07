module AutoPhase

using NMR.DataSet
using Optim

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
    return sum([k<0.0?k*k:0.0 for k in x])
end

# this is the minimisation target for automatic phase correction
function goalfun(x,spect,γ=1.e-5)
    piv=(spect.istart+spect.istop)/2
    c=PhaseCorrect(spect,Ph0=x[1],Ph1=x[2],Pivot=piv);
    return entropy(derivative(c))+γ*penalty(real.(val(c)));
end

function AutoPhaseCorrect(spect::Data1D;verbose=false)
  piv=(spect.istart+spect.istop)/2
  result=optimize(x->goalfun(x,spect),[π/2,-0.1],NelderMead());
  if verbose print(result) end;
  pc=Optim.minimizer(result);
  return PhaseCorrect(spect,Ph0=pc[1],Ph1=pc[2],Pivot=piv);
end



end
