module AutoPhase

using NMR.DataSet

export entropy

ent(x) = x*log(x)

"""
`entropy(spect::Data1D)` computes the entropy in the real part of an
NMR spectrum as defined by Chen et al. in
*Journal of Magnetic Resonance* **158** (2002) 164–168.
This quantity can be optimised with respect to zero- and first-order
phase correction for automatic (unsupervised) phase correction.
"""
function entropy(spect::Data1D)
    return sum(ent.(abs.(real.(spect.dat))))/sum(abs.(real.(spect.dat)))
end


end
