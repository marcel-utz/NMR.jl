# Home

## Contents
```@contents
```

## Example documentation
```@docs
    NMR
    FourierTransform
```

## Here is a code example
```@example
ENV["GKSwstype"] = "100" # hide
import NMR
import Plots # hide

f=NMR.readBrukerFID("../../test/data/10/fid")
d=NMR.Data1D(f,0.0,1.0)
NMR.plot(real(d))
Plots.savefig("plot.svg"); nothing # hide
```
![](plot.svg)

```@example
import NMR # hide
import Plots # hide
ENV["GKSwstype"] = "100" # hide
f=NMR.readBrukerFID("../../test/data/10/fid")
d=NMR.Data1D(f,0.0,1.0)
spect=NMR.FourierTransform(d,PPM=700, CTR=4.76)
NMR.plot(real(spect),xaxis=:flip)
Plots.savefig("plot-spec.svg"); nothing # hide
```
![](plot-spec.svg)