# User Manual

## Getting started

`NMR.jl` is a Julia package. It is installed using the Julia
package management system:
```@example
using Pkg
Pkg.add(url="https://github.com/marcel-utz/NMR.jl")
```

Once the install finishes, you are ready to go.

## Loading NMR Data
To start processing, you need to load some NMR data.
There are commands to load raw data from various spectrometer
vendors. Typically, the data is loaded raw. This means that
only the time-domain spectral data is returned, and you
have to manually convert it to a data object that can be further
processed.

Here is an example:
```@example brukerExpl
ENV["GKSwstype"] = "100" # hide
import NMR
import Plots # hide

f=NMR.readBrukerFID("../../test/data/10/fid")
````
`f`returns an array with the complex data points contained in the 
`fid` file. To convert this to useable time-domain data, 
you need to convert it into a `Data1D` object. 

```@docs
    Data1D
```

We assume
that the data represents a free induction decay with a total
duration of 1 second.

```@example brukerExpl
d=NMR.Data1D(f,0.0,1.0)
NMR.plot(real(d))
Plots.savefig("plot.svg"); nothing # hide
```
![](plot.svg)

```@example brukerExpl
spect=NMR.FourierTransform(d,PPM=700, CTR=4.76)
NMR.plot(real(spect),xaxis=:flip)
Plots.savefig("plot-spec.svg"); nothing # hide
```
![](plot-spec.svg)