# User Manual

## Getting started

`NMR.jl` is a Julia package. It is installed using the Julia
package management system:
```@julia
import Pkg
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
```
`f`returns an array with the complex data points contained in the 
`fid` file. To convert this to useable time-domain data, 
you need to convert it into a `Data1D` object. 

Bruker NMR systems store the acquisition parameters in a separate file,
which we can read into a Julia dictionary as follows:
```@example brukerExpl
acqus=NMR.readBrukerParameterFile("../../test/data/10/acqus")
```
To convert the raw fid data into a `Data1D` object which we can process, 
we need some of these parameters. The time step between subsequent points
in the FID is given by the inverse of the spectral width. Moreover, Bruker
FIDs actually begin *before* ``t=0``. We have to remove these points:

```@example brukerExpl
dwellTime=1/acqus["SW_h"]
f=f[74:end]
d=NMR.Data1D(f,0.0,length(f)*dwellTime)
NMR.plot(real(d))
Plots.savefig("plot-fid.svg"); nothing # hide
```
![](plot-fid.svg)

To convert the time-domain data into a spectrum, we use `NMR.FourierTransform`. By default,
FourierTransform interprets the time domain in the `Data1D` object in seconds, and produces
another `Data1D` object with a horizontal axis in Hz:
```@example brukerExpl
spect=NMR.FourierTransform(d)
NMR.plot(real(spect),xaxis=:flip,xlabel="frequency [Hz]")
Plots.savefig("plot-spec.svg"); nothing # hide
```
![](plot-spec.svg)

To obtain a spectrum with a horizontal axis in ppm chemical shift, we have 
to indicate the conversion. This is done by giving the number of Hz per ppm 
as a parameter. In our case, the spectrum was acquired on a 700 MHz spectrometer.
One ppm therefore corresponds to 700 Hz. The precise factor
is contained in the Bruker parameter `"SFO1"`. Also, we can calibrate the horizontal axis
by indicating the chemical shift at the centre of the spectrum:
```@example brukerExpl
spect=NMR.FourierTransform(d,PPM=acqus["SFO1"],CTR=4.76)
NMR.plot(real(spect),xaxis=:flip,xlabel="1H Chem. Shift [ppm]")
Plots.savefig("plot-specHz.svg"); nothing # hide
```
![](plot-specHz.svg)

## Phase Correction
The above spectrum still shows artefacts. To clean it up, we need
to correct the phase. This can either be done manually,
by supplying `PH0`and `PH1` arguments to `FourierTransform`, or
we can resort to automatic phase correction:

```@example brukerExpl
spect = NMR.AutoPhaseCorrectChen(spect)
NMR.plot(real(spect),xaxis=:flip,xlabel="1H Chem. Shift [ppm]")
Plots.savefig("plot-specHz-ap.svg"); nothing # hide
```
![](plot-specHz-ap.svg)

## Cutting Regions
The above spectrum covers a wide range without any signals, which
is of no interest. We can chop out the central, important part:
```@example brukerExpl
s2=NMR.cut(spect,-0.5,7.0)
NMR.plot(real(s2),xaxis=:flip,xlabel="1H Chem. Shift [ppm]")
Plots.savefig("plot-S2.svg"); nothing # hide
```
![](plot-S2.svg)