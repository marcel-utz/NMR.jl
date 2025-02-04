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

The `NMR.jl` comes with example data that can be used to play
with its functions. The data can be found through a dictionary in the module
`Examples`:
```@example brukerExpl
import NMR
import NMR.Examples

data=Examples.Data["HCC cell culture media spectra"]
```
### Bruker data
We can read raw data from an "fid" file in this set as follows:
```@example brukerExpl
ENV["GKSwstype"] = "100" # hide
import Plots # hide
f=NMR.readBrukerFID(data["files"][1]*"/fid")
```
`f`returns an array with the complex data points contained in the 
`fid` file. To convert this to useable time-domain data, 
you need to convert it into a `Data1D` object. 

Bruker NMR systems store the acquisition parameters in a separate file,
which we can read into a Julia dictionary as follows:
```@example brukerExpl
acqus=NMR.readBrukerParameterFile(data["files"][1]*"/acqus")
```
To convert the raw fid data into a `Data1D` object which we can process, 
we need some of these parameters. The time step between subsequent points
in the FID is given by the inverse of the spectral width. Moreover, Bruker
FIDs actually begin *before* ``t=0`` due to the 
group delay of the digital filter. The number of 
points to remove is found in the `acqus` parameter
`"GRPDLY"`:

```@example brukerExpl
dwellTime=1/acqus["SW_h"]
f=f[acqus["GRPDLY"]:end]
d=NMR.Data1D(f,0.0,length(f)*dwellTime)
NMR.plot(real(d))
Plots.savefig("plot-fid.svg"); nothing # hide
```
![](plot-fid.svg)

To convert the time-domain data into a spectrum, we use `NMR.FourierTransform`. By default,
FourierTransform interprets the time domain in the `Data1D` object in seconds, and produces
another `Data1D` object with a horizontal axis in Hz:
```@example brukerExpl
import Plots
spect=NMR.FourierTransform(d)
Plots.plot(real(spect),xaxis=:flip,xlabel="frequency [Hz]")
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
import Plots
spect=NMR.FourierTransform(d,PPM=acqus["SFO1"],CTR=4.76)
Plots.plot(real(spect),xaxis=:flip,xlabel="1H Chem. Shift [ppm]")
Plots.savefig("plot-specHz.svg"); nothing # hide
```
![](plot-specHz.svg)

### JEOL data
We can also import data from a JEOL spectrometer. In this case, the parameters and the
data are combined in the same (binary) file. JEOL `.jdf` files have a rather elaborate
structure. Fortunately, there is a function that parses this, and returns the parameters as
a dictionary, as well as the raw data:
```@example JEOLExpl
import Plots # hide
import NMR
import NMR.Examples

data=Examples.Data["DMEM cell culture medium"]
s=open(data["files"][1])
header,params,data = NMR.readJEOL(s)
close(s)
Plots.plot(data)
Plots.savefig("plot-JEOL-raw.svg");nothing # hide
``` 
![](plot-JEOL-raw.svg)
 
## Phase Correction
The above spectrum still shows artefacts. To clean it up, we need
to correct the phase. This can either be done manually,
by supplying `PH0`and `PH1` arguments to `FourierTransform`, or
we can resort to automatic phase correction:
```@example brukerExpl
spect = NMR.AutoPhaseCorrectChen(spect)
Plots.plot(real(spect),xaxis=:flip,xlabel="1H Chem. Shift [ppm]")
Plots.savefig("plot-specHz-ap.svg"); nothing # hide
```
![](plot-specHz-ap.svg)

## Cutting Regions
The above spectrum covers a wide range without any signals, which
is of no interest. We can chop out the central, important part:
```@example brukerExpl
s2=NMR.cut(spect,-0.5,9.0)
Plots.plot(real(s2),xaxis=:flip,xlabel="1H Chem. Shift [ppm]")
Plots.savefig("plot-S2.svg"); nothing # hide
```
![](plot-S2.svg)

## Computing Integrals
In Chemistry, it is customary to show the integral of the spectral signal
along with the spectrum. This makes the intensity of the peaks directly visible
as a step height for each signal. The integrated spectrum can be computed
like this:
```@example brukerExpl
intSpect = NMR.integrate(s2)
Plots.plot(real(s2),xaxis=:flip,xlabel="1H Chem. Shift [ppm]",label="spectrum")
Plots.plot!(20*intSpect,label="integral")
Plots.savefig("plot-S3.svg"); nothing # hide
```
![](plot-S3.svg)

## Baseline Correction
In the above plot, it can be seen that the integral line has regions where it decreases
from left to right. This results from corresponding regions in the spectrum that have negative baseline
intensity. The baseline of the NMR spectrum, while visually quite acceptable, is numerically
not perfectly adjusted. The baseline of the real part of the spectrum can be computed
as follows:
```@example brukerExpl
Plots.plot(s2,xaxis=:flip)
Plots.plot!(NMR.medianBaseline(real(s2),wdw=512),linewidth=4.0,ylims=2.0e9*[-1,10])
Plots.savefig("plot-S4.svg"); nothing # hide
```
![](plot-S4.svg)

With this baseline is subtracted from the spectrum before integration, a much cleaner
integral curve is obtained:
```@example brukerExpl
spectBc = real(s2)-NMR.medianBaseline(real(s2),wdw=512)
Plots.plot(spectBc,xaxis=:flip,ylims=2e10*[-1,10])
intSpect = NMR.integrate(spectBc,flip=true)
Plots.plot!(10*intSpect)
Plots.savefig("plot-S5.svg"); nothing # hide
```
![](plot-S5.svg)

## Detecting Peaks
The locations, heights, and widths of the peaks in the spectrum can be automatically
detected using the `NMR.peaks()` function. It returns a data structure that contains
all relevant information:
```@example brukerExpl
pks=NMR.peaks(spectBc,threshold=1e2)
```
This information can then be used for further analysis, for example for comparison
with a database, or for plotting. Here is an example:
```@example brukerExpl
Plots.plot(pks.deconvolution,xaxis=:flip,xlims=[2.9,4.1],ylims=6.0e11*[-1,10],label="deconvolution",yaxis=false,grid=false,legend=Symbol(:outer,:bottomright),minorticks=true,xlabel="1H Chem. Shift [ppm]")
Plots.plot!(spectBc+2.0e12,label="spectrum") 
Plots.plot!(pks.positions,5.5e12*ones(length(pks.positions)),seriestype=:scatter,marker=:vline,label="peaks")
Plots.savefig("plot-S6.svg"); nothing # hide
```
![](plot-S6.svg)


