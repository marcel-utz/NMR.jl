@doc raw"""
    module NMR

provides routines for processing 1D NMR spectra, as well as for simple Hilbert space
spin dynamics simulations. 

### Version History

#### v0.6.5
- added `medianBaseline()` for baseline correction
- added an entry for DSS in the HMDB table

"""
module NMR

import Pkg

   __precompile__(true);

   include("DataSet.jl")
   include("AutoPhase.jl")
   include("Scaling.jl")
   include("SimpleGraphics.jl")
   include("SimplePlot.jl")
   include("Contour.jl")
   include("LinearPredict.jl")
   include("MapPlot.jl")
   include("PauliMatrix.jl")
   include("SpinSim.jl")
   include("Peaks.jl")
   include("HMDB.jl")

   include("Bruker.jl")
   include("Varian.jl")

  function __init__()
    println("Module NMR v0.6.5\n(c)mu 2018-2021");
  end

end
