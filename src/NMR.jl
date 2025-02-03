@doc raw"""
    module NMR

provides routines for processing 1D NMR spectra, as well as for simple Hilbert space
spin dynamics simulations. 

## Version History

### v0.9.0
- added JEOL data reading
- general cleanup of the code
- AutoPhaseChen() now works properly
- comlete documentation
- prepare transition to 1.0.0

### v0.8.0
- added GISSMO library interface

### v0.7.1 
- added readBrukerParameterFile()

#### v0.6.5
- added `medianBaseline()` for baseline correction
- added an entry for DSS in the HMDB table

"""
module NMR

import Pkg

   __precompile__(true);

   include("Data1D.jl")
   include("DataSet.jl")
   include("Examples.jl")
   include("AutoPhase.jl")
   include("LinearPredict.jl")
   include("PauliMatrix.jl")
   include("SpinSim.jl")
   include("Peaks.jl")
   include("HMDB.jl")
   include("GISSMO.jl")
   include("ProcessBNF.jl")
   include("Bruker.jl")
   include("Varian.jl")
   include("JEOL.jl")

  function __init__()
    println("Module NMR v0.9.0\n(c)mu 2018-2025");
  end

end
