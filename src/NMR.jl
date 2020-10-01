module NMR

import Pkg

   __precompile__(false);

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
    println("Module NMR\n(c)mu 2018,2019,2020");
  end

end
