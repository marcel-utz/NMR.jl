module NMR

import Pkg

#   __precompile__(true);

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

   include("NeuralNet.jl")

  function __init__()
    println("Module NMR $(Pkg.installed()["NMR"])\n(c)mu 2018,2019");
  end

end