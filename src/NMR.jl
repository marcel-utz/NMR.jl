module NMR

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

   include("Bruker.jl")
   include("Varian.jl")

  function __init__()
    println(STDERR,"Module NMR v0.2\n(c)mu 2018");
  end

end
