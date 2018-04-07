module NMR

   include("Scaling.jl")
   include("SimpleGraphics.jl")
   include("SimplePlot.jl")
   include("Contour.jl")
   include("LinearPredict.jl")
   include("MapPlot.jl")
   include("PauliMatrix.jl")
   include("SpinSim.jl")
   include("Bruker.jl")
   include("DataSet.jl")
   include("AutoPhase.jl")

   export SimpleGraphics
   export Scaling
   export Contour
   export MapPlot
   export PauliMatrix
   export SpinSim
   export LinearPredict
   export SimplePlot
   export Bruker
   export AutoPhase

end
