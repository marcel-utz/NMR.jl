module NMR

   include("Scaling.jl")
   include("SimpleGraphics.jl")
   include("SimplePlot.jl")
   include("Contour.jl")
   include("LinearPredict.jl")
   include("MapPlot.jl")
   include("PauliMatrix.jl")
   include("SpinSim.jl")

   export SimpleGraphics
   export Scaling
   export Contour
   export MapPlot
   export PauliMatrix
   export SpinSim
   export LinearPredict
   export SimplePlot

end
