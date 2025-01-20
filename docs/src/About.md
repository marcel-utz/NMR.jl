# About NMR.jl

`NMR.jl` is a library of tools for the processing, plotting, and interpretation
of NMR data. The project started out as an internal quick-and-dirty set of tools
in the Utz research laboratory at the University of Southampton. At the time,
Julia was in a very early stage, and its use was experimental. In the meantime,
a sizeable number of research projects have been carried out in the group using
early versions of NMR.jl, and it seemed like a good idea to complete the tools
and the documentation, and make the package available to the community.

!!! note "Applications"
    `NMR.jl` is intended to be applicable to any kind of NMR data, and aims at
    implementing a broad set of features to enable even advanced NMR data processing.
    Some of its functionality is designed for convenience, making everyday tasks 
    accessible with a reasonable default set of parameters.
    Power users may want to use the more low-level routines, which are underlying
    the high-level interface.
