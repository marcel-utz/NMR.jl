*A set of tools for processing, plotting, and interpreting NMR data*

# Overview

## Package Features

- Reading and processing NMR data from JEOL, Bruker, and Varian/Agilent Spectrometers
- Simulations of Spin dynamics based on efficient sparse matrix representations of the spin Hamiltonian and density operator
- Some specialised and advanced tools for data interpretation in the context of metabolomics

## Contents
```@contents
```

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


## Manual Outline
If you would like to jump right in, read the Getting started section. For a more in-depth and complete
documentation, refer to the Manual. Finally, a complete list and documentation of every function 
can be found under API.

## Feedback
`NMR.jl`is still under active development, and we would appreciate your feedback, including feature requests,
bug reports, and general comments. Please contact `marcel.utz@kit.edu` by email.

## Citing `NMR.jl`
If you publish your work and have been using NMR.jl to process, present, and/or interpret your data, we would
appreciate if you could acknowledge this by citing our work. The reference is TBA.

