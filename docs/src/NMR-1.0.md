# Roadmap to Version 1.0

This is a roadmap towards publishing `NMR.jl` first for internal use, and then 
to the wider community.

In its current state, `NMR.jl`is a package of tools that have been designed for
only internal use. Most of it is not properly documented. The data structures
have not been thought through carefully, but have mostly just grown as things went
along.

For professional use, a few requirements are imperative:
- good documentation, including a set of example data
- clean integration into the Julia universe
    - Data types
    - Plots
    - Input and output; IJulia etc
- well thought-out `NMRData` type
- seamless integration of spectral analysis and simulation tools
- loading multiple vendor data formats
- proper handling of metadata

The current roadmap for development is shown in the following table:

| Version | Planned features  | Target release date |
|----------|-------------------|----------------------|
| 0.9.0   | cleanup, document basic functionality, how-to documentation for simple spectral processing | 22.1.2025|
| 0.10.0  | new data structures, transform, plots   | 31.1.2025
| 0.11.0  | new vendor and metadata support | 28.2.2025 |
| 0.12.0  | integration of AI processing    | 31.3.2025 |
| 0.13.0  | complete documentation          | 15.4.2025 |

