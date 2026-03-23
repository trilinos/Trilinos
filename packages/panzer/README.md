# Panzer 

A partial differential equation assembly engine for strongly
coupled complex multiphysics systems.

Panzer provides global tools for finite element analysis. It handles continuous and discontinuous high-order compatible finite elements, as implemented in [Intrepid2](https://github.com/trilinos/Trilinos/blob/master/packages/intrepid2/README.md) on unstructured meshes. Panzer relies on [Phalanx](https://github.com/trilinos/Trilinos/blob/master/packages/phalanx/README.md) to manage with efficiency and flexibility the assembly of complex problems. Panzer also enables the solution of nonlinear problems, by interfacing with several Trlinos linear and nonlinear solvers. It computes derivatives and sensitivities through automatic differentiation ([Sacado](https://github.com/trilinos/Trilinos/blob/master/packages/sacado/README.md)). It supports [Tpetra](https://github.com/trilinos/Trilinos/blob/master/packages/tpetra/README.md) data structures and achieves performance portability through the [Kokkos](https://kokkos.org/) programming model.

For an introductory tutorial to what Panzer is, how it works, and how you use it, check out this 
[presentation](https://trilinos.github.io/pdfs/siamCseTalk.pdf). 
You can examine the associated code by [cloning Trilinos](https://github.com/trilinos/trilinos "git clone git@github.com:trilinos/Trilinos") 
and then poking around in `Trilinos/packages/panzer/adapters-stk/tutorial/siamCse17/`.

## Subpackage Layout

core: Basic utilities for all subpackages
dof-mgr: degree of freedom manager
disc-fe: finite element discretization tools
adapters-stk: stk adapters


## Documentation

Panzer is part of the [Trilinos Project](https://trilinos.github.io), and additional information (e.g., examples, tutorials, and source code documentation) is available through [Panzer's Doxygen webpages](https://trilinos.github.io/docs/panzer/index.html).

## Questions?

Contact the lead developers:

- **Panzer Team**:     GitHub handle @trilinos/panzer
- **Roger Pawlowski**: GitHub handle [rppawlo](https://github.com/rppawlo), email: rppawlo@sandia.gov
- **Eric C. Cyr**:     GitHub handle: [eric-c-cyr](https://github.com/eric-c-cyr) or eccyr@sandia.gov

## Copyright and License

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For Panzer-specific copyright and license details, refer to the [panzer/COPYRIGHT](COPYRIGHT) and [panzer/LICENSE](LICENSE) files located in the `panzer` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.
