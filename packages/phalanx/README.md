# Phalanx

A Partial Differential Equation Field Evaluation 
Kernel for Flexible Management of Complex Dependency Chains.

Phalanx is a local field evaluation kernel specifically designed for general partial differential equation solvers. 
The main goal of Phalanx is to decompose a complex problem into a number of simpler problems with managed dependencies to support rapid development and extensibility of the PDE code. 
Through the use of template metaprogramming concepts, Phalanx supports arbitrary user defined data types and evaluation types. This allows for unprecedented flexibility for direct integration with user applications and 
provides extensive support for embedded technology such as automatic differentiation for sensitivity analysis and uncertainty quantification.

## Documentation

Phalanx is part of the [Trilinos Project](https://trilinos.github.io), and additional information (e.g., examples, tutorials, and source code documentation) is available through [Phalanx's Doxygen webpages](https://trilinos.github.io/docs/phalanx/index.html).

### Locally Building Doxygen Pages

To generate the Doxygen HTML Documentation:

1. Go to the doc directory:
> cd phalanx/doc
2. Build the documentation:
> doxygen Doxyfile
3. Open the main web page in your browser:
phalanx/doc/html/index.html

## Questions?

Contact the lead developers:

- **Phalanx team**: GitHub handle: @trilinos/phalanx
- **Roger Pawlowski**: GitHub handle [rppawlo](https://github.com/rppawlo), email: rppawlo@sandia.gov

## Copyright and License

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For Phalanx-specific copyright and license details, refer to the [phalanx/COPYRIGHT](COPYRIGHT) and [phalanx/LICENSE](LICENSE) files located in the `phalanx` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.
