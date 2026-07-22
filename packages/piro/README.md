# Piro: Strategy package for embedded analysis capabilitites

Piro is the top-level, unifying package of the Embedded Nonlinear Analysis Capability area. The main purpose of the package is to provide driver classes for the common uses of Trilinos nonlinear analysis tools. These drivers all can be constructed similarly, with a ModelEvaluator and a ParameterList, to make it simple to switch between different types of analysis. 
They also all inherit from the same base classes (reponse-only model evaluators) so that the resulting analysis can in turn driven by non-intrusive analysis routines.

## Documentation

Piro is part of the [Trilinos Project](https://trilinos.github.io), and additional information (e.g., examples, tutorials, and source code documentation) is available through [Piro's Doxygen webpages](https://trilinos.github.io/docs/piro/index.html).

## Questions?

Contact the lead developers:

* **Piro team**:    GitHub handle: @trilinos/piro
* **Mauro Perego**: GitHub handle: [mperego](https://github.com/mperego?) or mperego@sandia.gov

## Copyright and License

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For Piro-specific copyright and license details, refer to the [piro/COPYRIGHT](COPYRIGHT) and [piro/LICENSE](LICENSE) files located in the `piro` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.
