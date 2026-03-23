# Thyra: Interfaces and Support for Abstract Numerical Algorithms

The Thyra package contains a set of interfaces and supporting code
that defines basic interoperability mechanisms between different
types of abstract numerical algorithm (ANA) software. The foundations
for all Thyra interfaces are the mathematical concepts of vectors,
vector spaces, and linear operators. All other ANA interfaces and
support software are built on these fundamental operator/vector
interfaces.

## Documentation

Thyra is part of the [Trilinos Project](https://trilinos.github.io),
and additional information (e.g., examples, tutorials, and source
code documentation) is available through 
[Thyra's Doxygen webpages](https://trilinos.github.io/docs/thyra/index.html).

The following document describes the basic ideas behind Thyra and
provides an overview of the operator/vector interfaces:

 *   Bartlett, Roscoe. _Thyra Linear Operators and Vectors: Overview
     of Interfaces and Support Software for the Development and
     Interoperability of Abstract Numerical Algorithms._ SAND2007-5984,
     Sandia National Laboratories, 2007
     [[PDF](https://bartlettroscoe.github.io/publications/ThyraOverview2007.pdf)]

The primary Thyra ANA interfaces are broadly layered as followed:

 *   **Operator/vector interfaces**
     *   [Thyra::VectorBase](https://trilinos.github.io/docs/thyra/class_thyra_1_1_vector_base.html)
     *   [Thyra::VectorSpaceBase](https://trilinos.github.io/docs/thyra/class_thyra_1_1_vector_space_base.html)
     *   [Thyra::LinearOpBase](https://trilinos.github.io/docs/thyra/class_thyra_1_1_linear_op_base.html)
     *   [Thyra::MultiVectorBase](https://trilinos.github.io/docs/thyra/class_thyra_1_1_multi_vector_base.html)
 *   **Operator solve interfaces**
     *   [Thyra::PreconditionerBase](https://trilinos.github.io/docs/thyra/class_thyra_1_1_preconditioner_base.html)
     *   [Thyra::PreconditionerFactoryBase](https://trilinos.github.io/docs/thyra/class_thyra_1_1_preconditioner_factory_base.html)
     *   [Thyra::LinearOpWithSolveBase](https://trilinos.github.io/docs/thyra/class_thyra_1_1_linear_op_with_solve_base.html)
     *   [Thyra::LinearOpWithSolveFactoryBase](https://trilinos.github.io/docs/thyra/class_thyra_1_1_linear_op_with_solve_factory_base.html)
 *   **Nonlinear Interfaces**
     *   [Thyra::ModelEvaluator](https://trilinos.github.io/docs/thyra/class_thyra_1_1_model_evaluator.html)
     *   [Thyra::NonlinearSolverBase](https://trilinos.github.io/docs/thyra/class_thyra_1_1_nonlinear_solver_base.html)

A few important points about Thyra interfaces are:

 *   All interfaces are expressed as abstract C++ base classes (i.e. object-oriented)
 *   All interfaces are templated on a Scalar data type (i.e. generic)
 *   All memory management is performed using [Teuchos](https://trilinos.github.io/docs/teuchos/namespace_teuchos.html) memory management classes involving no raw C++ pointers (see below)

For each of these sets of interfaces, the Thyra package also a set of general adapter and general support software. See the Thyra package documentation for more details.

There are several Trilinos packages that implement ANAs, and/or can accept input as ANA object, and/or can provide implementations of ANA objects for other ANAs to use. Information related to Thyra and ANAs can be found below:

 *   **[Thyra](thyra.html)**
     *   Adapters
         *   [Thyra Adapters](https://trilinos.github.io/docs/thyra/group___tpetra___thyra___op___vec__adapters__grp.html):
             *   Utility classes and functions for mapping between Thyra wrapper objects, and Tpetra linear algebra objects

Thyra Coding and Documentation Guidelines

*   Bartlett, Roscoe. _Thyra Coding and Documentation Guidelines (TCDG)._ Sandia National Laboratories [[PDF](https://bartlettroscoe.github.io/publications/ThyraCodingGuideLines.pdf)]

## Questions?

Contact the lead developers:

- **Thyra Team**:         GitHub handle @trilinos/thyra
- **Roscoe A. Bartlett**: GitHub handle: [bartlettroscoe](https://github.com/bartlettroscoe) or rabartl@sandia.gov

## Copyright and License

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For Thyra-specific copyright and license details, refer to the [thyra/COPYRIGHT](COPYRIGHT) and [thyra/LICENSE](LICENSE) files located in the `thyra` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.
