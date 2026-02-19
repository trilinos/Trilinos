# MiniTensor Package

MiniTensor is a library for the use, manipulation, algebra and optimization of small vectors and tensors and problems that depend on them. It was first conceived and developed for the implementation of complex material models in Albany. Its purpose is to provide a compact representation of vector and tensor expressions. Its emphasis is on ease of use, and accurate algorithms, specifically those used for the development of constitutive models in finite deformation solid mechanics. It is smaller, but more focused than Blitz++, TVMet or Eigen. It is used in both research (Albany/LCM) and production (Sierra/SM) environments.

**Features**:

*    Vectors, second-order to fourth-order tensors.
*    Static (for production) and dynamic (for research) storage.
*    Basic manipulation, linear algebra, geometry and mechanics.
*    Unconstrained (native, ROL) and constrained (ROL) optimization.
*    Solution of nonlinear systems of equations.
*    Emphasis in accurate algorithms.
*    Fully templated, plays well with Sacado.

## Documentation

MiniTensor is part of the [Trilinos Project](https://trilinos.github.io), and additional information (e.g., examples, tutorials, and source code documentation) is available through [MiniTensor's Doxygen webpages](https://trilinos.github.io/docs/minitensor/index.html).

## Questions?

Contact the lead developers:

- **MiniTensor Team**: GitHub handle @trilinos/minitensor
- **Alejandro Mota**:  (GitHub handle: [lxmota](https://github.com/lxmota) or amota@sandia.gov)

## Copyright and License

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For MiniTensor-specific copyright and license details, refer to the [minitensor/COPYRIGHT](COPYRIGHT) and [minitensor/LICENSE](LICENSE) files located in the `minitensor` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.
