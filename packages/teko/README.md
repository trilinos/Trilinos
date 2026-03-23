# Teko: A package for block and physics based preconditioning

_Teko means “fuse” in Greek. This is suggestive of what Teko does, bringing together multiple physics to form one preconditioner. Additionally, the original impetus for Teko was for the development of preconditioners for simulation of magnetohydrodynamics and fusion reactors._

Teko is a package for development and implementation of block preconditioners. This includes support for manipulation and setup of block operators. Furthermore tools exist to support decomposition of a fully coupled operator. Additionally, facilities that allow the construction of approximate inverse operators using the full complement of available preconditioners and solvers are available in Teko. Finally, a small number of generic block preconditioners has been implemented in Teko, including block Jacobi, and block Gauss-Seidel. For the Navier-Stokes equation, Teko has implementations of SIMPLE, PCD and LSC. 
For details on these methods see [Stabilization and Scalable Block Preconditioning for the Navier-Stokes Equations](http://dx.doi.org/10.1016/j.jcp.2011.09.001) and the references therein.

## Documentation

Teko is part of the [Trilinos Project](https://trilinos.github.io), and additional information (e.g., examples, tutorials, and source code documentation) is available through [Teko's Doxygen webpages](https://trilinos.github.io/docs/teko/index.html).

## Questions?

Contact the lead developers:

- **Teko team**:        GitHub handle: @trilinos/teko
- **Malachi Phillips**: GitHub handle: [MalachiTimothyPhillips](https://github.com/MalachiTimothyPhillips) or malphil@sandia.gov
- **Christian Glusa**:  GitHub handle: [cgcgcg](https://github.com/cgcgcg) or caglusa@sandia.gov
 
## Copyright and License

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For Teko-specific copyright and license details, refer to the [teko/COPYRIGHT](COPYRIGHT) and [teko/LICENSE](LICENSE) files located in the `teko` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.
