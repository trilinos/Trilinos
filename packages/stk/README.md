# STK: Sandia Toolkit

The STK modules provide infrastructure to support the development of
computational engineering applications.

STK is composed of several modules (sub-packages):
* STK Util: utilities such as parallel communication, command line parsing, etc.
* STK Topology: definitions of mesh-entity (elements, sides, etc) node orderings, side orderings.
* STK Mesh: parallel unstructured mesh
* STK IO: reading/writing of STK Mesh to/from Exodus files
* STK Coupling: support for MPMD coupling of MPI applications
* STK Search: geometric proximity bounding-box search
* STK Transfer: copy solution field values between meshes
* STK Balance: parallel load partitioning and dynamic rebalancing
* STK Middle Mesh: common refinement of two surface meshes
* STK SIMD: a general interface to vector instructions such as SSE, AVX, etc.
* STK ExprEval: string function expression evaluation

Some STK modules depend on others, but in general it is not necessary to use all
of them together.
For example, applications can use STK Search and STK Transfer without using STK
Mesh. Also, STK Util does not depend on any other STK modules.

## Documentation

[STK Manual (pdf)](https://trilinos.github.io/pdfs/STKManual_2024-07-12-final.pdf)

## Questions?

Contact the developers:

- **STK team**:    STK-NGPTeam@sandia.gov

## Copyright and License

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For Trilinos developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.
