
![Galeri Logo](doc/galeri.png)

# Galeri: Finite Element and Matrix Generation Package

_Loosely translated, the Greek term “Galeri” simply means **“Gallery”**._

The Trilinos package **Galeri** contains a suite of utilities and
classes to generate a variety of (distributed) linear systems.
Galeri’s functionalities are very close to that of the MATLAB’s
gallery() function.  Galeri works with Tpetra or Xpetra linear
algebra objects.

Several well-know finite element and finite difference matrices can
be generated using only a few simple code lines.  Consider the
Tpetra-based example below, which generates a matrix corresponding
to a <tt>nx * ny</tt> Cartesian grid, divided across <tt>mx * my</tt>
processors (so that <tt>Comm.NumProc() = mx * my</tt>).:

```
    int main(int argv, char* argc[])
    {
    #ifdef HAVE_MPI
        MPI_Init(&argv, &argc);
        Tpetra_MpiComm Comm(MPI_COMM_WORLD);
    #else
        Tpetra_SerialComm Comm;
    #endif

    Teuchos::ParameterList GaleriList;
    GaleriList.set("nx", 10 * Comm.NumProc());
    GaleriList.set("ny", 10);

    GaleriList.set("mx", Comm.NumProc());
    GaleriList.set("my", 1);

    Tpetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
    Tpetra_RowMatrix* Matrix = CreateCrsMatrix("Laplace2D", Map, GaleriList);
    ...
    #ifdef HAVE_MPI
        MPI_Finalize();
    #endif
    }
```

And that’s it! (The code snippet above only misses the header files;
see the examples in the distribution for a compilable code.) A list
of supported matrices is reported
[here](https://trilinos.github.io/docs/galeri/gl__gallery_crs_matrix.html).
The supported Tpetra_Map’s are instead
[here](https://trilinos.github.io/docs/galeri/gl__gallery_maps.html).

Galeri also contains a nice and simple finite element code, to be
used for scalar and vector elliptic-type equations using Galerkin
and SUPG discretization techniques, on both 2D and 3D unstructured
grids, composed by triangles, quads, tetrahedra and hexahedra.

Galeri’s technical documentation is maintained using Doxygen; click
[here](https://trilinos.github.io/docs/galeri/index.html) to access
the latest Doxygen documentation.

## Overview  
Understanding, validating, using and developing algorithms and
software tools for distributed linear algebra solvers, like Krylov
accelerators and preconditioners, requires input data. This data
should have three critical properties:

1.  reflect real-life applications to a sufficient extent so that we can use them to predict performance;
2.  be sufficiently simple to be amenable to mathematical analysis;
3.  it is possible to write generators that easily provide instances of these problems.

The goal of the Galeri package is exactly to produce these input
data for linear solvers.

As regards real-life applications, we have here selected PDE-type
problems. Although Galeri can generate some non-PDE linear systems,
surely there is a strong focus on PDE applications. Among them, the
Laplace problem is probably the most widely studied application,
and a broad range of results are available. For most solvers and
preconditioners, often there is a sharp convergence bound that can
be used to validate a particular class, or to get a sense of the
performances of a new method.

Galeri’s Matrix generation capabilities can help to make examples
shorter and easier to read, since they can produce challenging
linear systems in a few code lines. Therefore, the attention of the
example’s reader is not distracted by complicated instructions
aiming to build this or that linear system. Since most linear algebra
packages of Trilinos use Galeri, users will quickly become familiar
with it.

Finally, Galeri is very convenient for software testing. A large
variety of problems can be produced, to understand or validate the
performances of a given class. Galeri is used in the testing of
[MueLu](https://trilinos.github.io/docs/muelu/index.html),
[Ifpack2](https://trilinos.github.io/docs/ifpack2/index.html), and
[Zoltan2](https://trilinos.github.io/docs/zoltan2/index.html)

## Documentation

Galeri is part of the [Trilinos Project](https://trilinos.github.io),
and additional information (e.g., examples, tutorials, and source
code documentation) is available through
[Galeri's Doxygen webpages](https://trilinos.github.io/docs/galeri/index.html).

## Questions?

Contact the lead developers:

- **Galeri team**: GitHub handle: @trilinos/galeri
- **Jonathan Hu**: GitHub handle: [jhux2](https://github.com/jhux2) or jhu@sandia.gov

## Copyright and License

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For Galeri-specific copyright and license details, refer to the [galeri/COPYRIGHT](COPYRIGHT) and [galeri/LICENSE](LICENSE) files located in the `galeri` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.
