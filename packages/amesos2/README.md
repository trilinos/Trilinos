# Amesos2: Templated Direct Sparse Solver Package


Amesos2 is a package for solving sparse linear systems using direct solvers. It differs from Amesos in that it is templated on the scalar and index types.

Amesos2 provides two native solvers:

*   KLU2 (Trilinos 11.12 and higher)
*   Basker (Trilinos 11.14 and higher)

KLU2 is enabled by default (12.4 and higher releases). Basker needs to be enabled during configure time if needed. All other solvers must be installed as third-party libraries (TPLs). Currently, Amesos2 supports these TPLs:

*   SuperLU
*   SuperLU-MT
*   SuperLU-Dist (exclude 9.0 & 9.1 - known segfault issues)
*   MUMPS
*   Pardiso
*   LAPACK (for testing purposes)
*   CSS (MKL)
*   Cholmod
*   Umfpack
*   cuSOLVER
*   Tacho
*   ShyLU-Basker
*   KLU2 (Trilinos 11.12 and higher)

Note that Amesos2 can only support scalar types that are supported by the desired third-party package (typically, float, double, complex, and double complex).

## Publications

If you use Amesos2 in your applications, please cite Amesos2 using the following publication:

*   [Amesos2 and Belos: Direct and iterative solvers for large sparse linear systems](http://dx.doi.org/10.3233/SPR-2012-0352 "Download 11.10"), Eric Bavier, Mark Hoemmen, Sivasankaran Rajamanickam, Heidi Thornquist, Scientific Programming, 2012, Volume 20, Issue 3.

Amesos2 is part of the 2nd generation, templated Trilinos stack, so it is templated both on the scalar type and the index type. We are adding support for more solvers. Please contact the developers if you have a specific request.

## Questions? 

Contact lead developers:

* **Amesos2 team**           (GitHub handle: @trilinos/amesos2)
* **Ichitaro Yamazaki**      (GitHub handle: [@iyamazaki](https://github.com/iyamazaki) or iyamaza@sandia.gov)

Other contributors:

* **Nathan Ellingwood**      (GitHub handle: [ndellingwood](https://github.com/ndellingwood)
* **Siva Rajamanickam**      (GitHub handle: [srajama1](https://github.com/srajama1) or srajama@sandia.gov)
* **Lisa Claus** (STRUMPACK) (lclaus@lbl.gov)

## Copyright and License

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For Amesos2-specific copyright and license details, refer to the [amesos2/COPYRIGHT](COPYRIGHT) and [amesos2/LICENSE](LICENSE) files located in the `amesos2` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.
