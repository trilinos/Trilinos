Requirements
############

Kokkos Kernels only requirement is a compatible installation of the Kokkos Core programming model.
In general the current version of Kokkos Kernels will work with the current version of Kokkos Core and the previous minor version of Kokkos Core.
The compiler/toolchain requirements are the same as Kokkos Core, in fact you should build both libraries with the same compiler and toolchain.
The current Kokkos Core requirements can be found `here <https://kokkos.org/kokkos-core-wiki/requirements.html>`_.

Third Party Libraries (TPLs)
============================

Kokkos Kernels provides access to numerous TPLs that can be declared at configuration time.
For vendor libraries we recommend to use the version that is compatible with the compiler/toolchain version you are using.
For BLAS libraries simply use a library that supports the BLAS specification.
