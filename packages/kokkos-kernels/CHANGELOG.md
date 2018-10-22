# Change Log

## [2.7.00](https://github.com/kokkos/kokkos-kernels/tree/2.7.00) (2018-05-24)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/2.6.00...2.7.00)

**Implemented enhancements:**

- Tests: add capability to build a unit test standalone [\#233](https://github.com/kokkos/kokkos-kernels/issues/233)
- Make KokkosKernels work without KOKKOS\_ENABLE\_DEPRECATED\_CODE [\#223](https://github.com/kokkos/kokkos-kernels/issues/223)
- Replace KOKKOS\_HAVE\_\* FLAGS with KOKKOS\_ENABLE\_\* [\#219](https://github.com/kokkos/kokkos-kernels/issues/219)
- Add team-based scal, mult, update, nrm2 [\#214](https://github.com/kokkos/kokkos-kernels/issues/214)
- Add team based abs [\#209](https://github.com/kokkos/kokkos-kernels/issues/209)
- Generated CPP files moving includes inside the ifdef's [\#199](https://github.com/kokkos/kokkos-kernels/issues/199)
- Implement BlockCRS in Kokkoskernels [\#184](https://github.com/kokkos/kokkos-kernels/issues/184)
- Spgemm hash promotion [\#171](https://github.com/kokkos/kokkos-kernels/issues/171)
- Batched BLAS enhancement [\#170](https://github.com/kokkos/kokkos-kernels/issues/170)
- Document & check CMAKE\_CXX\_USE\_RESPONSE\_FILE\_FOR\_OBJECTS=ON in CUDA build [\#148](https://github.com/kokkos/kokkos-kernels/issues/148)

**Fixed bugs:**

- Update drivers in perf\_tests/graph to use Kokkos::initialize\(\) [\#200](https://github.com/kokkos/kokkos-kernels/issues/200)
- unit tests failing/hanging on Volta [\#188](https://github.com/kokkos/kokkos-kernels/issues/188)
- Inner TRSM: SIMD build error; manifests in Ifpack2 [\#183](https://github.com/kokkos/kokkos-kernels/issues/183)
- d2\_graph\_color doesn't have a default coloring mechanism [\#168](https://github.com/kokkos/kokkos-kernels/issues/168)
- Unit tests do not build with Serial backend [\#154](https://github.com/kokkos/kokkos-kernels/issues/154)


## [2.6.00](https://github.com/kokkos/kokkos-kernels/tree/2.6.00) (2018-03-07)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/2.5.00...2.6.00)

**Implemented enhancements:**

- Spgemm hash promotion [\#171](https://github.com/kokkos/kokkos-kernels/issues/171)
- Batched BLAS enhancement [\#170](https://github.com/kokkos/kokkos-kernels/issues/170)

**Fixed bugs:**

- d2\_graph\_color doesn't have a default coloring mechanism [\#168](https://github.com/kokkos/kokkos-kernels/issues/168)
- Build error when MKL TPL is enabled [\#135](https://github.com/kokkos/kokkos-kernels/issues/135)


## [2.5.00](https://github.com/kokkos/kokkos-kernels/tree/2.5.00) (2017-12-15)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/0.10.03...2.5.00)

**Implemented enhancements:**

- KokkosBlas:   Add GEMM interface  [\#105](https://github.com/kokkos/kokkos-kernels/issues/105)
- KokkosBlas:   Add GEMM default Kernel [\#125](https://github.com/kokkos/kokkos-kernels/issues/125)
- KokkosBlas:   Add GEMV that wraps BLAS \(and cuBLAS\) [\#16](https://github.com/kokkos/kokkos-kernels/issues/16)
- KokkosSparse: Make SPMV test not print GBs of output if something goes wrong.  [\#111](https://github.com/kokkos/kokkos-kernels/issues/111)
- KokkosSparse: ETI SpGEMM and Gauss Seidel and take it out of Experimental namespace [\#74](https://github.com/kokkos/kokkos-kernels/issues/74)
- BuildSystem:  Fix Makesystem to correctly build library after aborted install [\#104](https://github.com/kokkos/kokkos-kernels/issues/104)
- BuildSystem:  Add option ot generate\_makefile.bash to define memoryspaces for instantiation [\#89](https://github.com/kokkos/kokkos-kernels/issues/89)
- BuildSystem:  generate makefile tpl option [\#66](https://github.com/kokkos/kokkos-kernels/issues/66)
- BuildSystem:  Add a simpler compilation script, README update etc [\#96](https://github.com/kokkos/kokkos-kernels/issues/96)

**Fixed bugs:**

- Internal Compiler Error GCC in GEMM [\#129](https://github.com/kokkos/kokkos-kernels/issues/129)
- Batched Team LU: bug for small team\_size [\#110](https://github.com/kokkos/kokkos-kernels/issues/110)
- Compiler BUG in IBM XL pragma unrolling [\#92](https://github.com/kokkos/kokkos-kernels/issues/92)
- Fix Blas TPL enables build [\#77](https://github.com/kokkos/kokkos-kernels/issues/77)
- Batched Gemm Failure [\#73](https://github.com/kokkos/kokkos-kernels/issues/73)
- CUDA 7.5 \(GCC 4.8.4\) build errors [\#72](https://github.com/kokkos/kokkos-kernels/issues/72)
- Cuda BLAS tests fail with UVM if CUDA\_LAUNCH\_BLOCKING=1 is not defined on Kepler [\#51](https://github.com/kokkos/kokkos-kernels/issues/51)
- CrsMatrix: sumIntoValues and replaceValues incorrectly count the number of valid column indices. [\#11](https://github.com/kokkos/kokkos-kernels/issues/11)
- findRelOffset test assumes UVM [\#32](https://github.com/kokkos/kokkos-kernels/issues/32)

## [0.10.03](https://github.com/kokkos/kokkos-kernels/tree/0.10.03) (2017-09-11)

**Implemented enhancements:**

- KokkosSparse: Fix unused variable warnings in spmv\_impl\_omp, spmv Test and graph color perf\_test [\#63](https://github.com/kokkos/kokkos-kernels/issues/63)
- KokkosBlas:  dot: Add unit test [\#15](https://github.com/kokkos/kokkos-kernels/issues/15)
- KokkosBlas:  dot: Add special case for multivector \* vector \(or vector \* multivector\) [\#13](https://github.com/kokkos/kokkos-kernels/issues/13)
- BuildSystem: Make KokkosKernels build independently of Trilinos [\#1](https://github.com/kokkos/kokkos-kernels/issues/1)
- BuildSystem: Fix ETI System not to depend on Tpetra ETI [\#5](https://github.com/kokkos/kokkos-kernels/issues/5)
- BuildSystem: Change CMake to work with new ETI system [\#19](https://github.com/kokkos/kokkos-kernels/issues/19)
- BuildSystem: Fix TpetraKernels names to KokkosKernels [\#4](https://github.com/kokkos/kokkos-kernels/issues/4)
- BuildSystem: Trilinos/KokkosKernels reports no ETI in almost any circumstance [\#29](https://github.com/kokkos/kokkos-kernels/issues/29)
- General:     Kokkos::ArithTraits\<double\>::nan\(\) is very slow [\#35](https://github.com/kokkos/kokkos-kernels/issues/35)
- General:     Design and Define New UnitTest infrastructure [\#28](https://github.com/kokkos/kokkos-kernels/issues/28)
- General:     Move Tpetra::Details::OrdinalTraits to KokkosKernels [\#22](https://github.com/kokkos/kokkos-kernels/issues/22)
- General:     Rename files and NameSpace to KokkosKernels [\#12](https://github.com/kokkos/kokkos-kernels/issues/12)
- General:      PrepareStandalone: Get rid of Teuchos usage [\#2](https://github.com/kokkos/kokkos-kernels/issues/2)
- General:      Fix warning with char being either signed or unsigned in ArithTraits [\#60](https://github.com/kokkos/kokkos-kernels/issues/60)
- Testing:      Make all tests run with -Werror [\#68](https://github.com/kokkos/kokkos-kernels/issues/68)

**Fixed bugs:**

- SPGEMM Test Fails for Cuda when compiled through Trilinos  [\#49](https://github.com/kokkos/kokkos-kernels/issues/49)
- Fix ArithTraits min for floating points [\#47](https://github.com/kokkos/kokkos-kernels/issues/47)
- Pthread ETI error [\#25](https://github.com/kokkos/kokkos-kernels/issues/25)
- Fix CMake Based ETI for Threads backend [\#46](https://github.com/kokkos/kokkos-kernels/issues/46)
- KokkosKernels\_ENABLE\_EXPERIMENTAL causes build error  [\#59](https://github.com/kokkos/kokkos-kernels/issues/59)
- ArithTraits warnings in CUDA build [\#71](https://github.com/kokkos/kokkos-kernels/issues/71)
- Graph coloring build warnings [\#3](https://github.com/kokkos/kokkos-kernels/issues/3)




\* *This Change Log was automatically generated by [github_changelog_generator](https://github.com/skywinder/Github-Changelog-Generator)*
