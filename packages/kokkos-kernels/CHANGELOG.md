# Change Log

## [3.2.00](https://github.com/kokkos/kokkos-kernels/tree/3.2.00) (2020-08-19)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/3.1.01...3.2.00)

**Implemented enhancements:**

- Add CudaUVMSpace specializations for cuBLAS IAMAX and SCAL [\#758](https://github.com/kokkos/kokkos-kernels/issues/758)
- Add wiki examples [\#735](https://github.com/kokkos/kokkos-kernels/issues/735)
- Support complex\_float, complex\_double in cuSPARSE SPMV wrapper [\#726](https://github.com/kokkos/kokkos-kernels/issues/726)
- Add performance tests for trmm and trtri [\#711](https://github.com/kokkos/kokkos-kernels/issues/711)
- SpAdd requires output values to be zero-initialized, but this shouldnt be needed [\#694](https://github.com/kokkos/kokkos-kernels/issues/694)
- SpAdd doesnt merge entries correctly [\#685](https://github.com/kokkos/kokkos-kernels/issues/685)
- cusparse SpMV merge algorithm [\#670](https://github.com/kokkos/kokkos-kernels/issues/670)
- TPL support for SpMV [\#614](https://github.com/kokkos/kokkos-kernels/issues/614)
- Add two BLAS/LAPACK calls needed by: Sptrsv supernode \#552 [\#589](https://github.com/kokkos/kokkos-kernels/issues/589)
- HashmapAccumulator has several unused members, misnamed parameters [\#508](https://github.com/kokkos/kokkos-kernels/issues/508)

**Fixed bugs:**

- Nightly test failure: spgemm unit tests failing on White \(Power8\) [\#780](https://github.com/kokkos/kokkos-kernels/issues/780)
- supernodal does not build with UVM enabled [\#633](https://github.com/kokkos/kokkos-kernels/issues/633)

## [3.1.1](https://github.com/kokkos/kokkos-kernels/tree/3.1.1) (2020-05-04)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/3.1.00...3.1.1)

** Fixed bugs:** 

- KokkosBatched QR PR breaking nightly tests [\#691](https://github.com/kokkos/kokkos-kernels/issues/691)

## [3.1.00](https://github.com/kokkos/kokkos-kernels/tree/3.1.00) (2020-04-14)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/3.0.00...3.1.00)

**Implemented enhancements:**

- Two-stage & Classical Gauss-Seidel [\#672](https://github.com/kokkos/kokkos-kernels/issues/672)
- Test transpose utilities [\#664](https://github.com/kokkos/kokkos-kernels/issues/664)
- cuSPARSE spmv wrapper doesn't actually use 'mode' [\#650](https://github.com/kokkos/kokkos-kernels/issues/650)
- Distance-2 improvements [\#625](https://github.com/kokkos/kokkos-kernels/issues/625)
- FindMKL module: which mkl versions to prioritize [\#480](https://github.com/kokkos/kokkos-kernels/issues/480)
- Add SuperLU as optional CMake TPL [\#545](https://github.com/kokkos/kokkos-kernels/issues/545)
- Revamp the ETI system [\#460](https://github.com/kokkos/kokkos-kernels/issues/460)

**Fixed bugs:**

- 2-stage GS update breaking cuda/10+rdc build [\#673](https://github.com/kokkos/kokkos-kernels/issues/673)
- Why CrsMatrix::staticcrsgraph\_type uses execution\_space and not device\_type? [\#665](https://github.com/kokkos/kokkos-kernels/issues/665)
- TRMM and TRTRI build failures with clang/7+cuda9+Cuda\_OpenMP and gcc/5.3+OpenMP [\#657](https://github.com/kokkos/kokkos-kernels/issues/657)
- cuSPARSE spmv wrapper doesn't actually use 'mode' [\#650](https://github.com/kokkos/kokkos-kernels/issues/650)
- Block Gauss-Seidel test fails when cuSPARSE is enabled [\#648](https://github.com/kokkos/kokkos-kernels/issues/648)
- cuda uvm test failures without launch blocking - expected behavior? [\#636](https://github.com/kokkos/kokkos-kernels/issues/636)
- graph\_color\_d2\_symmetric\_double\_int\_int\_TestExecSpace seg faults in cuda/10.1 + Volta nightly test on kokkos-dev-2 [\#634](https://github.com/kokkos/kokkos-kernels/issues/634)
- Build failures on kokkos-dev with clang/7.0.1 cuda/9.2 and blas/cublas/cusparse tpls [\#629](https://github.com/kokkos/kokkos-kernels/issues/629)
- Distance-2 improvements [\#625](https://github.com/kokkos/kokkos-kernels/issues/625)
- trsv - internal compiler error with intel/19 [\#607](https://github.com/kokkos/kokkos-kernels/issues/607)
- complex\_double misalignment still breaking SPGEMM [\#598](https://github.com/kokkos/kokkos-kernels/issues/598)
- PortableNumericCHASH can't align shared memory  [\#587](https://github.com/kokkos/kokkos-kernels/issues/587)
- Remove all references to Kokkos::Impl::is\_same [\#586](https://github.com/kokkos/kokkos-kernels/issues/586)
- Can I run KokkosKernels spgemm with float or int32 type? [\#583](https://github.com/kokkos/kokkos-kernels/issues/583)
- Kokkos Blas: gemv segfaults [\#443](https://github.com/kokkos/kokkos-kernels/issues/443)
- Generated kokkos-kernels file names are too long and are crashing cloning Trilinos on Windows [\#395](https://github.com/kokkos/kokkos-kernels/issues/395)


## [3.0.00](https://github.com/kokkos/kokkos-kernels/tree/3.0.00) (2020-01-27)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/2.9.00...3.0.00)

**Implemented enhancements:**

- BuildSystem: Standalone Modern CMake support [\#491](https://github.com/kokkos/kokkos-kernels/pull/491)
- Cluster GS and SGS: add cluster gauss-seidel implementation [\#455](https://github.com/kokkos/kokkos-kernels/pull/455)
- spiluk: Add sparse ILUK implementation [\#459](https://github.com/kokkos/kokkos-kernels/pull/459)
- BLAS gemm: Dot-based GEMM Cuda optimization for C = betaC + alphaA^TB - [\#490]https://github.com/kokkos/kokkos-kernels/pull/490)
- Sorting utilities: [\#461](https://github.com/kokkos/kokkos-kernels/pull/461)
- SGS: Support multiple rhs in SGS efficiently [\#488](https://github.com/kokkos/kokkos-kernels/issues/488)
- BLAS trsm: Add support and interface for trsm [\#513](https://github.com/kokkos/kokkos-kernels/issues/513)
- BLAS iamax: Implement iamax [\#87](https://github.com/kokkos/kokkos-kernels/issues/87)
- BLAS gesv: [\#449](https://github.com/kokkos/kokkos-kernels/issues/449)
- sptrsv supernodal: Add supernodal sparse triangular solver [\#552](https://github.com/kokkos/kokkos-kernels/pull/552)
- sptrsv: Add cusparse tpl support for sparse triangular solve, cudagraphs to fallback [\#555](https://github.com/kokkos/kokkos-kernels/pull/555)
- KokkosGraph: Output colors assigned during graph coloring [\#444](https://github.com/kokkos/kokkos-kernels/issues/444)
- MatrixReader: Full matrix market support [\#466](https://github.com/kokkos/kokkos-kernels/pull/466)

**Fixed bugs:**

- gemm: Fix bug for complex types in fallback impl [\#550](https://github.com/kokkos/kokkos-kernels/pull/550)
- gemv: Fix degenerate matrix cases [\#514](https://github.com/kokkos/kokkos-kernels/pull/514)
- spgemm: Fix cuda build with complex\_double misaligned shared memory access [\#500](https://github.com/kokkos/kokkos-kernels/issues/500)
- spgemm: Wrong team size heuristic used for SPGEMM when Kokkos deprecated=OFF [\#474](https://github.com/kokkos/kokkos-kernels/issues/474)
- dot: Improve accuracy for float and complex_float [\#574](https://github.com/kokkos/kokkos-kernels/issues/574)
- SpMV Struct: Fix bug with intel\_17\_0\_1 [\#456](https://github.com/kokkos/kokkos-kernels/issues/456)
- readmtx: Fix invalid read due to loop condition [\#453](https://github.com/kokkos/kokkos-kernels/issues/453)
- spgemm: Fix hashmap accumulator bug yielding crashes and wrong results [\#402](https://github.com/kokkos/kokkos-kernels/issues/402)
- KokkosGraph: Fix distance-1 graph coloring segfault [\#275](https://github.com/kokkos/kokkos-kernels/issues/275)
- UniformMemoryPool: does not re-initialize chunks that are freed [\#530](https://github.com/kokkos/kokkos-kernels/issues/530)


## [2.9.00](https://github.com/kokkos/kokkos-kernels/tree/2.9.00) (2019-06-24)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/2.8.00...2.9.00)

**Implemented enhancements:**

- KokkosBatched: Add specialization for float2, float4 and double4 [\#427](https://github.com/kokkos/kokkos-kernels/pull/427)
- KokkosBatched: Reduce VectorLength (16 to 8) [\#432](https://github.com/kokkos/kokkos-kernels/pull/432)
- KokkosBatched: Remove experimental name space for batched blas [\#371](https://github.com/kokkos/kokkos-kernels/issues/371)
- Capability: Initial sparse triangular solve capability [\#435](https://github.com/kokkos/kokkos-kernels/pull/435)
- Capability: Add support for MAGMA GESV TPL [\#409](https://github.com/kokkos/kokkos-kernels/pull/409)
- cuBLAS: Add CudaUVMSpace specializations for GEMM [\#397](https://github.com/kokkos/kokkos-kernels/issues/397)

**Fixed bugs:**

- Deprecated Code Fixes [\#411](https://github.com/kokkos/kokkos-kernels/issues/411)
- BuildSystem: Compilation error on rzansel [\#401](https://github.com/kokkos/kokkos-kernels/issues/401)

## [2.8.00](https://github.com/kokkos/kokkos-kernels/tree/2.8.00) (2019-02-05)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/2.7.24...2.8.00)

**Implemented enhancements:**

- Capability, Tests: C++14 Support and Testing [\#351](https://github.com/kokkos/kokkos-kernels/issues/351)
- Capability: Batched getrs [\#332](https://github.com/kokkos/kokkos-kernels/issues/332)
- More Kernel Labels for KokkosBlas [\#239](https://github.com/kokkos/kokkos-kernels/issues/239)
- Name all parallel kernels and regions [\#124](https://github.com/kokkos/kokkos-kernels/issues/124)

**Fixed bugs:**

- BLAS TPL: BLAS underscore mangling [\#369](https://github.com/kokkos/kokkos-kernels/issues/369)
- BLAS TPL, Complex: Promotion 2.7.24 broke MV unit tests in Tpetra with complex types [\#360](https://github.com/kokkos/kokkos-kernels/issues/360)
- GEMM: GEMM uses wrong function for computing shared memory allocation size [\#368](https://github.com/kokkos/kokkos-kernels/issues/368)
- BuildSystem: BLAS TPL macro not properly enabled with MKL BLAS [\#347](https://github.com/kokkos/kokkos-kernels/issues/347)
- BuildSystem: make clean - errors [\#353](https://github.com/kokkos/kokkos-kernels/issues/353)
- Compiler Workaround: Internal compiler error in KokkosBatched::Experimental::TeamGemm [\#349](https://github.com/kokkos/kokkos-kernels/issues/349)
- KokkosBlas: Some KokkosBlas kernels assume default execution space [\#14](https://github.com/kokkos/kokkos-kernels/issues/14)

## [2.7.24](https://github.com/kokkos/kokkos-kernels/tree/2.7.24) (2018-11-04)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/2.7.00...2.7.24)

**Implemented enhancements:**

- Enhance test\_all\_sandia script to set scalar and ordinal types [\#315](https://github.com/kokkos/kokkos-kernels/issues/315)
- Batched getri need [\#305](https://github.com/kokkos/kokkos-kernels/issues/305)
- Deterministic Coloring [\#271](https://github.com/kokkos/kokkos-kernels/issues/271)
- MKL - guard minor version for MKL v. 18 [\#268](https://github.com/kokkos/kokkos-kernels/issues/268)
- TPL Support for all BLAS functions using CuBLAS [\#247](https://github.com/kokkos/kokkos-kernels/issues/247)
- Add L1 variant to multithreaded Gauss-Seidel [\#240](https://github.com/kokkos/kokkos-kernels/issues/240)
- Multithreaded Gauss-Seidel does not support damping [\#221](https://github.com/kokkos/kokkos-kernels/issues/221)
- Guard 1-phase SpGEMM in Intel MKL  [\#217](https://github.com/kokkos/kokkos-kernels/issues/217)
- generate makefile with-spaces option  [\#98](https://github.com/kokkos/kokkos-kernels/issues/98)
- Add MKL version check [\#7](https://github.com/kokkos/kokkos-kernels/issues/7)

**Fixed bugs:**

- Perf test failures w/ just CUDA enabled [\#257](https://github.com/kokkos/kokkos-kernels/issues/257)
- Wrong signature for axpy blas functions [\#329](https://github.com/kokkos/kokkos-kernels/issues/329)
- Failing unit tests with float - unit test error checking issue [\#322](https://github.com/kokkos/kokkos-kernels/issues/322)
- cuda.graph\_graph\_color\* COLORING\_VBD test failures with cuda/9.2 + gcc/7.2 on White [\#317](https://github.com/kokkos/kokkos-kernels/issues/317)
- KokkosBatched::Experimental::SIMD\<T\> does not build with T=complex\<float\> [\#316](https://github.com/kokkos/kokkos-kernels/issues/316)
- simple test program fails using 3rdparty Eigen library [\#309](https://github.com/kokkos/kokkos-kernels/issues/309)
- KokkosBlas::dot is broken for complex, due to incorrect assumptions about Fortran ABI [\#307](https://github.com/kokkos/kokkos-kernels/issues/307)
- strides bug in kokkos tpl interface.  [\#292](https://github.com/kokkos/kokkos-kernels/issues/292)
- Failing spgemm unit test with MKL [\#289](https://github.com/kokkos/kokkos-kernels/issues/289)
- Fix the block\_pcg perf-test  when offsets are size\_t [\#287](https://github.com/kokkos/kokkos-kernels/issues/287)
- spotcheck warnings from kokkos  [\#284](https://github.com/kokkos/kokkos-kernels/issues/284)
- Linking error in tpl things [\#282](https://github.com/kokkos/kokkos-kernels/issues/282)
- Build failure with clang 3.9.0 [\#281](https://github.com/kokkos/kokkos-kernels/issues/281)
- CMake modification for TPLs. [\#276](https://github.com/kokkos/kokkos-kernels/issues/276)
- KokkosBatched warnings [\#259](https://github.com/kokkos/kokkos-kernels/issues/259)
- KokkosBatched contraction length bug [\#258](https://github.com/kokkos/kokkos-kernels/issues/258)
- Small error in KokkosBatched\_Gemm\_Serial\_Imp.hpp with SerialGemm\<Trans::Transpose,\*,\*\> [\#147](https://github.com/kokkos/kokkos-kernels/issues/147)

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
