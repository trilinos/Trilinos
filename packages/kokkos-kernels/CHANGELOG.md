# Change Log

## [3.6.01](https://github.com/kokkos/kokkos-kernels/tree/3.6.01) (2022-05-23)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/3.6.00...3.6.01)

### Bug Fixes and Improvements:

- Improve spiluk numeric phase to avoid race conditions and processing in chunks [\#1390](https://github.com/kokkos/kokkos-kernels/pull/1390)
- Improve sptrsv symbolic phase performance (level scheduling) [\#1380](https://github.com/kokkos/kokkos-kernels/pull/1380)
- Restore BLAS-1 MV paths for 1 column [\#1354](https://github.com/kokkos/kokkos-kernels/pull/1354)
- Fix check that view has const type [\#1370](https://github.com/kokkos/kokkos-kernels/pull/1370)
- Fix check that view has const type part 2 [\#1394](https://github.com/kokkos/kokkos-kernels/pull/1394)

## [3.6.00](https://github.com/kokkos/kokkos-kernels/tree/3.6.00) (2022-02-18)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/3.5.00...3.6.00)

### Features: 

#### Batched Sparse Linear algebra
- Kokkos Kernels is adding a new component to the library: batched sparse linear algebra.
- Similarly to the current dense batched algorithms, the new algorithms are called from
- the GPU and provide Team and TeamVector level of parallelism, SpMV also provides a Serial
- call on GPU.

- Add Batched CG and Batched GMRES [\#1155](https://github.com/kokkos/kokkos-kernels/pull/1155)
- Add Jacobi Batched preconditioner [\#1219](https://github.com/kokkos/kokkos-kernels/pull/1219)

#### Bsr and Tensor core algorithm for sparse linear algebra
- After introducing the BsrMatrix in release 3.5.0 new algorithms are now supporting this format.
- For release 3.6.0 we are adding matrix-vector (matvec) multiplication and Gauss-Seidel as well as an
- implementation of matvec that leverages tensor cores on Nvidia GPUs. More kernels are expected to
- support the Bsr format in future releases.

- Add Spmv for BsrMatrix [\#1255](https://github.com/kokkos/kokkos-kernels/pull/1255)
- Add BLAS to SpMV operations for BsrMatrix [\#1297](https://github.com/kokkos/kokkos-kernels/pull/1297)
- BSR format support in block Gauss-Seidel [\#1232](https://github.com/kokkos/kokkos-kernels/pull/1232)
- Experimental tensor-core SpMV for BsrMatrix [\#1090](https://github.com/kokkos/kokkos-kernels/pull/1090)

#### Improved AMD math libraries support
- rocBLAS and rocSPARSE TPLs are now officially supported, they can be enabled at configure time.
- Initial kernels that can call rocBLAS are GEMV, GEMM, IAMAX and SCAL, while rocSPARSE can be
- called for matrix-vector multiplication. Further support for TPL calls can be requested on slack
- and by GitHub issues.

- Tpl rocBLAS and rocSPARSE [\#1153](https://github.com/kokkos/kokkos-kernels/pull/1153)
- Add rocBLAS GEMV wrapper [\#1201](https://github.com/kokkos/kokkos-kernels/pull/1201)
- Add rocBLAS wrappers for GEMM, IAMAX, and SCAL [\#1230](https://github.com/kokkos/kokkos-kernels/pull/1230)
- SpMV: adding support for rocSPARSE TPL [\#1221](https://github.com/kokkos/kokkos-kernels/pull/1221)

#### Additional new features
- bhalf: Unit test Batched GEMM [\#1251](https://github.com/kokkos/kokkos-kernels/pull/1251) 
-   and demostrate GMRES example convergence with bhalf_t (https://github.com/kokkos/kokkos-kernels/pull/1300)
- Stream interface: adding stream support in GEMV and GEMM [\#1131](https://github.com/kokkos/kokkos-kernels/pull/1131)
- Improve double buffering batched gemm performance [\#1217](https://github.com/kokkos/kokkos-kernels/pull/1217)
- Allow choosing coloring algorithm in multicolor GS [\#1199](https://github.com/kokkos/kokkos-kernels/pull/1199)
- Batched: Add armpl dgemm support [\#1256](https://github.com/kokkos/kokkos-kernels/pull/1256)

### Deprecations:
- Deprecation warning: SpaceAccessibility move out of impl, see #1140 [\#1141](https://github.com/kokkos/kokkos-kernels/pull/1141)

### Backends and Archs Enhancements:

#### SYCL:
- Full Blas support on SYCL [\#1270](https://github.com/kokkos/kokkos-kernels/pull/1270)
- Get sparse tests enabled and working for SYCL [\#1269](https://github.com/kokkos/kokkos-kernels/pull/1269)
- Changes to make graph run on SYCL [\#1268](https://github.com/kokkos/kokkos-kernels/pull/1268)
- Allow querying free/total memory for SYCL [\#1225](https://github.com/kokkos/kokkos-kernels/pull/1225)
- Use KOKKOS_IMPL_DO_NOT_USE_PRINTF instead of printf in kernels [\#1162](https://github.com/kokkos/kokkos-kernels/pull/1162)

#### HIP:
- Work around hipcc size_t/int division with remainder bug [\#1262](https://github.com/kokkos/kokkos-kernels/pull/1262)

#### Other Improvements:
- Replace std::abs with ArithTraits::abs [\#1312](https://github.com/kokkos/kokkos-kernels/pull/1312)
- Batched/dense: Add Gemm_DblBuf LayoutLeft operator [\#1299](https://github.com/kokkos/kokkos-kernels/pull/1299)
- KokkosKernels: adding variable that returns version as a single number [\#1295](https://github.com/kokkos/kokkos-kernels/pull/1295)
- Add KOKKOSKERNELS_FORCE_SIMD macro (Fix #1040) [\#1290](https://github.com/kokkos/kokkos-kernels/pull/1290)
- Rename KOKKOS_IF_{HOST,DEVICE} -> KOKKOS_IF_ON_{HOST,DEVICE} [\#1278](https://github.com/kokkos/kokkos-kernels/pull/1278)
- Algo::Level{2,3}::Blocked::mb() [\#1265](https://github.com/kokkos/kokkos-kernels/pull/1265)
- Batched: Use SerialOpt2 for 33 to 39 square matrices [\#1261](https://github.com/kokkos/kokkos-kernels/pull/1261)
- Prune extra dependencies [\#1241](https://github.com/kokkos/kokkos-kernels/pull/1241)
- Improve double buffering batched gemm perf for matrix sizes >64x64 [\#1239](https://github.com/kokkos/kokkos-kernels/pull/1239)
- Improve graph color perf test [\#1229](https://github.com/kokkos/kokkos-kernels/pull/1229)
- Add custom implementation for strcasecmp [\#1227](https://github.com/kokkos/kokkos-kernels/pull/1227)
- Replace __restrict__ with KOKKOS_RESTRICT [\#1223](https://github.com/kokkos/kokkos-kernels/pull/1223)
- Replace array reductions in BLAS-1 MV reductions [\#1204](https://github.com/kokkos/kokkos-kernels/pull/1204)
- Update MIS-2 and aggregation [\#1143](https://github.com/kokkos/kokkos-kernels/pull/1143)
- perf_test/blas/blas3: Update SHAs for benchmarking [\#1139](https://github.com/kokkos/kokkos-kernels/pull/1139)

### Implemented enhancements BuildSystem
- Bump ROCm version 4.2 -> 4.5 in nightly Jenkins CI build [\#1279](https://github.com/kokkos/kokkos-kernels/pull/1279)
- scripts/cm_test_all_sandia: Add A64FX ci checks [\#1276](https://github.com/kokkos/kokkos-kernels/pull/1276)
- github/workflows: Add osx CI [\#1254](https://github.com/kokkos/kokkos-kernels/pull/1254)
- Update SYCL compiler version in CI [\#1247](https://github.com/kokkos/kokkos-kernels/pull/1247)
- Do not set Kokkos variables when exporting CMake configuration [\#1236](https://github.com/kokkos/kokkos-kernels/pull/1236)
- Add nightly CI check for SYCL [\#1190](https://github.com/kokkos/kokkos-kernels/pull/1190)
- Update cmake minimum version to 3.16 [\#866](https://github.com/kokkos/kokkos-kernels/pull/866)

### Incompatibilities:
- Kokkos::Impl: removing a few more instances of throw_runtime_exception [\#1320](https://github.com/kokkos/kokkos-kernels/pull/1320)
- Remove Kokkos::Impl::throw_runtime_exception from Kokkos Kernels [\#1294](https://github.com/kokkos/kokkos-kernels/pull/1294)
- Remove unused memory space utility [\#1283](https://github.com/kokkos/kokkos-kernels/pull/1283)
- Clean up Kokkos header includes [\#1282](https://github.com/kokkos/kokkos-kernels/pull/1282)
- Remove private Kokkos header include (Cuda/Kokkos_Cuda_Half.hpp) [\#1281](https://github.com/kokkos/kokkos-kernels/pull/1281)
- Avoid using #ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_* macro guards [\#1266](https://github.com/kokkos/kokkos-kernels/pull/1266)
- Rename enumerator Impl::Exec_{PTHREADS -> THREADS} [\#1253](https://github.com/kokkos/kokkos-kernels/pull/1253)
- Remove all references to the Kokkos QThreads backend [\#1238](https://github.com/kokkos/kokkos-kernels/pull/1238)
- Replace more occurences of Kokkos::Impl::is_view [\#1234](https://github.com/kokkos/kokkos-kernels/pull/1234)
- Do not use Kokkos::Impl::is_view [\#1214](https://github.com/kokkos/kokkos-kernels/pull/1214)
- Replace Kokkos::Impl::if_c -> std::conditional [\#1213](https://github.com/kokkos/kokkos-kernels/pull/1213)

### Bug Fixes:
- Fix bug in spmv_mv_bsrmatrix() for Ampere GPU arch [\#1315](https://github.com/kokkos/kokkos-kernels/pull/1315)
- Fix std::abs calls for rocBLAS/rocSparse [\#1310](https://github.com/kokkos/kokkos-kernels/pull/1310)
- cast literal 0 to fragment scalar type [\#1307](https://github.com/kokkos/kokkos-kernels/pull/1307)
- Fix 1303: maintain correct #cols on A in twostage [\#1304](https://github.com/kokkos/kokkos-kernels/pull/1304)
- Add dimension checking to generic spmv interface [\#1301](https://github.com/kokkos/kokkos-kernels/pull/1301)
- Add missing barriers to TeamGMRES, fix vector len [\#1285](https://github.com/kokkos/kokkos-kernels/pull/1285)
- Examples: fixing some issues related to type checking [\#1267](https://github.com/kokkos/kokkos-kernels/pull/1267)
- Restrict BsrMatrix specialization for AMPERE and VOLTA to CUDA [\#1242](https://github.com/kokkos/kokkos-kernels/pull/1242)
- Fix compilation errors for multi-vectors in kk_print_1Dview() [\#1231](https://github.com/kokkos/kokkos-kernels/pull/1231)
- src/batched: Fixes #1224 [\#1226](https://github.com/kokkos/kokkos-kernels/pull/1226)
- Fix SpGEMM crashing on empty rows [\#1220](https://github.com/kokkos/kokkos-kernels/pull/1220)
- Fix issue #1212 [\#1218](https://github.com/kokkos/kokkos-kernels/pull/1218)
- example/gmres: Specify half_t namespace [\#1208](https://github.com/kokkos/kokkos-kernels/pull/1208)
- Check that ordinal types are signed [\#1188](https://github.com/kokkos/kokkos-kernels/pull/1188)
- Fixing a couple of small issue with tensor core spmv [\#1185](https://github.com/kokkos/kokkos-kernels/pull/1185)
- Fix #threads setting in pcg for OpenMP [\#1182](https://github.com/kokkos/kokkos-kernels/pull/1182)
- SpMV: fix catch all case to avoid compiler warnings [\#1179](https://github.com/kokkos/kokkos-kernels/pull/1179)
- using namespace should be scoped to prevent name clashes [\#1177](https://github.com/kokkos/kokkos-kernels/pull/1177)
- using namespace should be scoped to prevent name clashes, see issue #1170 [\#1171](https://github.com/kokkos/kokkos-kernels/pull/1171)
- Fix bug with mkl impl of spgemm [\#1167](https://github.com/kokkos/kokkos-kernels/pull/1167)
- Add missing $ to KOKKOS_HAS_TRILINOS in sparse_sptrsv_superlu check [\#1160](https://github.com/kokkos/kokkos-kernels/pull/1160)
- Small fixes to spgemm, and plug gaps in testing [\#1159](https://github.com/kokkos/kokkos-kernels/pull/1159)
- SpMV: mismatch in #ifdef check and kernel specialization [\#1151](https://github.com/kokkos/kokkos-kernels/pull/1151)
- Fix values dimension for block sparse matrices [\#1147](https://github.com/kokkos/kokkos-kernels/pull/1147)

## [3.5.00](https://github.com/kokkos/kokkos-kernels/tree/3.5.00) (2021-10-19)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/3.4.01...3.5.00)

**Features:**
- Batched serial SVD [\#1107](https://github.com/kokkos/kokkos-kernels/pull/1107)
- Batched: Add BatchedDblBufGemm [\#1095](https://github.com/kokkos/kokkos-kernels/pull/1095)
- feature/gemv rps test -- RAJAPerf Suite Version of the BLAS2 GEMV Test [\#1085](https://github.com/kokkos/kokkos-kernels/pull/1085)
- Add new bsrmatrix [\#1077](https://github.com/kokkos/kokkos-kernels/pull/1077)
- Adding Kokkos GMRES example [\#1028](https://github.com/kokkos/kokkos-kernels/pull/1028)
- Add fast two-level mode N GEMV (#926) [\#939](https://github.com/kokkos/kokkos-kernels/pull/939)
- Batched: Add BatchedGemm interface [\#935](https://github.com/kokkos/kokkos-kernels/pull/935)
- OpenMPTarget: adding ETI and CMake logic for OpenMPTarget backend [\#886](https://github.com/kokkos/kokkos-kernels/pull/886)

**Implemented enhancements Algorithms and Archs:**
- Use float as accumulator for GEMV on half_t (Fix #1081) [\#1082](https://github.com/kokkos/kokkos-kernels/pull/1082)
- Supernodal SpTRSV: add option to use MAGMA TPL for TRTRI [\#1069](https://github.com/kokkos/kokkos-kernels/pull/1069)
- Updates for running GMRES example with half precision [\#1067](https://github.com/kokkos/kokkos-kernels/pull/1067)
- src/blas/impl: Explicitly cast to LHS type for ax [\#1073](https://github.com/kokkos/kokkos-kernels/pull/1073)
- Update BatchedGemm interface to match design proposal [\#1054](https://github.com/kokkos/kokkos-kernels/pull/1054)
- Move dot-based GEMM out of TPL CUBLAS [\#1050](https://github.com/kokkos/kokkos-kernels/pull/1050)
- Adding ArmPL option to spmv perf_test [\#1038](https://github.com/kokkos/kokkos-kernels/pull/1038)
- Add (right) preconditioning to GMRES [\#1078](https://github.com/kokkos/kokkos-kernels/pull/1078)
- Supernodal SpTRSV: perform TRMM only if TPL CuBLAS is enabled [\#1027](https://github.com/kokkos/kokkos-kernels/pull/1027)
- Supernodal SpTRSV: support SuperLU version < 5 [\#1012](https://github.com/kokkos/kokkos-kernels/pull/1012)
- perf_test/blas/blas3: Add dgemm armpl experiment [\#1005](https://github.com/kokkos/kokkos-kernels/pull/1005)
- Supernodal SpTRSV: run TRMM on device for setup [\#983](https://github.com/kokkos/kokkos-kernels/pull/983)
- Merge pull request #951 from vqd8a/move_sort_ifpack2riluk [\#972](https://github.com/kokkos/kokkos-kernels/pull/972)
- Point multicolor GS: faster handling of long/bulk rows [\#993](https://github.com/kokkos/kokkos-kernels/pull/993)
- Make CRS sorting utils work with unmanaged [\#963](https://github.com/kokkos/kokkos-kernels/pull/963)
- Add sort and make sure using host mirror on host memory in kspiluk_symbolic [\#951](https://github.com/kokkos/kokkos-kernels/pull/951)
- GEMM: call GEMV instead in certain cases [\#948](https://github.com/kokkos/kokkos-kernels/pull/948)
- SpAdd performance improvements, better perf test, fix mtx reader columns [\#930](https://github.com/kokkos/kokkos-kernels/pull/930)

**Implemented enhancements BuildSystem:**
- Automate documentation generation [\#1116](https://github.com/kokkos/kokkos-kernels/pull/1116)
- Move the batched dense files to specific directories [\#1098](https://github.com/kokkos/kokkos-kernels/pull/1098)
- cmake: Update SUPERLU tpl option for Tribits [\#1066](https://github.com/kokkos/kokkos-kernels/pull/1066)
- cmake/Modules: Allow user to use MAGMA_DIR from env [\#1007](https://github.com/kokkos/kokkos-kernels/pull/1007)
- Supernodal SpTRSV: update TPLs requirements [\#997](https://github.com/kokkos/kokkos-kernels/pull/997)
- cmake: Add MAGMA TPL support [\#982](https://github.com/kokkos/kokkos-kernels/pull/982)
- Host only macro: adding macro to check for any device backend [\#940](https://github.com/kokkos/kokkos-kernels/pull/940)
- Prevent redundant spmv kernel instantiations (reduce library size) [\#937](https://github.com/kokkos/kokkos-kernels/pull/937)
- unit-test: refactor infrastructure to remove most *.cpp [\#906](https://github.com/kokkos/kokkos-kernels/pull/906)

**Implemented enhancements Other:**
- Allow reading integer mtx files into floating-point matrices [\#1100](https://github.com/kokkos/kokkos-kernels/pull/1100)
- Warnings: remove -Wunused-parameter warnings in Kokkos Kernels [\#962](https://github.com/kokkos/kokkos-kernels/pull/962)
- Clean up CrsMatrix raw pointer constructor [\#949](https://github.com/kokkos/kokkos-kernels/pull/949)
- unit_test/batched: Remove *_half fns from gemm unit tests [\#943](https://github.com/kokkos/kokkos-kernels/pull/943)
- Move sorting functionality out of Impl:: [\#932](https://github.com/kokkos/kokkos-kernels/pull/932)

**Incompatibilities:**
- Deprecation warning: SpaceAccessibility move out of impl [\#1141](https://github.com/kokkos/kokkos-kernels/pull/1141)
- Rename CUDA_SAFE_CALL to KOKKOS_IMPL_CUDA_SAFE_CALL [\#1130](https://github.com/kokkos/kokkos-kernels/pull/1130)
- Workaround error with intel [\#1128](https://github.com/kokkos/kokkos-kernels/pull/1128)
- gmres: disable examples for builds with ibm/xl [\#1123](https://github.com/kokkos/kokkos-kernels/pull/1123)
- CrsMatrix: deprecate constructor without ncols input [\#1115](https://github.com/kokkos/kokkos-kernels/pull/1115)
- perf_test/blas/blas3: Disable simd verify for cuda/10.2.2 [\#1093](https://github.com/kokkos/kokkos-kernels/pull/1093)
- Replace impl/Kokkos_Timer.hpp includes with Kokkos_Timer.hpp [\#1074](https://github.com/kokkos/kokkos-kernels/pull/1074)
- Remove deprecated ViewAllocateWithoutInitializing [\#1058](https://github.com/kokkos/kokkos-kernels/pull/1058)
- src/sparse: spadd resolve deprecation warnings [\#1053](https://github.com/kokkos/kokkos-kernels/pull/1053)
- Give full namespace path for D2 coloring [\#999](https://github.com/kokkos/kokkos-kernels/pull/999)
- Fix -Werror=deprecated errors with c++20 standard [\#964](https://github.com/kokkos/kokkos-kernels/pull/964)
- Deprecation: a deprecated function is called in the SpADD perf_test [\#954](https://github.com/kokkos/kokkos-kernels/pull/954)

**Enabled tests:**
- HIP: enabling all unit tests [\#968](https://github.com/kokkos/kokkos-kernels/pull/968)
- Fix build and add CI coverage for LayoutLeft=OFF [\#965](https://github.com/kokkos/kokkos-kernels/pull/965)
- Enable SYCL tests [\#927](https://github.com/kokkos/kokkos-kernels/pull/927)
- Fixup HIP nightly builds [\#907](https://github.com/kokkos/kokkos-kernels/pull/907)

**Fixed Bugs:**
- Fix SpGEMM for Nvidia Turing/Ampere [\#1118](https://github.com/kokkos/kokkos-kernels/pull/1118)
- Fix #1111: spmv tpl instantiations [\#1112](https://github.com/kokkos/kokkos-kernels/pull/1112)
- Fix C's numCols in spadd simplified interface [\#1102](https://github.com/kokkos/kokkos-kernels/pull/1102)
- Fix #1089 (failing batched UTV tests) [\#1096](https://github.com/kokkos/kokkos-kernels/pull/1096)
- Blas GEMM: fix early exit logic, see issue #1088 [\#1091](https://github.com/kokkos/kokkos-kernels/pull/1091)
- Fix #1048: handle mode C spmv correctly in serial/openmp [\#1084](https://github.com/kokkos/kokkos-kernels/pull/1084)
- src/batched: Fix multiple definitions of singleton [\#1072](https://github.com/kokkos/kokkos-kernels/pull/1072)
- Fix host accessing View in non-host space [\#1057](https://github.com/kokkos/kokkos-kernels/pull/1057)
- Fix559: Intel 18 has trouble with pointer in ternary expr [\#1042](https://github.com/kokkos/kokkos-kernels/pull/1042)
- Work around team size AUTO issue on kepler [\#1020](https://github.com/kokkos/kokkos-kernels/pull/1020)
- Supernodal SpTrsv: fix out-of-bound error [\#1019](https://github.com/kokkos/kokkos-kernels/pull/1019)
- Some fixes for MAGMA TPL and gesv [\#1008](https://github.com/kokkos/kokkos-kernels/pull/1008)
- Merge pull request #981 from Tech-XCorp/4005-winllvmbuild [\#984](https://github.com/kokkos/kokkos-kernels/pull/984)
- This is a PR for 4005 vs2019build, which fixes a few things on Windows [\#981](https://github.com/kokkos/kokkos-kernels/pull/981)
- Fix build for no-ETI build [\#977](https://github.com/kokkos/kokkos-kernels/pull/977)
- Fix invalid mem accesses in new GEMV kernel [\#961](https://github.com/kokkos/kokkos-kernels/pull/961)
- Kokkos_ArithTraits.hpp: Fix isInf and isNan with complex types [\#936](https://github.com/kokkos/kokkos-kernels/pull/936)

## [3.4.01](https://github.com/kokkos/kokkos-kernels/tree/3.4.01) (2021-05-19)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/3.4.00...3.4.01)

**Fixed Bugs:**
- Windows: Fixes for Windows [\#981](https://github.com/kokkos/kokkos-kernels/pull/981)
- Sycl: ArithTraits fixes for Sycl [\#959](https://github.com/kokkos/kokkos-kernels/pull/959)
- Sparse: Added code to allow KokkosKernels coloring to accept partial colorings [\#938](https://github.com/kokkos/kokkos-kernels/pull/938)
- Sparse: Include sorting within spiluk [\#972](https://github.com/kokkos/kokkos-kernels/pull/972)
- Sparse: Fix CrsMatrix raw pointer constructor [\#971](https://github.com/kokkos/kokkos-kernels/pull/971)
- Sparse: Fix spmv Serial beta==-1 code path [\#947](https://github.com/kokkos/kokkos-kernels/pull/947)

## [3.4.00](https://github.com/kokkos/kokkos-kernels/tree/3.4.00) (2021-04-25)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/3.3.01...3.4.00)

**Features:**
- SYCL: adding ETI and CMake logic for SYCL backend [\#924](https://github.com/kokkos/kokkos-kernels/pull/924)

**Implemented enhancements Algorithms and Archs:**
- Two-stage GS: add damping factors [\#921](https://github.com/kokkos/kokkos-kernels/pull/921)
- Supernodal SpTRSV, improve symbolic performance [\#899](https://github.com/kokkos/kokkos-kernels/pull/899)
- Add MKL SpMV wrapper [\#895](https://github.com/kokkos/kokkos-kernels/pull/895)
- Serial code path for spmv [\#893](https://github.com/kokkos/kokkos-kernels/pull/893)

**Implemented enhancements BuildSystem:**
- Cmake: Update ArmPL support [\#901](https://github.com/kokkos/kokkos-kernels/pull/901)
- Cmake: Add ARMPL TPL support [\#880](https://github.com/kokkos/kokkos-kernels/pull/880)
- IntelClang guarding __assume_aligned with !defined(__clang__) [\#878](https://github.com/kokkos/kokkos-kernels/pull/878)

**Implemented enhancements Other:**
- Add static_assert/throw in batched eigendecomp [\#931](https://github.com/kokkos/kokkos-kernels/pull/931)
- Workaround using new/delete in kernel code [\#925](https://github.com/kokkos/kokkos-kernels/pull/925)
- Blas perf_test updates [\#892](https://github.com/kokkos/kokkos-kernels/pull/892)

**Fixed bugs:**
- Fix ctor CrsMat mirror with CrsGraph mirror [\#918](https://github.com/kokkos/kokkos-kernels/pull/918)
- Fix nrm1, removed cublas nrminf, improved blas tests [\#915](https://github.com/kokkos/kokkos-kernels/pull/915)
- Fix and testing coverage mainly in graph coarsening [\#910](https://github.com/kokkos/kokkos-kernels/pull/910)
- Fix KokkosSparse for nightly test failure [\#898](https://github.com/kokkos/kokkos-kernels/pull/898)
- Fix view types across ternary operator [\#894](https://github.com/kokkos/kokkos-kernels/pull/894)
- Make work_view_t typedef consistent [\#885](https://github.com/kokkos/kokkos-kernels/pull/885)
- Fix supernodal SpTRSV build with serial+openmp+cuda [\#884](https://github.com/kokkos/kokkos-kernels/pull/884)
- Construct SpGEMM C with correct ncols [\#883](https://github.com/kokkos/kokkos-kernels/pull/883)
- Matrix Converter: fixing issue with deallocation after Kokkos::fininalize [\#882](https://github.com/kokkos/kokkos-kernels/pull/882)
- Fix >1024 team size error in sort_crs_* [\#872](https://github.com/kokkos/kokkos-kernels/pull/872)
- Fixing seg fault with empty matrix in kspiluk [\#871](https://github.com/kokkos/kokkos-kernels/pull/871)

## [3.3.01](https://github.com/kokkos/kokkos-kernels/tree/3.3.01) (2021-01-18)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/3.3.00...3.3.01)

**Fixed Bugs:**
- With CuSparse enabled too many variants of SPMV were instantiated even if not requested. Up to 1GB executable size increase.

## [3.3.00](https://github.com/kokkos/kokkos-kernels/tree/3.3.00) (2020-12-16)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/3.2.01...3.3.00)

**Implemented enhancements:**
- Add permanent RCM reordering interface, and a basic serial implementation [\#854](https://github.com/kokkos/kokkos/pull/#854)
- Half\_t explicit conversions [\#849](https://github.com/kokkos/kokkos/pull/#849)
- Add batched gemm performance tests [\#838](https://github.com/kokkos/kokkos/pull/#838)
- Add HIP support to src and perf\_test [\#828](https://github.com/kokkos/kokkos/pull/#828)
- Factor out coarsening [\#827](https://github.com/kokkos/kokkos/pull/#827)
- Allow enabling/disabling components at configuration time [\#823](https://github.com/kokkos/kokkos/pull/#823)
- HIP: CMake work on tests and ETI  [\#820](https://github.com/kokkos/kokkos/pull/#820)
- HIP: KokkosBatched - hip specialization [\#812](https://github.com/kokkos/kokkos/pull/#812)
- Distance-2 maximal independent set [\#801](https://github.com/kokkos/kokkos/pull/#801)
- Use batched TRTRI & TRMM for Supernode-sptrsv setup [\#797](https://github.com/kokkos/kokkos/pull/#797)
- Initial support for half precision [\#794](https://github.com/kokkos/kokkos/pull/#794)

**Fixed bugs:**
- Fix issue with HIP and Kokkos\_ArithTraits [\#844](https://github.com/kokkos/kokkos/pull/#844)
- HIP: fixing round of issues on AMD [\#840](https://github.com/kokkos/kokkos/pull/#840)
- Throw an exception if BLAS GESV is not enabled [\#837](https://github.com/kokkos/kokkos/pull/#837)
- Fixes -Werror for gcc with c++20 [\#836](https://github.com/kokkos/kokkos/pull/#836)
- Add fallback condition to use spmv\_native when cuSPARSE does not work [\#834](https://github.com/kokkos/kokkos/pull/#834)
- Fix install testing refactor for inline builds [\#811](https://github.com/kokkos/kokkos/pull/#811)
- HIP: fix ArithTraits to support HIP backend [\#809](https://github.com/kokkos/kokkos/pull/#809)
- cuSPARSE 11: fix spgemm and spmv\_struct\_tunning compilation error [\#804](https://github.com/kokkos/kokkos/pull/#804)

**Incompatibilities:**
- Remove pre-3.0 deprecated code [\#825](https://github.com/kokkos/kokkos/pull/#825)

## [3.2.01](https://github.com/kokkos/kokkos-kernels/tree/3.2.01) (2020-11-17)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/3.2.00...3.2.01)

**Fixed bugs:**

- Cpp14 Fixes: [\#790](https://github.com/kokkos/kokkos-kernels/pull/790)

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

## [3.1.01](https://github.com/kokkos/kokkos-kernels/tree/3.1.01) (2020-05-04)
[Full Changelog](https://github.com/kokkos/kokkos-kernels/compare/3.1.00...3.1.01)

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
