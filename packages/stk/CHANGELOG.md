# CHANGELOG

5.27.3    (STK_VERSION 5270300) 12/16/25
  stk_mesh: initial GPU Mesh mod: NgpMesh::batch_change_entity_parts
            and NgpMesh::update_bulk_data
  stk_mesh: Fixed behavior of multi-state field rotation. Now gives
            correct behavior even for persistent copies of NgpField
            and FieldData instances.
  stk_mesh: Field-BLAS: deprecate field_eamin and field_eamax
  stk_util: For MPI operations that use device-aware MPI, such as the
            stk-mesh parallel_sum and friends, an environment variable
            allows users to explicitly override and use host MPI:
            $ export STK_USE_HOST_MPI_OVERRIDE=true
  all tests: fixes related to googletest being removed from Trilinos
            Users need to separately install googletest
            Use -DGTest_DIR=/path/to/googletest-installation
            Also use -DTPL_ENABLE_gtest
  all:      Requiring/using c++20

5.27.2-02    (STK_VERSION 5270202) 11/14/25
  stk_util: deprecate STK_FUNCTION/STK_INLINE_FUNCTION; prefer direct usage
            of KOKKOS_FUNCTION/KOKKOS_INLINE_FUNCTION
  stk_mesh: Change default access tag for Field::data<>() to ReadOnly.
            Also misc bug fixes related to shell-edge connectivity,
            multi-state field initialization, etc.

5.27.2-01    (STK_VERSION 5270201) 11/05/25
  stk_doc_tests: fix build error when Seacas/Ioss is on but MPI is off
  stk_mesh: add GPU-compatible APIs for parallel_sum and friends.
  stk_util: add Default constructor for MPITag.

5.27.2       (STK_VERSION 5270200) 11/04/25
  stk_expreval: fix Nvidia/cuda compile warning: "calling a constexpr __host__ function("min")
           from a __host__ __device__ function("time_space_normal") is not allowed. "
  stk_mesh: restore ability to build with -DSTK_USE_DEVICE_MESH on non-GPU builds
            (very useful for debugging)

5.27.1-01    (STK_VERSION 5270101) 10/24/25
  stk_unit_test_utils: Fix missing-header build-error encountered by Albany

5.27.1    (STK_VERSION 5270100) 10/21/25
  stk_util: deprecate unused/unimplemented 'malloc_used()' and also the
            HeapAlloc specialization of stk::diag::MetricTraits.
  stk_mesh: Added ability to rename Part (MetaData::rename)

5.25.7    (STK_VERSION 5250700) 9/15/2025
  stk_util: add missing include of <cstdlib> in stk_util/util/AlignedAllocator.hpp
            add missing include of <cstddef> in stk_util/parallel/CommBuffer.hpp
            fix error in stk_expreval/Function.hpp, replacing DBL_MIN with
            std::numeric_limits<double>::min()
  stk_mesh: add optional (off by default) STK_UNIFIED_MEMORY which allows for using
            the same memory allocation, host and device, for field data. This is
            targeted at the ATS-4 platform.

5.25.6-01 (STK_VERSION 5250601) 8/25/2025
  stk_mesh: mostly minor refinements to the new Field Data API
            most stk sub-modules now converted to new Field APIs
            removed deprecated header stk_mesh/base/GetBuckets.hpp
  stk_util: deleted unused stk_util/parallel/DebugTool.hpp
            deleted unused stk_util/util/Pool.hpp

5.25.6 (STK_VERSION 5250600) 7/14/2025
  stk_mesh: Fix Release/Debug hybrid build inconsistency in new
  Field data access.

5.25.4 (STK_VERSION 5250400) 7/9/2025
  stk_tools: new disconnect-elem capability
  stk_mesh: bug-fix for late-created multi-state field, NgpField
  stk_mesh: refinements to new Field Data Access APIs

5.25.3 (STK_VERSION 5250300) 6/10/2025
  stk_unit_test_utils: fix cmake error, some files that depend on
                      stk_search_util needed to be excluded
                      if stk_search_util is disabled.

5.25.2 (STK_VERSION 5250200) 6/6/2025
  stk_transfer_util: fix Lapack detection (only call find_package
                     if TPL_LAPACK_LIBRARIES is not already set.)
  stk_mesh: added new Field data access APIs
  stk_transfer_util: added SimpleTransfer and related classes

5.25.1-02 (STK_VERSION 5250102) 5/27/2025
  stk_mesh: add overload of for_each_entity_run(NgpMesh.. takes bucketIds
  stk_mesh: misc cleanups, remove un-needed structure (fixed-elem-conn)
            from Bucket, remove fmwk-aux-relation from BulkData.
  stk_transfer_util: fix lapack dependency issue in cmake

5.25.1-01 (STK_VERSION 5250101) 5/15/2025
  stk_unit_test_utils: fix build error if STKSearchUtil is disabled.

5.25.1 (STK_VERSION 5250100) 5/05/2025
  stk_util: fix command-line-parser bug with partial flag matching
  stk_tools: stk_block_extractor now preserves case in exodus variable/field names
             stk_block_extractor also now preserves long name lengths
  stk_balance: stk_balance now preserves name lengths
  stk_mesh: Calling code must now initialize Kokkos before creating STK Mesh
            and must not finalize Kokkos until after destroying STK Mesh.
  stk_util: convenience functions stk::initialize and stk::finalize added
            in stk_util/parallel/Parallel.hpp. These functions can be used to
            initialize/finalize Kokkos and MPI.
  stk_mesh: NgpMesh now has local_ids
  stk_mesh: fix MacOS build error for BulkData::declare_entities
  stk_mesh: NGP field-data is bucketized
  stk_mesh: parallel_sum_including_ghosts can work on device or host

5.23.8-03 (STK_VERSION 5230803) 4/22/2025
  stk_tools: fix compiler error in pmesh lib when SEACASNemesis
             is not enabled.
  stk_transfer: Change from libstk_transfer.a to libstk_transfer_lib.a
                In preparation for the coming-soon stk_transfer executable.
                This is consistent with stk_balance_lib and stk_balance exe
                and avoids duplicate target names.

5.23.8-02 (STK_VERSION 5230802) 4/10/2025
  stk_util: CommSparse can switch underlying comm scheme to pre-post recvs
            instead of the default which is sends and probes. This can be
            set at run-time with an environment variable:
                export STK_UNKNOWN_PATTERN_EXCHANGER=Prepost
  stk_mesh: Ghost comm info is now symmetric. This produces a change in
            the procs returned from BulkData::comm_procs(entity, procs). Now
            the ghost-receiver procs know about each other. Previously they
            only knew about the entity owner.

5.23.8-01 (STK_VERSION 5230801) 3/27/2025
  stk_util: Remove unused diag/Resource2.h, diag/String.hpp
  stk_util: Fix size_t issue by including <cstddef> in parallel/ReceiveCounter.hpp
  stk_util: Fix FP err-handling issue for MacOS
  stk_balance: Deprecate stk_balance_m2n executable. (Functionality is available
               in stk_balance executable.)

5.23.8 (STK_VERSION 5230800) 3/10/2025
  stk_mesh: Contains reversion of bucketized NGP fields.
  stk_topology: Fix incorrect deprecation macro.

5.23.5 (STK_VERSION 5230500) 2/11/2025
  stk_mesh: fix calls to member template 'volatile_fast_shared_comm_map' to use '.template ...' syntax
  stk_util: make CommSparse work with messages > 2 GB with newer coupling versions

5.23.2 (STK_VERSION 5230200) 12/11/2024
  misc fixes for AMD/ROCm (ATS-4)
  stk_mesh: speedup for device-field multi-state rotation
            reduce stacksize (sizeof(DeviceMesh)) from ~2900 to ~470
  stk_search: misc fixes
  stk_io: add query for existence of fields on database

5.21.6-1 (STK_VERSION 5210601) 10/31/2024
  stk_mesh, stk_search: more fixes for HIP unified and Cuda no-uvm builds

5.21.6 (STK_VERSION 5210600) 10/25/2024
  stk_search: fix build-error (instantiation error for morton_lbvh_search) for gcc 13.2
  stk_util: added parallel/OutputStreams.hpp
    - which includes the functions outputP0(), output(), output_flush(), set_outputP0(..), reset_default_output_streams().

5.21.5-2 (STK_VERSION 5210502) 10/07/2024
  stk_search: Fixed HIP sort error.
  stk_mesh: add multi-field NGP-FieldBLAS field_fill

5.21.5-1 (STK_VERSION 5210501) 9/27/2024
  stk_mesh: deprecate BulkData::relation_exist

5.21.5 (STK_VERSION 5210500) 9/25/2024
   general: Fixed MI300A unified memory build errors.
   stk_search: Turned off sorted results by default.

5.21.4-1 (STK_VERSION 5210401) 9/04/2024
   Fix cmake configuration errors that occurred on AMD MI300A platform

5.21.4 (STK_VERSION 5210400) 8/29/2024
   minor fixes, no signficant API changes

5.21.3-1 (STK_VERSION 5210301) 8/19/2024
   stk_mesh: fix ~65K limitation on per-bucket size of upward-connectivity tables
    - This is an implementation detail that only arises if user sets large bucket
      capacities such that 'entities-per-bucket'*'avg-upward-connections-per-entity' > 65K.
      By default buckets are capped at 512 and this limit won't be reached in practice.
      - Was using 16bit index type, now defaults to 32bit index type, but is settable
        via cmake-settable option 'STK_ENABLE_16BIT_UPWARDCONN_INDEX_TYPE=ON/OFF'

5.21.3 (STK_VERSION 5210300) 8/12/2024
   general: compile-warnings/errors fixed for gcc 12 and arm
   stk_mesh: BulkData::change_entity_owner now returns a bool
    - indicates whether any entities actually changed owner.
   stk_mesh: continue "simple field" transition
    - remove some previously-deprecated functions/behaviors
    - 'use_simple_field()' calls now unnecessary, will be deprecated in future
   stk_mesh: deleted previously-deprecated types/methods
    - DeviceMeshIndex and usage in DeviceMesh/HostMesh
    - NgpFieldBase::rotate_multistate_data()
    - DeviceMesh::host_get_entity(..)
   stk_search: changing template parameters for input Views
    - now less restrictive on how users' views were declared

5.21.1 (STK_VERSION 5210100) 6/17/2024
   stk_mesh: Deprecated BulkData::find_permutation and BulkData::check_permutation
             (replaced by free-functions in FindPermutation.hpp)
   stk_mesh: Added destroy_relations free-function in DestroyRelations.hpp
   stk_mesh: Added overloads of field_fill, field_copy, field_axpbyz in NgpFieldBLAS.hpp

5.19.4 (STK_VERSION 5190401) 5/29/2024
   stk_search: fixed bug in morton: (accessing device view on host)
   stk_search: fixed implementations to respect execution-space
   stk_mesh: change default bucket capacity from 512 to 16 (now grows as needed)

5.19.4 (STK_VERSION 5190400) 5/18/2024
   Added field_fill with Selector in stk_mesh/base/NgpFieldBLAS.hpp
   Added field_copy with and without Selector in stk_mesh/base/NgpFieldBLAS.hpp
   More work on stk-search coarse_search APIs and implementation correctness.
   Minor cleanups in cmake support.
   Fixed stk_search KDTREE OpenMP build

5.19.2 (STK_VERSION 5190201) 4/22/2024
   Fixed broken stk_search install targets

5.19.2 (STK_VERSION 5190200) 4/17/2024
   Added bool argument to stk::mesh::BulkData::update_field_data_states()
   Initial creation of STK_VERSION macro

