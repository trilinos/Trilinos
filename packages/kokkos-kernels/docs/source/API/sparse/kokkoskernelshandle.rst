KokkosKernels::Experimental::KokkosKernelsHandle
################################################

.. toctree::
   :maxdepth: 1
   :hidden:

   handle_get_create_destroy

Defined in header ``KokkosKernels_Handle.hpp``

.. code:: cppkokkos

  template <class size_type_, class lno_t_, class scalar_t_, class ExecutionSpace,
            class TemporaryMemorySpace, class PersistentMemorySpace>
  class KokkosKernelsHandle;

``KokkosKernels::KokkosKernelsHandle`` implements an opaque handle that stores algorithm specific sub-handles, verbosity parameter and performance optimization parameters.

Template Parameters
===================

:size_type\_: type associated with the row map of a matrix and in general with offset and size in memory.

:lno_t\_: type associated with the column indices of a CrsMatrix, it is commonly used to label entries in a set.

:scalar_t\_: type associated with the values stored in a matrix and to perform floating point operations in general.

:ExecutionSpace: type of the execution space where kernels are executed.

:TemporaryMemorySpace: a memory space associated with buffers and temporary allocations in kernels.

:PersistentMemorySpace: a memory space associated with permanent data generated during a kernel, for instance to allocate view that can be returned to the user.

General Member Types
====================

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Member type
     - Definition

   * - HandleExecSpace
     - Default execution space associated with the handle.

   * - HandleTempMemorySpace
     - Memory space typically with faster access for temporary data and buffers (in practice often the same as ``HandlePersistentMemorySpace``)

   * - HandlePersistentMemorySpace
     - Default memory space for data structures allocation

   * - size_type
     - potentially unsigned integer type, often the type used to store row offsets of a CrsMatrix row map.

   * - const_size_type
     - ``const size_type``.

   * - nnz_lno_t
     - signed integer type, often the type used to store column indices in a CrsMatrix or colors of vertices in graph coloring.

   * - const_nnz_lno_t
     - ``const nnz_lno_t``.

   * - nnz_scalar_t
     - floating point type typically used to store values in a CrsMatrix.

   * - const_nnz_scalar_t
     - ``const nnz_scalar_t``.

Sub-handle Member Types
=======================

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Member type
     - Definition

   * - GraphColoringHandleType
     - Type of the associated distance one graph coloring handle.

   * - GraphColorDistance2HandleType
     - Type of the associated distance two graph coloring handle.

   * - GaussSeidelHandleType
     - Type of the associated Gauss-Seidel handle.

   * - PointGaussSeidelHandleType
     - Type of the associated point Gauss-Seidel handle, sub-class of GaussSeidelHandleType.

   * - ClusterGaussSeidelHandleType
     - Type of the associated cluster Gauss-Seidel handle, sub-class of GaussSeidelHandleType.

   * - TwoStageGaussSeidelHandleType
     - Type of the associated two stages Gauss-Seidel handle, sub-class of GaussSeidelHandleType.

   * - TwoStageGaussSeidelSPTRSVHandleType
     - Type of the associated SpTRSV sub-handle of the two stages Gauss-Seidel handle.

   * - SPGEMMHandleType
     - Type of the associated SpGEMM handle.

   * - SPADDHandleType
     - Type of the associated SpADD handle.

   * - SPTRSVHandleType
     - Type of the associated SpTRSV handle.

   * - SPILUKHandleType
     - Type of the associated sparse ILU(K) handle.

   * - PAR_ILUTHandleType
     - Type of the associated parallel ILUt handle.

   * - GMRESHandleType
     - Type of the associated GMRES handle.

View Member Types
=================

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Member type
     - Definition

   * - in_scalar_nnz_view_t
     - View of rank 1 holding ``scalar`` values and allocated in ``HandleTempMemorySpace``.

   * - row_lno_temp_work_view_t
     - alias for ``size_type_temp_work_view_t`` below.

   * - size_type_temp_work_view_t
     - View of rank 1 holding ``size_type`` values and allocated in ``HandleTempMemorySpace``.

   * - row_lno_persistent_work_view_t
     - alias for ``size_type_persistent_work_view_t`` below.

   * - size_type_persistent_work_view_t
     - View of rank 1 holding ``size_type`` values and allocated in ``HandlePersistentMemorySpace``.

   * - row_lno_persistent_work_host_view_t
     - ``host_mirror_type`` type of ``row_lno_persistent_work_host_view_t`` type.

   * - size_type_persistent_work_host_view_t
     - ``host_mirror_type`` type of ``size_type_persistent_work_view_t`` type.

   * - scalar_temp_work_view_t
     - View of rank 1 holding ``scalar`` values and allocated in ``HandleTempMemorySpace``.

   * - scalar_persistent_work_view_t
     - View of rank 1 holding ``scalar`` values and allocated in ``HandlePersistentMemorySpace``.

   * - scalar_persistent_work_view2d_t
     - View of rank 2 holding ``scalar`` values and allocated in ``HandlePersistentMemorySpace``.

   * - nnz_lno_temp_work_view_t
     - View of rank 1 holding ``nnz_lno_t`` values and allocated in ``HandleTempMemorySpace``.

   * - nnz_lno_persistent_work_view_t
     - View of rank 1 holding ``nnz_lno_t`` values and allocated in ``HandlePersistentMemorySpace``.

   * - nnz_lno_persistent_work_host_view_t
     - ``host_mirror_type`` type of ``nnz_lno_persistent_work_host_view_t`` type.

   * - bool_persistent_view_t
     - View of rank 1 holding ``bool`` values and allocated in ``HandlePersistentMemorySpace``.

   * - bool_temp_view_t
     - View of rank 1 holding ``bool`` values and allocated in ``HandleTempMemorySpace``.

.. note::

   ``bool_persistent_view_t`` and ``bool_temp_view_t`` only appear in ``spgemm``

.. comment: we could deprecate these typedef and move them to spgemm. In general, I am not convinced that we really need to have all these type definitions here vs in the algorithms where they are actually used unless the handle itself holds data of the type defined.

Data Members
============

All data members of this class are private.

Member Functions
================

Generic functions
-----------------

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Name
     - Definition

   * - constructor
     - Both default and copy constructors are supported by ``KokkosKernelsHandle``.

   * - destructor
     - Destructs the handle and all the sub-handles.

   * - set_verbose
     - Set the verbosity level (true or false).

   * - get_verbose
     - Retrieve the verbosity level.

   * - set_team_work_size
     - Set the team work size.

   * - get_set_team_work_size
     - Retrieve the current team work size.

   * - get_team_work_size
     - Get a team work size based on previously set value and input parameters.

   * - get_handle_exec_space
     - Retrieve the execution space (a Kokkos Kernels execution space is an enum not a Kokkos execution space).

   * - set_dynamic_scheduling
     - Set the dynamic scheduling mode.

   * - is_dynamic_scheduling
     - Querry if the current scheduling is dynamic.

   * - set_shmem_size
     - Set the shmem size for kernels use.

   * - get_shmem_size
     - Get the current shmem size.

   * - set_suggested_vector_size
     - Set the suggested vector size.

   * - get_set_suggested_vector_size
     - Get the current suggested vector size.

   * - get_suggested_vector_size
     - Get a suggested vector size based on the one currently set and on input parameters.

   * - set_suggested_team_size
     - Set the suggested team size.

   * - get_set_suggested_team_size
     - Get the current suggested team size.

   * - get_suggested_team_size
     - Get a suggested team size based on the one currently set and on input parameters.

.. note::

   ``get_execution_space`` returns an enum that is defined in our Impl namespace and that has a name that can be confusing for the user.

.. comment: We might want to avoid this and rework that feature.

Sub-handle create/get/destroy
-----------------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Name
     - Definition

   * - :ref:`get_spadd_handle <handle_spadd_get>`
     - Get a pointer to the SpADD handle.

   * - :ref:`create_spadd_handle <handle_spadd_create>`
     - Construct a new SpADD handle after destroying the current one.

   * - :ref:`destroy_spadd_handle <handle_spadd_destroy>`
     - Destroy the currently owned SpADD handle.

   * - :ref:`get_spgemm_handle <handle_spgemm_get>`
     - Get a pointer to the SpGEMM handle.

   * - :ref:`create_spgemm_handle <handle_spgemm_create>`
     - Construct a new SpGEMM handle after destroying the current one.

   * - :ref:`destroy_spgemm_handle <handle_spgemm_destroy>`
     - Destroy the currently owned SpGEMM handle.

   * - :ref:`get_graph_coloring_handle <handle_coloring_get>`
     - Get a pointer to the graph coloring handle.

   * - :ref:`create_graph_coloring_handle <handle_coloring_create>`
     - Construct a new graph coloring handle from input parameters after destroying the current one.

   * - :ref:`destroy_graph_coloring_handle <handle_coloring_destroy>`
     - Destroy the currently owned graph coloring handle.

   * - :ref:`get_distance2_graph_coloring_handle <handle_d2coloring_get>`
     - Get a pointer to the distance two graph coloring handle.

   * - :ref:`create_distance2_graph_coloring_handle <handle_d2coloring_create>`
     - Construct a new distance two graph coloring handle from input parameters after destroying the current one.

   * - :ref:`destroy_distance2_graph_coloring_handle <handle_d2coloring_destroy>`
     - Destroy the currently owned distance two graph coloring handle.

   * - :ref:`get_gmres_handle <handle_gmres_get>`
     - Get a pointer to the GMRES handle.

   * - :ref:`create_gmres_handle <handle_gmres_create>`
     - Construct a new GMRES handle after destroying the current one.

   * - :ref:`destroy_gmres_handle <handle_gmres_destroy>`
     - Destroy the currently owned GMRES handle.

   * - :ref:`get_spiluk_handle <handle_spiluk_get>`
     - Get a pointer to the sparse ILUK handle.

   * - :ref:`create_spiluk_handle <handle_spiluk_create>`
     - Construct a new sparse ILUK handle after destroying the current one.

   * - :ref:`destroy_spiluk_handle <handle_spiluk_destroy>`
     - Destroy the currently owned sparse ILUK handle.

   * - :ref:`get_par_ilut_handle <handle_parilut_get>`
     - Get a pointer to the parallel ILUt handle.

   * - :ref:`create_par_ilut_handle <handle_parilut_create>`
     - Construct a new parallel ILUt handle after destroying the current one.

   * - :ref:`destroy_par_ilut_handle <handle_parilut_destroy>`
     - Destroy the currently owned parallel ILUt handle.

   * - :ref:`get_sptrsv_handle <handle_sptrsv_get>`
     - Get a pointer to the sptrsv handle

   * - :ref:`create_sptrsv_handle <handle_sptrsv_create>`
     - Destroy the current sptrsv handle and create a new one with the provided input parameters.

   * - :ref:`destroy_sptrsv_handle <handle_sptrsv_destroy>`
     - Destroy the sptrsv handle.

   * - :ref:`get_gs_handle <handle_gs_get>`
     - Get a pointer to the Gauss-Seidel handle

   * - :ref:`create_gs_handle <handle_gs_create>`
     - Destroy the current Gauss-Seidel handle and create a new one with the provided input parameters.

   * - :ref:`destroy_gs_handle <handle_gs_destroy>`
     - Destroy the Gauss-Seidel handle.

Gauss-Seidel functions
----------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Name
     - Definition

   * - get_point_gs_handle
     - Attempts to get a point Gauss-Seidel handle.

   * - get_cluster_gs_handle
     - Attempts to get a cluster Gauss-Seidel handle.

   * - get_twostage_gs_handle
     - Attempts to get a two-stage Gauss-Seidel handle.

   * - set_gs_set_num_outer_sweeps
     - Two stage Gauss-Seidel number of outer sweeps.

   * - set_gs_set_num_inner_sweeps
     - Two stage Gauss-Seidel number of inner sweeps.

   * - set_gs_set_inner_damp_factor
     - Two stage Gauss-Seidel inner damping factor.

   * - set_gs_twostage
     - Set the two stage variant to be used with this handle.

   * - set_gs_twostage_compact_form
     - Set handle to use a compact or classic recurence form.

   * - get_gs_sptrsvL_handle
     - Get a pointer to the Gauss-Seidel inner sptrsvL handle.

   * - get_gs_sptrsvU_handle
     - Get a pointer to the Gauss-Seidel inner sptrsvU handle.

   * - create_gs_sptrsvL_handle
     - Destroy the current Gauss-Seidel inner sptrsvL handle and create a new one with the provided input parameters.

   * - create_gs_sptrsvU_handle
     - Destroy the current Gauss-Seidel inner sptrsvU handle and create a new one with the provided input parameters.

   * - destroy_gs_sptrsvL_handle
     - Destroy the Gauss-Seidel inner sptrsvL handle.

   * - destroy_gs_sptrsvU_handle
     - Destroy the Gauss-Seidel inner sptrsvU handle.


`Additional functions available only when supernodal sptrsv is enabled`

.. list-table::
   :widths: 30 70
   :header-rows: 0
   :align: left

   * - set_sptrsv_verbose
     - Set the verbosity level.

   * - set_sptrsv_perm
     - Set the permutation vector (as pointer) stored in ``SPTRSVHandle``

   * - set_sptrsv_supernodes
     - Set the number of supernodes,  host view of supercols (map from supernodes to columns), and etree (as pointer) stored in ``SPTRSVHandle``

   * - set_sptrsv_diag_supernode_sizes
     - Set the supernode sizes unblocked and blocked stored in ``SPTRSVHandle``

   * - set_sptrsv_unit_diagonal
     - Set the unit diagonal flag stored in ``SPTRSVHandle``

   * - set_sptrsv_merge_supernodes
     - Set the merge supernodes flag stored in ``SPTRSVHandle``

   * - set_sptrsv_invert_diagonal
     - Set invert diagonal flag stored in ``SPTRSVHandle``

   * - get_sptrsv_invert_diagonal
     - Get invert diagonal flag stored in ``SPTRSVHandle``

   * - set_sptrsv_invert_offdiagonal
     - Set invert offdiagonal flag stored in ``SPTRSVHandle``

   * - get_sptrsv_invert_offdiagonal
     - Get invert offdiagonal flag stored in ``SPTRSVHandle``

   * - set_sptrsv_etree
     - Set the etree (as pointer) stored in ``SPTRSVHandle``

   * - set_sptrsv_column_major
     - Set the flag if data is column major stored in ``SPTRSVHandle``

   * - is_sptrsv_lower_tri
     - Get boolean for is_lower_tri stored in ``SPTRSVHandle``

   * - is_sptrsv_column_major
     - Get boolean for is_column_major stored in ``SPTRSVHandle``

   * - set_sptrsv_trmm_on_device
     - Set the flag for trmm_on_device stored in ``SPTRSVHandle``



Example
=======
