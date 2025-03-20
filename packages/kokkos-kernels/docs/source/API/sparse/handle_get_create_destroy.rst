KokkosKernels::Experimental::KokkosKernelsHandle::{get,create,destroy}
######################################################################

Defined in header ``KokkosKernels_Handle.hpp``

.. code:: cppkokkos

  KernelHandeType* get_kernel_handle();
  void create_kernel_handle(Args...);
  void destroy_kernel_handle();

A generic set of ``get``, ``create`` and ``destroy`` member functions from the ``KokkosKernelsHandle`` to manage a sub-handle associated with the ``Kernel`` algorithm.

SpADD
=====

.. _handle_spadd_get:

get
---

.. code:: cppkokkos

  SPADDHandleType *get_spadd_handle();

Returns a pointer to the spadd sub-handle owned by this ``KokkosKernelsHandle``.

.. _handle_spadd_create:

create
------

.. code:: cppkokkos

  void create_spadd_handle(bool input_sorted = false, bool input_merged = false);

Destroy any previous spadd sub-handle and create a new owned one specifying if the associated input matrices is sorted and/or merged.

Parameters
^^^^^^^^^^

:input_sorted: Specify if the input matrices are sorted, by default they are assumed un-sorted.

:input_merged: Specify if the input matrices are merged, by default they are assumed un-merged.

.. _handle_spadd_destroy:

destroy
-------

.. code:: cppkokkos

  void destroy_spadd_handle();

Destroys the spadd sub-handle owned by this ``KokkosKernelsHandle``.

SpGEMM
======

.. _handle_spgemm_get:

get
---

.. code:: cppkokkos

  SPGEMMHandleType *get_spgemm_handle();

Returns a pointer to the spgemm sub-handle owned by this ``KokkosKernelsHandle``.

.. _handle_spgemm_create:

create
------

.. code:: cppkokkos

  void create_spgemm_handle(KokkosSparse::SPGEMMAlgorithm spgemm_algo = KokkosSparse::SPGEMM_DEFAULT);

Destroy any previous spgemm sub-handle and create a new owned one specifying which algorithm to use.

Parameters
^^^^^^^^^^

:spgemm_algo: Specify the algorithm to be used to perform the spgemm operation.

.. _handle_spgemm_destroy:

destroy
-------

.. code:: cppkokkos

  void destroy_spgemm_handle();

Destroys the spgemm sub-handle owned by this ``KokkosKernelsHandle``.

Graph Coloring
==============

.. _handle_coloring_get:

get
---

.. code:: cppkokkos

  GraphColoringHandleType *get_graph_coloring_handle();

.. _handle_coloring_create:

create
------

.. code:: cppkokkos

  void create_graph_coloring_handle(KokkosGraph::ColoringAlgorithm coloring_type = KokkosGraph::COLORING_DEFAULT);

.. _handle_coloring_destroy:

destroy
-------

.. code:: cppkokkos

  void destroy_graph_coloring_handle()

Distance 2 Graph Coloring
=========================

.. _handle_d2coloring_get:

get
---

.. code:: cppkokkos

  GraphColorDistance2HandleType *get_distance2_graph_coloring_handle();

.. _handle_d2coloring_create:

create
------

.. code:: cppkokkos

  void create_distance2_graph_coloring_handle(KokkosGraph::GraphColoringAlgorithmDistance2 coloring_type = KokkosGraph::COLORING_D2_DEFAULT);

.. _handle_d2coloring_destroy:

destroy
-------

.. code:: cppkokkos

  void destroy_distance2_graph_coloring_handle();

GMRES
=====

.. _handle_gmres_get:

get
---

.. code:: cppkokkos

  GMRESHandleType *get_gmres_handle();

.. _handle_gmres_create:

create
------

.. code:: cppkokkos

  void create_gmres_handle(const size_type m = 50, const typename GMRESHandleType::float_t tol = 1e-8, const size_type max_restart = 50);

.. note::

   I am not sure that using ``GMRESHandleType::float_t`` in the general handle is great, might want to use one of the public type defined in the general handle and eventually cast that to the type the sub-handle expect in the implementation of the create function.

.. _handle_gmres_destroy:

destroy
-------

.. code:: cppkokkos

  void destroy_gmres_handle();

SpILUK
======

.. _handle_spiluk_get:

get
---

.. code:: cppkokkos

  SPILUKHandleType *get_spiluk_handle();

.. _handle_spiluk_create:

create
------

.. code:: cppkokkos

  void create_spiluk_handle(KokkosSparse::Experimental::SPILUKAlgorithm algm, size_type nrows, size_type nnzL,
                            size_type nnzU, size_type block_size = 0);

.. _handle_spiluk_destroy:

destroy
-------

.. code:: cppkokkos

  void destroy_spiluk_handle();

Parallel ILUt
=============

.. _handle_parilut_get:

get
---

.. code:: cppkokkos

  PAR_ILUTHandleType *get_par_ilut_handle();

.. _handle_parilut_create:

create
------

.. code:: cppkokkos

  void create_par_ilut_handle(const size_type max_iter                                            = 20,
                              const typename PAR_ILUTHandleType::float_t residual_norm_delta_stop = 1e-2,
                              const typename PAR_ILUTHandleType::float_t fill_in_limit            = 0.75,
                              const bool async_update = false, const bool verbose = false);

.. _handle_parilut_destroy:

destroy
-------

.. code:: cppkokkos

  void destroy_par_ilut_handle();

SpTRSV
======

.. _handle_sptrsv_get:

get
---

.. code:: cppkokkos

  SPTRSVHandleType *get_sptrsv_handle();

.. _handle_sptrsv_create:

create
------

.. code:: cppkokkos

  void create_sptrsv_handle(KokkosSparse::Experimental::SPTRSVAlgorithm algm, size_type nrows, bool lower_tri,
                            size_type block_size = 0);

.. _handle_sptrsv_destroy:

destroy
-------

.. code:: cppkokkos

  void destroy_sptrsv_handle();

Gauss-Seidel
============

.. _handle_gs_get:

get
---

.. code:: cppkokkos

  GaussSeidelHandleType *get_gs_handle();

.. _handle_gs_create:

create
------

.. code:: cppkokkos

  void create_gs_handle(const HandleExecSpace &handle_exec_space, int num_streams,
                        KokkosSparse::GSAlgorithm gs_algorithm            = KokkosSparse::GS_DEFAULT,
                        KokkosGraph::ColoringAlgorithm coloring_algorithm = KokkosGraph::COLORING_DEFAULT)
  void create_gs_handle(KokkosSparse::GSAlgorithm gs_algorithm            = KokkosSparse::GS_DEFAULT,
                        KokkosGraph::ColoringAlgorithm coloring_algorithm = KokkosGraph::COLORING_DEFAULT);
  void create_gs_handle(KokkosSparse::ClusteringAlgorithm clusterAlgo, nnz_lno_t hint_verts_per_cluster,
                        KokkosGraph::ColoringAlgorithm coloring_algorithm = KokkosGraph::COLORING_DEFAULT);

.. _handle_gs_destroy:

destroy
-------

.. code:: cppkokkos

  void destroy_gs_handle();

