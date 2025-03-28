KokkosSparse::spiluk_symbolic
#############################

Defined in header ``KokkosSparse_spiluk.hpp``

.. code:: cppkokkos

  template <typename KernelHandle, typename ARowMapType, typename AEntriesType, typename LRowMapType,
            typename LEntriesType, typename URowMapType, typename UEntriesType>
  void spiluk_symbolic(KernelHandle* handle, typename KernelHandle::const_nnz_lno_t fill_lev, ARowMapType& A_rowmap,
                       AEntriesType& A_entries, LRowMapType& L_rowmap, LEntriesType& L_entries, URowMapType& U_rowmap,
                       UEntriesType& U_entries, int nstreams = 1);

Performs the symbolic phase of an incomplete LU factorization with level of fill k.

.. math::

   A \approx L*U

Parameters
==========

:handle: an instance of ``KokkosKernels::KokkosKernelsHandle`` from which an spiluk_handle will be used to extract control parameters.

:fill_lev: the level of fill to be used in the ILU factorization.

:A_rowmap, A_entries: rowmap and entries describing the graph of the input CrsMatrix.

:L_rowmap, L_entries: rowmap and entries of the lower triangular ``L`` matrix.

:U_rowmap, U_entries: rowmap and entries of the upper triangular ``U`` matrix.

:nstreams: the number of streams to be used while performing the numeric phase of the factorization (default: 1, if not set).

Type Requirements
=================

Two main requirements are that the types of the ``rowmap`` and ``entries`` should be compatible with the type used to build a StaticCrsGraph or a CrsMatrix and they should be compatible with the types of the ``KernelsHandle``.

- ``A_rowmap``, ``A_entries``, ``L_rowmap``, ``L_entries``, ``U_rowmap`` and ``U_entries`` should all be Kokkos View of rank 1 and share the same ``device_type``.

  - ``Kokkos::is_view_v<ARowMapType> == true && ARowMapType::rank() == 1``
  - ``Kokkos::is_view_v<AEntriesType> == true && AEntriesType::rank() == 1``
  - ``Kokkos::is_view_v<LRowMapType> == true && LRowMapType::rank() == 1``
  - ``Kokkos::is_view_v<LEntriesType> == true && LEntriesType::rank() == 1``
  - ``Kokkos::is_view_v<URowMapType> == true && URowMapType::rank() == 1``
  - ``Kokkos::is_view_v<UEntriesType> == true && UEntriesType::rank() == 1``
  - ``std::is_same_v<typename LRowMapType::device_type, typename ARowMapType::device_type>``
  - ``std::is_same_v<typename LRowMapType::device_type, typename URowMapType::device_type>``
  - ``std::is_same_v<typename LEntriesType::device_type, typename AEntriesType::device_type>``
  - ``std::is_same_v<typename LEntriesType::device_type, typename UEntriesType::device_type>``
  - ``std::is_same_v<typename LRowMapType::device_type, typename LEntriesType::device_type>``

- ``A_rowmap``, ``A_entries``, ``L_rowmap``, ``L_entries``, ``U_rowmap`` and ``U_entries`` should have value types and execution spaces compatible with the ``KernelsHandle``

  - ``std::is_same_v<typename ARowMapType::non_const_value_type, typename KernelHandle::size_type>``
  - ``std::is_same_v<typename AEntriesType::non_const_value_type, typename KernelHandle::ordinal_type>``
  - ``std::is_same_v<typename LRowMapType::non_const_value_type, typename KernelHandle::size_type>``
  - ``std::is_same_v<typename LEntriesType::non_const_value_type, typename KernelHandle::ordinal_type>``
  - ``std::is_same_v<typename URowMapType::non_const_value_type, typename KernelHandle::size_type>``
  - ``std::is_same_v<typename UEntriesType::non_const_value_type, typename KernelHandle::ordinal_type>``
  - ``std::is_same_v<std::is_same<typename LRowMapType::device_type::execution_space, typename KernelHandle::SPILUKHandleType::execution_space>``
  - ``std::is_same_v<std::is_same<typename LEntriesType::device_type::execution_space, typename KernelHandle::SPILUKHandleType::execution_space>``

- ``L_rowmap``, ``L_entries``, ``U_rowmap`` and ``U_entries`` are storing output data that needs to be modifiable (non-const)

  - ``std::is_same_v<typename LRowMapType::value_type, typename LRowMapType::non_const_value_type>``
  - ``std::is_same_v<typename LEntriesType::value_type, typename LEntriesType::non_const_value_type>``
  - ``std::is_same_v<typename URowMapType::value_type, typename URowMapType::non_const_value_type>``
  - ``std::is_same_v<typename UEntriesType::value_type, typename UEntriesType::non_const_value_type>``

Example
=======

.. code:: cppkokkos

  #include <Kokkos_Core.hpp>
  #include <KokkosSparse_CrsMatrix.hpp>
  #include <KokkosSparse_spiluk.hpp>
  #include <KokkosKernels_IOUtils.hpp>

  int main(int argc, char* argv[]) {
    Kokkos::initialize();
    {

      using scalar_t  = double;
      using lno_t     = int;
      using size_type = int;
      using crsMat_t  = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Kokkos::DefaultExecutionSpace, void, size_type>;

      using graph_t         = typename crsmat_t::StaticCrsGraphType;
      using lno_view_t      = typename graph_t::row_map_type::non_const_type;
      using lno_nnz_view_t  = typename graph_t::entries_type::non_const_type;
      using scalar_view_t   = typename crsmat_t::values_type::non_const_type;

      using ViewVectorType  = Kokkos::View<scalar_t*>;
      using execution_space = typename ViewVectorType::device_type::execution_space;
      using memory_space    = typename ViewVectorType::device_type::memory_space;

      using KernelHandle    = KokkosKernels::Experimental::KokkosKernelsHandle <size_type, lno_t, scalar_t, execution_space, memory_space, memory_space>;

      // Read and fill matrix
      crsmat_t A        = KokkosKernels::Impl::read_kokkos_crst_matrix<crsmat_t>("mtx filename");
      graph_t  graph    = A.graph;
      const size_type N = graph.numRows();
      typename KernelHandle::const_nnz_lno_t fill_lev = lno_t(2) ;
      const size_type nnzA = A.graph.entries.extent(0);

      // Create KokkosKernelHandle with an spiluk algorithm, limited by configuration at compile-time and set via the handle
      // Some options: {SEQLVLSCHD_RP, SEQLVLSCHD_TP1}
      KernelHandle kh;

      //kh.create_spiluk_handle(KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_RP, N, EXPAND_FACT*nnzA*(fill_lev+1), EXPAND_FACT*nnzA*(fill_lev+1));
      kh.create_spiluk_handle(KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_TP1, N, EXPAND_FACT*nnzA*(fill_lev+1), EXPAND_FACT*nnzA*(fill_lev+1));

      auto spiluk_handle = kh.get_spiluk_handle();

      lno_view_t     L_row_map("L_row_map", N + 1);
      lno_nnz_view_t L_entries("L_entries", spiluk_handle->get_nnzL());
      scalar_view_t  L_values ("L_values",  spiluk_handle->get_nnzL());
      lno_view_t     U_row_map("U_row_map", N + 1);
      lno_nnz_view_t U_entries("U_entries", spiluk_handle->get_nnzU());
      scalar_view_t  U_values ("U_values",  spiluk_handle->get_nnzU());

      KokkosSparse::Experimental::spiluk_symbolic(&kh, fill_lev, A.graph.row_map, A.graph.entries, 
                                                  L_row_map, L_entries, U_row_map, U_entries);

      Kokkos::resize(L_entries, spiluk_handle->get_nnzL());
      Kokkos::resize(L_values,  spiluk_handle->get_nnzL());
      Kokkos::resize(U_entries, spiluk_handle->get_nnzU());
      Kokkos::resize(U_values,  spiluk_handle->get_nnzU());

      spiluk_handle->set_team_size(16);
	  
      KokkosSparse::Experimental::spiluk_numeric(&kh, fill_lev, 
                                                 A.graph.row_map, A.graph.entries, A.values, 
                                                 L_row_map, L_entries, L_values, U_row_map, U_entries, U_values );

      kh.destroy_spiluk_handle();
    }
    Kokkos::finalize();
  }


