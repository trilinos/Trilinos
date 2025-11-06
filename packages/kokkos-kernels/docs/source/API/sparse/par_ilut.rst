KokkosSparse::par_ilut_symbolic
###############################

Parallel threshold incomplete LU factorization ILU(t)

This file provides KokkosSparse::par_ilut.  This function performs a
local (no MPI) sparse ILU(t) on matrices stored in
compressed row sparse ("Crs") format. It is expected that symbolic
is called before numeric. The handle offers an async_update
flag that controls whether asynchronous updates are allowed while computing
L U factors. This is useful for testing as it allows for repeatable
(deterministic) results but may cause the algorithm to take longer (more
iterations) to converge. The par_ilut algorithm will repeat (iterate) until
max_iters is hit or the improvement in the residual from iter to iter drops
below a certain threshold.

This algorithm is described in the paper:
PARILUT - A New Parallel Threshold ILU Factorization - Anzt, Chow, Dongarra

==========

Defined in header: :code:`KokkosSparse_par_ilut.hpp`

.. code:: cppkokkos

  template <typename KernelHandle, typename ARowMapType, typename AEntriesType, typename LRowMapType,
            typename URowMapType>
  void par_ilut_symbolic(KernelHandle* handle, ARowMapType& A_rowmap, AEntriesType& A_entries, LRowMapType& L_rowmap,
                         URowMapType& U_rowmap);

par_ilut_symbolic performs the symbolic phase of the above algorithm.

The sparsity pattern of A will be analyzed and L_rowmap and U_rowmap will be
populated with the L (lower triangular) and U (upper triagular) non-zero
counts respectively. Having a separate symbolic phase allows for reuse when
dealing with multiple matrices with the same sparsity pattern. This routine
will set some values on handle for symbolic info (row count, nnz counts).

.. math::

   A \approx L*U

Parameters
==========

:handle: an instance of ``KokkosKernels::KokkosKernelsHandle`` from which an par_ilut_handle will be used to extract control parameters. It is expected that create_par_ilut_handle has been called on it

:A_rowmap, A_entries: rowmap and entries describing the graph of the input CrsMatrix (inputs)

:L_rowmap: rowmap of the lower triangular ``L`` matrix. (output)

:U_rowmap: rowmap of the upper triangular ``U`` matrix. (output)

Type Requirements
-----------------

- ``A_rowmap`` A rank 1 view of size type
- ``A_entries`` A rank 1 view of ordinal type
- ``L_rowmap`` A rank 1 view of (non const) size type
- ``U_rowmap`` A rank 1 view of (non const) size type

KokkosSparse::par_ilut_numeric
##############################

Defined in header: :code:`KokkosSparse_par_ilut.hpp`

.. code:: cppkokkos

  template <typename KernelHandle, typename ARowMapType, typename AEntriesType, typename AValuesType,
            typename LRowMapType, typename LEntriesType, typename LValuesType, typename URowMapType,
            typename UEntriesType, typename UValuesType>
  void par_ilut_numeric(KernelHandle* handle, ARowMapType& A_rowmap, AEntriesType& A_entries, AValuesType& A_values,
                        LRowMapType& L_rowmap, LEntriesType& L_entries, LValuesType& L_values, URowMapType& U_rowmap,
                        UEntriesType& U_entries, UValuesType& U_values)

Performs the numeric phase (for specific CSRs, not reusable) of the
par_ilut algorithm (described above). This is a non-blocking
function. It is expected that par_ilut_symbolic has already been
called for the provided KernelHandle.

.. math::

   A \approx L*U

Parameters
==========

:handle: an instance of ``KokkosKernels::KokkosKernelsHandle`` from which an par_ilut_handle will be used to extract control parameters. It is expected that create_par_ilut_handle has been called on it.

:A_rowmap, A_entries, A_values: The input CrsMatrix (inputs)

:L_rowmap, L_entries, L_values: The output L CrsMatrix (outputs)

:U_rowmap, U_entries, U_values: The output U CrsMatrix (outputs)

Type Requirements
-----------------

- ``A_rowmap`` A rank 1 view of size type
- ``A_entries`` A rank 1 view of ordinal type
- ``A_values`` A rank 1 view of scalar type
- ``L_rowmap`` A rank 1 view of size type
- ``L_entries`` A rank 1 view of ordinal type
- ``L_values`` A rank 1 view of scalar type
- ``U_rowmap`` A rank 1 view of size type
- ``U_entries`` A rank 1 view of ordinal type
- ``U_values`` A rank 1 view of scalar type


KokkosSparse::PAR_ILUTHandle
############################

Defined in header: :code:`KokkosSparse_par_ilut_handle.hpp`

.. code:: cppkokkos

  PAR_ILUTHandle(const size_type max_iter, const float_t residual_norm_delta_stop, const float_t fill_in_limit,
                 const bool async_update, const bool verbose);

The handle for the par_ilut algorithm, this should be created from a KernelHandle.

.. math::

   A \approx L*U

Parameters
==========

:max_iter: Hard cap on the number of par_ilut iterations
:residual_norm_delta_stop: When the change in residual from
                           iteration to iteration drops below
                           this, the algorithm will stop (even if
                           max_iters has not been hit). If this is set to
                           zero, computing residual step will be skipped which
                           can reduce overall memory use and speed up the individual
                           iterations (it will always do max_iter iterations though).
:fill_in_limit: The threshold for removing candidates
                from the intermediate L and U is set such
                that the resulting sparsity pattern has
                at most `fill_in_limit` times the number
                of non-zeros of the ILU(0)
                factorization. This selection is executed
                separately for both factors L and U. A higher fill limit
                (2 or 3) may be necessary for very sparse matrices to achieve a
                good preconditioner but this will increase the resources needed
                by par_ilut.
:async_update: Whether compute LU factors should do asychronous
               updates. When ON, the algorithm will usually converge
               faster but it makes the algorithm non-deterministic.
:verbose: Print information while executing par_ilut

Example
=======

.. code:: cppkokkos

  #include "Kokkos_Core.hpp"
  #include "KokkosSparse_par_ilut.hpp"

  using scalar_t  = double;
  using lno_t     = int;
  using size_type = int;

  int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
       using exe_space      = typename device::execution_space;
       using mem_space      = typename device::memory_space;
       using RowMapType     = Kokkos::View<size_type*, device>;
       using EntriesType    = Kokkos::View<lno_t*, device>;
       using ValuesType     = Kokkos::View<scalar_t*, device>;
       using sp_matrix_type = KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type>;
       using KernelHandle =
           KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, exe_space, mem_space, mem_space>;
       using float_t = typename KokkosKernels::ArithTraits<scalar_t>::mag_type;

       // Create a diagonally dominant sparse matrix to test:
       //  par_ilut settings max_iters, res_delta_stop, fill_in_limit, and
       //  async_update are all left as defaults
       constexpr auto n             = 5000;
       constexpr auto m             = 15;
       constexpr auto tol           = ParIlut::TolMeta<float_t>::value;
       constexpr auto numRows       = n;
       constexpr auto numCols       = n;
       constexpr auto diagDominance = 1;
       constexpr bool verbose       = false;

       size_type nnz = 10 * numRows;
       auto A        = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<sp_matrix_type>(
           numRows, numCols, nnz, 0, lno_t(0.01 * numRows), diagDominance);

       KokkosSparse::sort_crs_matrix(A);

       // Make kernel handles
       KernelHandle kh;
       kh.create_gmres_handle(m, tol);
       auto gmres_handle = kh.get_gmres_handle();
       gmres_handle->set_verbose(verbose);
       using GMRESHandle    = typename std::remove_reference<decltype(*gmres_handle)>::type;
       using ViewVectorType = typename GMRESHandle::nnz_value_view_t;

       kh.create_par_ilut_handle();
       auto par_ilut_handle = kh.get_par_ilut_handle();
       par_ilut_handle->set_verbose(verbose);

       // Pull out views from CRS
       auto row_map = A.graph.row_map;
       auto entries = A.graph.entries;
       auto values  = A.values;

       // Allocate L and U CRS views as outputs
       RowMapType L_row_map("L_row_map", numRows + 1);
       RowMapType U_row_map("U_row_map", numRows + 1);

       // Initial L/U approximations for A
       par_ilut_symbolic(&kh, row_map, entries, L_row_map, U_row_map);

       const size_type nnzL = par_ilut_handle->get_nnzL();
       const size_type nnzU = par_ilut_handle->get_nnzU();

       EntriesType L_entries("L_entries", nnzL);
       ValuesType L_values("L_values", nnzL);
       EntriesType U_entries("U_entries", nnzU);
       ValuesType U_values("U_values", nnzU);

       par_ilut_numeric(&kh, row_map, entries, values, L_row_map, L_entries, L_values, U_row_map, U_entries, U_values);
    }
    Kokkos::finalize();
    return 0;
  }
