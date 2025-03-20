KokkosSparse::spadd_symbolic
############################

Defined in header ``KokkosSparse_spadd.hpp``

.. code:: cppkokkos

  template <typename ExecSpace, typename KernelHandle, typename alno_row_view_t_,
            typename alno_nnz_view_t_, typename blno_row_view_t_,
            typename blno_nnz_view_t_, typename clno_row_view_t_>
  void spadd_symbolic(const ExecSpace &exec, KernelHandle *handle,
                      typename KernelHandle::const_nnz_lno_t m,  //same type as column indices
                      typename KernelHandle::const_nnz_lno_t n,
		      const alno_row_view_t_ a_rowmap,
                      const alno_nnz_view_t_ a_entries, const blno_row_view_t_ b_rowmap,
		      const blno_nnz_view_t_ b_entries, clno_row_view_t_ c_rowmap)

  template <typename KernelHandle, typename... Args>
  void spadd_symbolic(KernelHandle *handle, Args... args)

  template <typename ExecSpace, typename KernelHandle, typename AMatrix, typename BMatrix,
            typename CMatrix>
  void spadd_symbolic(const ExecSpace &exec, KernelHandle *handle, const AMatrix &A,
                      const BMatrix &B, CMatrix &C);

  template <typename KernelHandle, typename AMatrix, typename BMatrix, typename CMatrix>
  void spadd_symbolic(KernelHandle *handle, const AMatrix &A, const BMatrix &B, CMatrix &C);

Performs the symbolic phase of the addition of two sparse matrices.

.. math::

   C = \beta*C + \alpha*(A+B)


1. Computes the row map of ``C`` using the resources of ``exec`` and store the total number of non-zeros of ``C`` in the handle.
2. Same as 1. but use the resources of ``KernelHandle::HandleExecSpace{}``.
3. Compute the row map of ``C``, create and allocate views to hold column indices and values of ``C`` without initializing them using the resources of ``exec``. Construct matrix ``C`` from the row map and the newly created column indices and values views.
4. Same as 3. but use the resources of ``KernelHandle::HandleExecSpace{}``.

Parameters
==========

:exec: execution space instance.

:handle: spadd kernels handle obtained from an instance of ``KokkosKernels::KokkosKernelsHandle``.

:m, n: number of rows and column of the matrices ``A``, ``B`` and ``C``.

:a_rowmap, b_rowmap: row maps of the input matrices ``A`` and ``B``.

:a_entries, b_entries: column indices of the entries in each row of ``A`` and ``B``.

:c_rowmap: row map of the output matrix ``C``.

:A, B, C: three crs matrices.

Type Requirements
-----------------

Checks are performed in the :doc:`KokkosKernelsHandle <kokkoskernelshandle>` that was used to obtain the spadd handle.


Example
=======

.. code:: cppkokkos

  #include "Kokkos_Core.hpp"

  #include "KokkosKernels_default_types.hpp"
  #include "KokkosSparse_spadd.hpp"

  #include "KokkosKernels_Test_Structured_Matrix.hpp"

  using Scalar  = default_scalar;
  using Ordinal = default_lno_t;
  using Offset  = default_size_type;
  using Layout  = default_layout;

  int main(int argc, char* argv[]) {
    Kokkos::initialize();

    using device_type = typename Kokkos::Device<
        Kokkos::DefaultExecutionSpace,
        typename Kokkos::DefaultExecutionSpace::memory_space>;
    using execution_space = typename device_type::execution_space;
    using memory_space    = typename device_type::memory_space;
    using matrix_type =
        typename KokkosSparse::CrsMatrix<Scalar, Ordinal, device_type, void,
                                       Offset>;

    int return_value = 0;

    {
      // The mat_structure view is used to generate a matrix using
      // finite difference (FD) or finite element (FE) discretization
      // on a cartesian grid.
      // Each row corresponds to an axis (x, y and z)
      // In each row the first entry is the number of grid point in
      // that direction, the second and third entries are used to apply
      // BCs in that direction.
      Kokkos::View<Ordinal* [3], Kokkos::HostSpace> mat_structure(
          "Matrix Structure", 2);
      mat_structure(0, 0) = 10;  // Request 10 grid point in 'x' direction
      mat_structure(0, 1) = 1;   // Add BC to the left
      mat_structure(0, 2) = 1;   // Add BC to the right
      mat_structure(1, 0) = 10;  // Request 10 grid point in 'y' direction
      mat_structure(1, 1) = 1;   // Add BC to the bottom
      mat_structure(1, 2) = 1;   // Add BC to the top

      matrix_type A =
          Test::generate_structured_matrix2D<matrix_type>("FD", mat_structure);
      matrix_type B =
          Test::generate_structured_matrix2D<matrix_type>("FE", mat_structure);
      matrix_type C;

      // Create KokkosKernelHandle
      using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<
          Offset, Ordinal, Scalar, execution_space, memory_space, memory_space>;
      KernelHandle kh;
      kh.create_spadd_handle(false);

      const Scalar alpha = 2.5;
      const Scalar beta  = 1.2;


      KokkosSparse::spadd_symbolic(&kh, A, B, C);
      KokkosSparse::spadd_numeric(&kh, alpha, A, beta, B, C);
      kh.destroy_spadd_handle();

      std::cout << "spadd was performed correctly!" << std::endl;
    }

    Kokkos::finalize();

    return return_value;
  }
