SPARSE -- KokkosKernels sparse interfaces
=========================================

crsmatrix
---------
.. doxygenclass::    KokkosSparse::CrsMatrix
    :members:

ccsmatrix
---------
.. doxygenclass::    KokkosSparse::CcsMatrix
    :members:

coomatrix
---------
.. doxygenclass::    KokkosSparse::CooMatrix
    :members:

crs2ccs
-------
.. doxygenfunction:: KokkosSparse::crs2ccs(OrdinalType nrows, OrdinalType ncols, SizeType nnz, ValViewType vals, RowMapViewType row_map, ColIdViewType col_ids)
.. doxygenfunction:: KokkosSparse::crs2ccs(KokkosSparse::CrsMatrix<ScalarType, OrdinalType, DeviceType, MemoryTraitsType, SizeType> &crsMatrix)

ccs2crs
-------
.. doxygenfunction:: KokkosSparse::ccs2crs(OrdinalType nrows, OrdinalType ncols, SizeType nnz, ValViewType vals, ColMapViewType col_map, RowIdViewType row_ids)
.. doxygenfunction:: KokkosSparse::ccs2crs(KokkosSparse::CcsMatrix<ScalarType, OrdinalType, DeviceType, MemoryTraitsType, SizeType> &ccsMatrix)

coo2crs
-------
.. doxygenfunction:: KokkosSparse::coo2crs(DimType m, DimType n, RowViewType row, ColViewType col, DataViewType data)
.. doxygenfunction:: KokkosSparse::coo2crs(KokkosSparse::CooMatrix<ScalarType, OrdinalType, DeviceType, MemoryTraitsType, SizeType> &cooMatrix)

crs2coo
-------
.. doxygenfunction:: KokkosSparse::crs2coo(OrdinalType, OrdinalType, SizeType, ValViewType, RowMapViewType, ColIdViewType)
.. doxygenfunction:: KokkosSparse::crs2coo(KokkosSparse::CrsMatrix<ScalarType, OrdinalType, DeviceType, MemoryTraitsType, SizeType> &crsMatrix)

spmv
----
.. doxygenfunction:: KokkosSparse::spmv(const ExecutionSpace& space, KokkosKernels::Experimental::Controls controls, const char mode[], const AlphaType &alpha, const AMatrix &A, const XVector &x, const BetaType &beta, const YVector &y)
.. doxygenfunction:: KokkosSparse::spmv(KokkosKernels::Experimental::Controls controls, const char mode[], const AlphaType &alpha, const AMatrix &A, const XVector &x, const BetaType &beta, const YVector &y)
.. doxygenfunction:: KokkosSparse::spmv(const ExecutionSpace& space, const char mode[], const AlphaType &alpha, const AMatrix &A, const XVector &x, const BetaType &beta, const YVector &y)
.. doxygenfunction:: KokkosSparse::spmv(const char mode[], const AlphaType &alpha, const AMatrix &A, const XVector &x, const BetaType &beta, const YVector &y)

trsv
----
.. doxygenfunction:: KokkosSparse::trsv

spgemm
------
.. doxygenfunction:: spgemm_symbolic(KernelHandle& kh, const AMatrix& A, const bool Amode, const BMatrix& B, const bool Bmode, CMatrix& C)
.. doxygenfunction:: spgemm_numeric(KernelHandle& kh, const AMatrix& A, const bool Amode, const BMatrix& B, const bool Bmode, CMatrix& C)
.. doxygenfunction:: spgemm(const AMatrix& A, const bool Amode, const BMatrix& B, const bool Bmode)

block_spgemm
------------
.. doxygenfunction:: block_spgemm_symbolic(KernelHandle& kh, const AMatrixType& A, const bool transposeA, const BMatrixType& B,const bool transposeB, CMatrixType& C)
.. doxygenfunction:: block_spgemm_numeric(KernelHandle& kh, const AMatrix& A, const bool Amode, const BMatrix& B, const bool Bmode, CMatrix& C)

gauss_seidel
------------
.. doxygenfunction:: create_gs_handle(KokkosSparse::GSAlgorithm gs_algorithm, KokkosGraph::ColoringAlgorithm coloring_algorithm)
.. doxygenfunction:: create_gs_handle(const HandleExecSpace&, int, KokkosSparse::GSAlgorithm gs_algorithm, KokkosGraph::ColoringAlgorithm coloring_algorithm)
.. doxygenfunction:: create_gs_handle(KokkosSparse::ClusteringAlgorithm, nnz_lno_t, KokkosGraph::ColoringAlgorithm)
.. doxygenfunction:: gauss_seidel_symbolic(const ExecutionSpace &space, KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map, lno_nnz_view_t_ entries, bool is_graph_symmetric)
.. doxygenfunction:: gauss_seidel_symbolic(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map, lno_nnz_view_t_ entries, bool is_graph_symmetric)
.. doxygenfunction:: gauss_seidel_numeric(const ExecutionSpace &space, KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, bool is_graph_symmetric)
.. doxygenfunction:: gauss_seidel_numeric(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, bool is_graph_symmetric)
.. doxygenfunction:: gauss_seidel_numeric(const ExecutionSpace &space, KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, scalar_nnz_view_t_ given_inverse_diagonal, bool is_graph_symmetric)
.. doxygenfunction:: gauss_seidel_numeric(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, scalar_nnz_view_t_ given_inverse_diagonal, bool is_graph_symmetric)
.. doxygenfunction:: symmetric_gauss_seidel_apply(const ExecutionSpace &space, KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec, bool init_zero_x_vector, bool update_y_vector, typename KernelHandle::nnz_scalar_t omega, int numIter)
.. doxygenfunction:: symmetric_gauss_seidel_apply(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec, bool init_zero_x_vector, bool update_y_vector, typename KernelHandle::nnz_scalar_t omega, int numIter)
.. doxygenfunction:: forward_sweep_gauss_seidel_apply(const ExecutionSpace &space, KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec, bool init_zero_x_vector, bool update_y_vector, typename KernelHandle::nnz_scalar_t omega, int numIter)
.. doxygenfunction:: forward_sweep_gauss_seidel_apply(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec, bool init_zero_x_vector, bool update_y_vector, typename KernelHandle::nnz_scalar_t omega, int numIter)
.. doxygenfunction:: backward_sweep_gauss_seidel_apply(const ExecutionSpace &space, KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec, bool init_zero_x_vector, bool update_y_vector, typename KernelHandle::nnz_scalar_t omega, int numIter)
.. doxygenfunction:: backward_sweep_gauss_seidel_apply(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec, bool init_zero_x_vector, bool update_y_vector, typename KernelHandle::nnz_scalar_t omega, int numIter)

block_gauss_seidel
------------------
.. doxygenfunction:: block_gauss_seidel_symbolic(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols, typename KernelHandle::const_nnz_lno_t block_size, lno_row_view_t_ row_map, lno_nnz_view_t_ entries, bool is_graph_symmetric)
.. doxygenfunction:: block_gauss_seidel_numeric(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols,typename KernelHandle::const_nnz_lno_t block_size, lno_row_view_t_ row_map, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, bool is_graph_symmetric)
.. doxygenfunction:: symmetric_block_gauss_seidel_apply(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols, typename KernelHandle::const_nnz_lno_t block_size, lno_row_view_t_ row_map, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec, bool init_zero_x_vector, bool update_y_vector, typename KernelHandle::nnz_scalar_t omega, int numIter)
.. doxygenfunction:: forward_sweep_block_gauss_seidel_apply(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols, typename KernelHandle::const_nnz_lno_t block_size, lno_row_view_t_ row_map, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec, bool init_zero_x_vector, bool update_y_vector, typename KernelHandle::nnz_scalar_t omega, int numIter)
.. doxygenfunction:: backward_sweep_block_gauss_seidel_apply(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows, typename KernelHandle::const_nnz_lno_t num_cols, typename KernelHandle::const_nnz_lno_t block_size,  lno_row_view_t_ row_map, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec, bool init_zero_x_vector, bool update_y_vector, typename KernelHandle::nnz_scalar_t omega, int numIter)  

par_ilut
--------
.. doxygenfunction:: par_ilut_symbolic(KernelHandle* handle, ARowMapType& A_rowmap, AEntriesType& A_entries, LRowMapType& L_rowmap, URowMapType& U_rowmap)
.. doxygenfunction:: par_ilut_numeric(KernelHandle* handle, ARowMapType& A_rowmap, AEntriesType& A_entries, AValuesType& A_values, LRowMapType& L_rowmap, LEntriesType& L_entries, LValuesType& L_values, URowMapType& U_rowmap, UEntriesType& U_entries, UValuesType& U_values)
.. doxygenclass::    KokkosSparse::Experimental::PAR_ILUTHandle
    :members:

gmres
-----
.. doxygenfunction:: gmres(KernelHandle* handle, AMatrix& A, BType& B, XType& X, Preconditioner<AMatrix>* precond)

sptrsv
------
.. doxygenfunction:: sptrsv_symbolic(const ExecutionSpace &space, KernelHandle *handle, lno_row_view_t_ rowmap, lno_nnz_view_t_ entries)
.. doxygenfunction:: sptrsv_symbolic(KernelHandle *handle, lno_row_view_t_ rowmap, lno_nnz_view_t_ entries)
.. doxygenfunction:: sptrsv_symbolic(ExecutionSpace &space, KernelHandle *handle, lno_row_view_t_ rowmap, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values)
.. doxygenfunction:: sptrsv_symbolic(KernelHandle *handle, lno_row_view_t_ rowmap, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values)
.. doxygenfunction:: sptrsv_solve(ExecutionSpace &space, KernelHandle *handle, lno_row_view_t_ rowmap, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, BType b, XType x)
.. doxygenfunction:: sptrsv_solve(KernelHandle *handle, lno_row_view_t_ rowmap, lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, BType b, XType x)
.. doxygenfunction:: sptrsv_solve(ExecutionSpace &space, KernelHandle *handle, XType x, XType b)
.. doxygenfunction:: sptrsv_solve(KernelHandle *handle, XType x, XType b)
.. doxygenfunction:: sptrsv_solve(ExecutionSpace &space, KernelHandle *handleL, KernelHandle *handleU, XType x, XType b)
.. doxygenfunction:: sptrsv_solve(KernelHandle *handleL, KernelHandle *handleU, XType x, XType b)
