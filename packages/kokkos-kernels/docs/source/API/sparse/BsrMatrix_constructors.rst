KokkosSparse::BsrMatrix<>::BsrMatrix
####################################

.. code:: cppkokkos

  BsrMatrix();

  template <typename SType, typename OType, class DType, class MTType, typename IType>
  explicit BsrMatrix(const BsrMatrix<SType, OType, DType, MTType, IType>& B);

  BsrMatrix(const std::string& label, const staticcrsgraph_type& graph, const OrdinalType& blockDimIn);

  BsrMatrix(const std::string& label, OrdinalType nrows, OrdinalType ncols, size_type annz, ScalarType* vals,
            OrdinalType* rows, OrdinalType* cols, OrdinalType blockdim, bool pad = false);

  BsrMatrix(const std::string& label,
            const OrdinalType nrows, const OrdinalType ncols, const size_type annz,
            const values_type& vals, const row_map_type& rows, const index_type& cols, const OrdinalType blockDimIn);

  BsrMatrix(const std::string& label, const OrdinalType& ncols, const values_type& vals,
            const staticcrsgraph_type& graph, const OrdinalType& blockDimIn);

  template <typename SType, typename OType, class DType, class MTType, typename IType>
  BsrMatrix(const KokkosSparse::CrsMatrix<SType, OType, DType, MTType, IType>& crs_mtx, const OrdinalType blockDimIn);

Constructs a BsrMatrix from specified inputs.

1. Default constructors with empty graph and values.
2. Copy constructor, it perform shallow copies of the input matrix.
3. Construct with a graph that will be shared, allocate values.
4. Constructor from COO matrix stored on Host. This constructor is fairly slow and mostly intended to help with testing and reading data from file. The ``pad`` input parameters is currently not used. 
5. Construct the matrix from views for the ``row_map``, ``entries`` and ``values`` using their copy constructors (shallow copies).
6. Construct the matrix from ``graph``, ``values`` and ``ncols`` using their respective copy constructors (shallow copies).
7. Constructor using a CrsMatrix with an appropriate block structure. The reduced graph and values of the BsrMatrix are allocated and constructed on host before copying them to the appropriate memory space.

..
   .. warning::

      Another question regarding the constructors, why are we not templating on the objects but rather on the underlying types: Ordinal, Scalar, MemoryTraits...

      Finally, we do not do any static asserts in the constructors which seems wrong... should we check that device is a Kokkos device, values are Views, graph is a StaticCrsGraph, etc...

Parameters
==========

:label: A label that can be used when views are constructed.

:B: An block sparse matrix that holds the input data.

:crs_mtx: An CrsMatrix that holds the input data for the constructor.

:graph: A graph of the structure of the matrix used as input.

:vals: The values to be held by the matrix.

:rows: The row offsets of the matrix rows.

:cols: The column indices of values in each row.

:nrows: The number of rows in the matrix.

:ncols: The number of columns non zero entries in the matrix, this value is an upper bound.

:annz: The number of point values in this Bsr (should be divisible by blockDimIn ^ 2)

:blockDimIn: The dimension of the blocks to be stored.

Type Requirements
-----------------

- `OrdinalType` must be a signed integer.
