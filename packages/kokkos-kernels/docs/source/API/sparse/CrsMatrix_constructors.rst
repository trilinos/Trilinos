KokkosSparse::CrsMatrix<>::CrsMatrix
####################################

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION CrsMatrix();

  template <typename InScalar, typename InOrdinal, class InDevice, class InMemTraits, typename InSizeType>
  KOKKOS_INLINE_FUNCTION CrsMatrix(const CrsMatrix<InScalar, InOrdinal, InDevice, InMemTraits, InSizeType>& B);

  template <typename InScalar, typename InOrdinal, typename InDevice, typename InMemTraits, typename InSizeType>
  CrsMatrix(const std::string&, const CrsMatrix<InScalar, InOrdinal, InDevice, InMemTraits, InSizeType>& mat_);

  template <typename InOrdinal, typename InLayout, typename InDevice, typename InMemTraits, typename InSizeType>
  CrsMatrix(const std::string& label,
            const Kokkos::StaticCrsGraph<InOrdinal, InLayout, InDevice, InMemTraits, InSizeType>& graph_,
            const OrdinalType& ncols);

  template <typename InOrdinal, typename InLayout, typename InDevice, typename InMemTraits, typename InSizeType>
  CrsMatrix(const std::string&, const OrdinalType& ncols, const values_type& vals,
            const Kokkos::StaticCrsGraph<InOrdinal, InLayout, InDevice, InMemTraits, InSizeType>& graph_);

  CrsMatrix(const std::string& /*label*/, OrdinalType nrows, OrdinalType ncols, size_type annz, ScalarType* val,
            OrdinalType* rowmap, OrdinalType* cols);

  CrsMatrix(const std::string& /* label */, const OrdinalType nrows, const OrdinalType ncols, const size_type annz,
            const values_type& vals, const row_map_type& rowmap, const index_type& cols);

Constructs a CrsMatrix from specified inputs.

1. Default constructors with empty graph and values.
2. Copy constructor, it performs shallow copies of the underlying data into the constructed CrsMatrix.
3. Copy constructor, does a deep copy of the ``mat_`` into the constructed CrsMatrix. ``mat_`` and the constructed CrsMatrix can be in different memory spaces.
4. Constructor from existing ``graph_``. It makes a shallow copy of the graph and initializes the ``values`` view to zeros and sets its label, assign the number of columns to ``ncols``.
5. Construct the matrix from ``graph_``, ``values`` and ``ncols`` using their respective copy constructors (shallow copies).
6. Constructor from raw pointers on host, the pointers are wrapped into unmanaged views that are then deep copied into device views.
7. Constructor using input views and copy constructs the underlying graph and values view.

..
   .. warning::

      A couple of constructors are marked as KOKKOS_INLINE_FUNCTION which means they are collable on device. Is this something that is really intended? If so should the corresponding destructor also be marked as KOKKOS_INLINE_FUNCTION so it can be called from a device?

      Another question regarding the constructors, why are we not templating on the objects but rather on the underlying types: Ordinal, Scalar, MemoryTraits...

      Finally, we do not do any static asserts in the constructors which seems wrong... should we check that device is a Kokkos device, values are Views, graph is a StaticCrsGraph, etc...

Parameters
==========

:label: A label that can be used when views are constructed. This label will appear in the label used for the underlying Kokkos Views.

:mat\_: An input matrix that holds the input data.

:graph\_: A graph of the structure of the matrix used as input.

:ncols: The number of columns non zero entries in the matrix, this value is an upper bound.

:vals: The values to be held by the matrix.

:rowmap: The row offsets of the matrix rows.

:cols: The column indices of values in each row.

Type Requirements
-----------------

- `OrdinalType` must be a signed integer.
