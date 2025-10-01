KokkosSparse::CcsMatrix<>::CcsMatrix
####################################

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION CcsMatrix();

  CcsMatrix(const std::string& /* label */, const OrdinalType nrows, const OrdinalType ncols, const size_type annz,
            const values_type& vals, const col_map_type& colmap, const index_type& rows);

Constructs a CcsMatrix from specified inputs.

1. Default constructors with empty graph and values.
2. Constructor which shallow-copies the provided views.

..
   .. warning::

      A couple of constructors are marked as KOKKOS_INLINE_FUNCTION which means they are collable on device. Is this something that is really intended? If so should the corresponding destructor also be marked as KOKKOS_INLINE_FUNCTION so it can be called from a device?

      Another question regarding the constructors, why are we not templating on the objects but rather on the underlying types: Ordinal, Scalar, MemoryTraits...

      Finally, we do not do any static asserts in the constructors which seems wrong... should we check that device is a Kokkos device, values are Views, graph is a StaticCcsGraph, etc...

Parameters
==========

:label: A label that can be used when views are constructed. This label will appear in the label used for the underlying Kokkos Views.

:nrows: The number of rows in the matrix.

:ncols: The number of columns in the matrix.

:annz: The number of nonzeros in the matrix. This should match the lengths of `vals` and `rows`.

:vals: The values to be held by the matrix.

:colmap: The offsets of the matrix columns.

:rows: The row indices of nonzero values in each column.

Type Requirements
-----------------

- `OrdinalType` must be a signed integer.
