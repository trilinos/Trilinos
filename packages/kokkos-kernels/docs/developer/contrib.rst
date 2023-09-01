Contributing
============

Comment Style
-------------
We follow doxygen style comments for both external (API) and internal members. See https://www.doxygen.nl/manual/docblocks.html for details.
Our documentation can be generated using the `-DKokkosKernels_ENABLE_DOCS:BOOL=ON` cmake flag; see `Building the Documentation`.

In general, we prefer that the prototype has the doxygen style comment rather than the definition. If there is no prototype, then the definition should have the doxygen style comment.

.. code-block::
    :caption: API Doxygen Style Example

        /// \brief Blocking wrapper for accessing a Kokkos View.
        /// \tparam ViewValueType The value type (Scalar or Vector) of each view element
        /// \tparam ViewType The view type
        /// \param v The view handle
        /// \param m The requested row index of v
        /// \param n The requested col index of v
        /// \return If m and n are within the extents of v, a valid element of v;
        ///         otherwise, the last element of v.
        ///
        template <class ViewValueType, class ViewType>
        KOKKOS_INLINE_FUNCTION ViewValueType
        access_view_bounds_check(ViewType v, int m, int n, const BoundsCheck::Yes &);

.. code-block::
    :caption: Type Doxygen Style Example

    /// \class CooMatrix
    ///
    /// \brief Coordinate format implementation of a sparse matrix.
    ///
    /// \tparam ScalarType The type of scalar entries in the sparse matrix.
    /// \tparam OrdinalType The type of index entries in the sparse matrix.
    /// \tparam Device The Kokkos Device type.
    /// "Coo" stands for "coordinate format".
    template <class ScalarType>
    class CooMatrix {
      public:
        //! Type of each value in the matrix
        using scalar_type = ScalarType;

      private:
        size_type m_num_rows, m_num_cols;

      public:
        //! The data in the matrix
        scalar_type data;

      /// \brief Default constructor; constructs an empty sparse matrix.
      KOKKOS_INLINE_FUNCTION
      CooMatrix() : m_num_rows(0), m_num_cols(0) {}

**NOTE:** To have vscode generate the "\\\\\\" style stubs:

1. install the C/C++ IntelliSense, debugging, and code browsing extension.

2. go to Settings, Extensions, C/C++, Doxygen Documentation Generator Settings, and ensure the setting for Doxdocgen is "\\\\\\".

3. place your cursor on the line above `template ...` and type "\\\\\\".

Including your documentation with directives
--------------------------------------------
Rather than have the documentation generation system default to generating documentation for the entire code base,
we opt-in to what we would like to include in the generated documentation. To opt-in, simply place the publicly facing
function signature or the class name in the appropriate ReStructuredText file. For example, to document a sparse
function and class open up kokkos-kernels/docs/developer/apidocs/sparse.rst:

.. code-block::
    :caption: Function signature example

    coo2crs
    -------
    .. doxygenfunction:: KokkosSparse::coo2crs(DimType, DimType, RowViewType, ColViewType, DataViewType)
    .. doxygenfunction:: KokkosSparse::coo2crs(KokkosSparse::CooMatrix<ScalarType, OrdinalType, DeviceType, MemoryTraitsType, SizeType> &cooMatrix)

Note that only the signature is required. One may specify the parameter names and any default values, but this is not required.

.. code-block::
    :caption: User defined type example

    coomatrix
    ---------
    .. doxygenclass::    KokkosSparse::CooMatrix
      :members:

For a full list of available directives, see https://breathe.readthedocs.io/en/latest/.

Library policies
----------------

System-specific functions
-------------------------
For portability, any system-specific function that is not in the C++ standard should not be invoked from kokkos-kernels.

Upcasting and downcasting
-------------------------
TODO

Blocking and non-blocking interfaces
------------------------------------
All the APIs are non-blocking unless:
1. A TPL is enabled
2. The result vector resides on the host and work is offloaded to a device

When a TPL is enabled, we follow the blocking semantics of the TPL interface.

If no TPLs are enabled, callers can avoid blocking calls by using any overload which accepts a result vector type as a template argument.