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