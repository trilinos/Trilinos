KokkosSparse::spiluk_numeric
############################

Defined in header ``KokkosSparse_spiluk.hpp``

.. code:: c++

  template <typename KernelHandle, typename ARowMapType, typename AEntriesType, typename AValuesType,
            typename LRowMapType, typename LEntriesType, typename LValuesType, typename URowMapType,
            typename UEntriesType, typename UValuesType>
  void spiluk_numeric(KernelHandle* handle, typename KernelHandle::const_nnz_lno_t fill_lev, ARowMapType& A_rowmap,
                      AEntriesType& A_entries, AValuesType& A_values, LRowMapType& L_rowmap, LEntriesType& L_entries,
                      LValuesType& L_values, URowMapType& U_rowmap, UEntriesType& U_entries, UValuesType& U_values)

Performs the numeric phase of an incomplete LU factorization with level of fill k.

.. math::

   A \approx L*U

Parameters
==========

:handle: an instance of ``KokkosKernels::KokkosKernelsHandle`` from which an spiluk_handle will be used to extract control parameters.

:fill_lev: the level of fill to be used in the ILU factorization.

:A_rowmap, A_entries, A_values: rowmap, entries and values describing the input CrsMatrix.

:L_rowmap, L_entries, L_values: rowmap, entries and values of the lower triangular ``L`` matrix.

:U_rowmap, U_entries, U_values: rowmap, entries and values of the upper triangular ``U`` matrix.

Type Requirements
=================

Two main requirements are that the types of the ``rowmap``, ``entries`` and ``values`` should be compatible with the type used to build a CrsMatrix and they should be compatible with the types of the ``KernelsHandle``.

- ``A_rowmap``, ``A_entries``, ``A_values``, ``L_rowmap``, ``L_entries``, ``L_values``, ``U_rowmap``, ``U_entries`` and ``U_values`` should all be Kokkos View of rank 1 and share the same ``device_type``.

  - ``Kokkos::is_view_v<ARowMapType> && ARowMapType::rank() == 1``
  - ``Kokkos::is_view_v<AEntriesType> && AEntriesType::rank() == 1``
  - ``Kokkos::is_view_v<AValuesType> && AValuesType::rank() ==1``
  - ``Kokkos::is_view_v<LRowMapType> && LRowMapType::rank() == 1``
  - ``Kokkos::is_view_v<LEntriesType> && LEntriesType::rank() == 1``
  - ``Kokkos::is_view_v<LValuesType> && LValuesType::rank() == 1``
  - ``Kokkos::is_view_v<URowMapType> && URowMapType::rank() == 1``
  - ``Kokkos::is_view_v<UEntriesType> && UEntriesType::rank() == 1``
  - ``Kokkos::is_view_v<UValuesType> && UValuesType::rank() == 1``
  - ``std::is_same_v<typename LRowMapType::device_type, typename ARowMapType::device_type>``
  - ``std::is_same_v<typename LRowMapType::device_type, typename URowMapType::device_type>``
  - ``std::is_same_v<typename LEntriesType::device_type, typename AEntriesType::device_type>``
  - ``std::is_same_v<typename LEntriesType::device_type, typename UEntriesType::device_type>``
  - ``std::is_same_v<typename LValuesType::device_type, typename AValuesType::device_type>``
  - ``std::is_same_v<typename LValuesType::device_type, typename UValuesType::device_type>``
  - ``std::is_same_v<typename LRowMapType::device_type, typename LEntriesType::device_type>``
  - ``std::is_same_v<typename LRowMapType::device_type, typename LValuesType::device_type>``

- ``A_rowmap``, ``A_entries``, ``L_rowmap``, ``L_entries``, ``U_rowmap`` and ``U_entries`` should have value types and execution spaces compatible with the ``KernelsHandle``

  - ``std::is_same_v<typename ARowMapType::non_const_value_type, typename KernelHandle::size_type>``
  - ``std::is_same_v<typename AEntriesType::non_const_value_type, typename KernelHandle::ordinal_type>``
  - ``std::is_same_v<typename AValuesType::value_type, typename KernelHandle::scalar_type>``
  - ``std::is_same_v<typename LRowMapType::non_const_value_type, typename KernelHandle::size_type>``
  - ``std::is_same_v<typename LEntriesType::non_const_value_type, typename KernelHandle::ordinal_type>``
  - ``std::is_same_v<typename LValuesType::value_type, typename KernelHandle::scalar_type>``
  - ``std::is_same_v<typename URowMapType::non_const_value_type, typename KernelHandle::size_type>``
  - ``std::is_same_v<typename UEntriesType::non_const_value_type, typename KernelHandle::ordinal_type>``
  - ``std::is_same_v<typename UValuesType::value_type, typename KernelHandle::scalar_type>``
  - ``std::is_same_v<std::is_same<typename LRowMapType::device_type::execution_space, typename KernelHandle::SPILUKHandleType::execution_space>``
  - ``std::is_same_v<std::is_same<typename LEntriesType::device_type::execution_space, typename KernelHandle::SPILUKHandleType::execution_space>``
  - ``std::is_same_v<std::is_same<typename LValuesType::device_type::execution_space, typename KernelHandle::SPILUKHandleType::execution_space>``

- ``L_entries``, ``L_values``, ``U_entries`` and ``U_values`` are storing output data that needs to be modifiable (non-const)

  - ``std::is_same_v<typename LEntriesType::value_type, typename LEntriesType::non_const_value_type>``
  - ``std::is_same_v<typename LValuesType::value_type, typename LValuesType::non_const_value_type>``
  - ``std::is_same_v<typename UEntriesType::value_type, typename UEntriesType::non_const_value_type>``
  - ``std::is_same_v<typename UValuesType::value_type, typename UValuesType::non_const_value_type>``

Example
=======

.. literalinclude:: ../../../../example/sparse/KokkosSparse_example_spiluk.cpp
  :language: c++
  :lines: 3-


