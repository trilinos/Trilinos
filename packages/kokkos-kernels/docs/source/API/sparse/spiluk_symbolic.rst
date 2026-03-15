KokkosSparse::spiluk_symbolic
#############################

Defined in header ``KokkosSparse_spiluk.hpp``

.. code:: c++

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

.. literalinclude:: ../../../../example/sparse/KokkosSparse_example_spiluk.cpp
  :language: c++
  :lines: 3-


