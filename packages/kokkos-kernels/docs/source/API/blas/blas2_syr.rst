KokkosBlas::syr
###############

Defined in header: :code:`KokkosBlas2_syr.hpp`

.. code:: c++

  template <class ExecutionSpace, class XViewType, class AViewType>
  void syr(const ExecutionSpace& space, const char trans[], const char uplo[],
           const typename AViewType::const_value_type& alpha, const XViewType& x,
	   const AViewType& A)

  template <class XViewType, class AViewType>
  void syr(const char trans[], const char uplo[],
           const typename AViewType::const_value_type& alpha, const XViewType& x,
           const AViewType& A)

Perform a symmetric rank-1 update of matrix ``A`` by vector ``x`` with scaling factor ``alpha``,

.. math::

   A = A + \alpha (x * x^T).

1. Iterate of the entries of ``A`` and update them with the ``alpha`` scaled product of the entries of ``x`` and ``x^T`` using resources from the ``space`` instance.
2. Iterate of the entries of ``A`` and update them with the ``alpha`` scaled product of the entries of ``x`` and ``x^T`` using resources from the default instance of ``typename AViewType::execution_space``.

The function will throw a runtime exception if ``A.extent(0) != X.extent(0) || A.extent(1) != X.extent(0)``.

Parameters
==========

:space: execution space instance.

:trans: "T" for transpose and "H" for conjugate transpose of ``x``, all characters after the first one are ignored.

:alpha: scaling parameter for the rank update.

:x: vector used to perform the rank update.

:A: matrix being updated.

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `XViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XViewType::memory_space>::accessible``

- `AViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename AViewType::memory_space>::accessible``
  - ``std::is_same_v<typename AViewType::value_type, typename AViewType::non_const_value_type>``

Example
=======

.. code:: cppkokkos

  #include <Kokkos_Core.hpp>
  #include <KokkosBlas2_syr.hpp>

  int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
      constexpr int M = 5;

      Kokkos::View<double**> A("A", M, M);
      Kokkos::View<double*> x("X", M);

      Kokkos::deep_copy(A, 1.0);
      Kokkos::deep_copy(x, 3.0);

      const double alpha = double(1.0);

      KokkosBlas::syr("T", "U", alpha, x, A);
    }
    Kokkos::finalize();
  }
