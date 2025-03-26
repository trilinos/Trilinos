KokkosBlas::iamax
#################

Defined in header: :code:`KokkosBlas1_iamax.hpp`

.. code:: c++

  template <class execution_space, class XVector>
  typename XVector::size_type iamax(const execution_space& space, const XVector& x);

  template <class XVector>
  typename XVector::size_type iamax(const XVector& x);

  template <class execution_space, class RV, class XMV>
  void iamax(const execution_space& space, const RV& R, const XMV& X);

  template <class RV, class XMV>
  void iamax(const RV& R, const XMV& X);

Finds the index of the entry in each column of ``X`` with the largest magnitude.

1. iterate over the entries in ``X``, identify the lowest-index entry with the largest magnitude, fence the ``space`` instance
2. iterate over the entries in ``X``, identify the lowest-index entry with the largest magnitude, fence the default instance associated with ``typename XVector::execution_space``
3. iterate over the entries in ``X``, store the lowest-index entry with the largest magnitude in each column of ``X`` using the resources of ``space``
4. iterate over the entries in ``X``, store the lowest-index entry with the largest magnitude in each column of ``X`` using the resources of the default instance associated with ``typename XMV::execution_space``

The function will throw a runtime exception if ``X.extent(1) != R.extent(0)``

The result (return value or stored value) from this function is undefined if the vector has no entries.

Parameters
==========

:space: an execution space instance

:R: view holding the indices of entry with the largest magnitude in each column of X

:X: view for which the we identify the index of the entry with the largest magnitude per column

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `XVector` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XVector::memory_space>::accessible == true``

- `RV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 or 1 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename RV::memory_space>::accessible == true``
  - ``std::is_same_v<typename RV::value_type, typename RV::non_const_value_type> == true``

- `XMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible == true``
  - ``RV::rank == XMV::rank - 1``

Example
=======

.. code:: c++

  #include <Kokkos_Core.hpp>
  #include <Kokkos_Random.hpp>
  #include <KokkosBlas1_iamax.hpp>
  
  int main(int argc, char* argv[]) {
    Kokkos::initialize();
    {
      int N = 100;
      if (argc >= 2) {
        N = atoi(argv[1]);
      }
  
      using ViewType  = Kokkos::View<double*>;
      using Scalar    = typename ViewType::non_const_value_type;
      using AT        = Kokkos::ArithTraits<Scalar>;
      using mag_type  = typename AT::mag_type;
      using size_type = typename ViewType::size_type;
  
      ViewType x("X", N);
  
      typename ViewType::HostMirror h_x = Kokkos::create_mirror_view(x);
  
      Kokkos::Random_XorShift64_Pool<typename ViewType::device_type::execution_space> rand_pool(13718);
      Kokkos::fill_random(x, rand_pool, Scalar(10));
  
      Kokkos::deep_copy(h_x, x);
  
      size_type max_loc = KokkosBlas::iamax(x);
  
      mag_type expected_result   = Kokkos::ArithTraits<mag_type>::min();
      size_type expected_max_loc = 0;
      for (int i = 0; i < N; i++) {
        mag_type val = AT::abs(h_x(i));
        if (val > expected_result) {
          expected_result  = val;
          expected_max_loc = i + 1;
        }
      }
  
      printf("Iamax of X: %i, Expected: %i\n", max_loc, expected_max_loc);
    }
    Kokkos::finalize();
  }

output:

.. code::

   Iamax of X: 60, Expected: 60
