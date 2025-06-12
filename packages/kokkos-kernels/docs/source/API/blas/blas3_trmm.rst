KokkosBlas::trmm
################

Defined in header: :code:`KokkosBlas3_trmm.hpp`

.. code:: c++

  // Version 1: Takes execution_space as argument
  template <class execution_space, class AViewType, class BViewType>
  void trmm(const execution_space& space, const char side[], const char uplo[],
            const char trans[], const char diag[],
	    typename BViewType::const_value_type& alpha,
	    const AViewType& A, const BViewType& B);

  // Version 2: Infers execution_space from AViewType
  template <class AViewType, class BViewType>
  void trmm(const char side[], const char uplo[], const char trans[], const char diag[],
            typename BViewType::const_value_type& alpha, const AViewType& A,
	    const BViewType& B)

Computes the product on the left or right of a triangular matrix ``A`` with a dense matrix ``B``, storing the result back in ``B``

.. math::

  B = alpha*op(A)*B\quad // Left \\\\
  B = alpha*B*op(A)\quad // Right


Implementation
=================

1. Version 1: check input control parameters, compute the matrix product using the resources of ``space``
2. Version 2: check input control parameters, compute the matrix product using the resources of the default instance of ``typename AViewType::execution_space``

The functions will throw a runtime exception if ``A.extent(0) != A.extent(1) || (uplo == 'L' ? B_m : B_n) != A_n`` or if the input control parameters are not supported, see Parameters section.

.. note::

   Currently we require ``A`` to be square but that is not what BLAS does, it only requires that ``A`` has compatible dimensions to be multipled with ``B``, otherwise ``A`` can be rectangular. Only the entries required for the storage of the triangular part of interest will be used.

Parameters
==========

:space: execution space instance.

:side: the side of ``B`` which will be multiplied by ``A``, valid values are ``L, l`` for multiplication on the left or ``R, r`` for multiplication on the right.

:uplo: whether the triangular matrix stored is an upper or lower factor, valid values are ``U, u`` for upper triangular factor or ``L, l`` for lower triangular factor.

:trans: operation applied to ``A`` while multiplying, valid values are ``N, n`` for no operator applied, ``T, t`` for transpose operator applied, ``C, c`` for conjugate transpose operator applied.

:diag: indicate if the diagonal entries of the triangular matrix are all unit, valid values are ``U, u`` for unit diagonal (actual entries are ingnored) or ``N, n`` for non-unit diagonal.

:alpha: scaling parameter for the matrix multiplication.

:A: the triangular matrix that will be multiplied to the dense matrix. Note that ``A`` does not have to actually be triangular, only the triangular part of interest (based on the value of ``uplo``) will be used.

:B: the dense matrix storing the result of the multiplication.

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `AViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename AViewType::memory_space>::accessible``

- `BViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename BViewType::memory_space>::accessible``
  - ``std::is_same_v<typename AViewType::value_type, typename AViewType::non_const_value_type>``

Example
=======

.. code:: cppkokkos

  #include<Kokkos_Core.hpp>
  #include<Kokkos_Random.hpp>
  #include<KokkosBlas3_trmm.hpp>

  int main(int argc, char* argv[]) {
    Kokkos::initialize();
    {
      int M = atoi(argv[1]);
      int N = atoi(argv[2]);
      int K = N;

      using ViewType = Kokkos::View<double**>;
      using Scalar   = typename ViewType::value_type;

      ViewType A("A",K,K);
      ViewType B("B",M,N);

      Kokkos::Random_XorShift64_Pool<typename ViewType::device_type::execution_space> rand_pool(13718);
      Kokkos::fill_random(A,rand_pool,Scalar(10));
      Kokkos::fill_random(B,rand_pool,Scalar(10));

      const Scalar alpha = 1.0;
   
      KokkosBlas::trmm("R","L","T","N",alpha,A,B);
    }
    Kokkos::finalize();
  }
