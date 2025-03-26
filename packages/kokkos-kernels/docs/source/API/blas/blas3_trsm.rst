KokkosBlas::trsm
################

Defined in header: :code:`KokkosBlas3_trsm.hpp`

.. code:: c++

  // Version 1: Takes execution_space as argument 
  template <class execution_space, class AViewType, class BViewType>
  void trsm(const execution_space& space, const char side[], const char uplo[],
            const char trans[], const char diag[],
            typename BViewType::const_value_type& alpha,
	    const AViewType& A, const BViewType& B);

  // Version 2: Infers execution_space from AViewType
  template <class AViewType, class BViewType>
  void trsm(const char side[], const char uplo[], const char trans[], const char diag[],
            typename BViewType::const_value_type& alpha, const AViewType& A,
	    const BViewType& B);

Solve a triangular matrix system of the form

.. math::

   A*X=alpha*B\quad    // left \\\\
   X*A=alpha*B\quad    // right

where ``A`` is a triangular matrix and ``X`` and ``B`` are vectors or matrices. The solution ``X`` is stored back in ``B``

Implementation
=================
1. Version 1: check input control parameters, solve the triangular system of equations using the resources of ``space``
2. Version 2: check input control parameters, solve the triangular system of equations using the resources of the default instance of ``typename AViewType::execution_space``


Parameters
==========

:space: execution space instance.

:side: control parameter specifying on which side the solver is applied, supported values are ``L, l`` for left side and ``R, r`` for right side.

:uplo: control parameter specifying if the triangular matrix is upper or lower triangular, supported values are ``U, u`` for upper and ``L, l`` for lower.

:trans: control parameter specifying what operation on the entires of ``A`` should be performed. Supported values are ``N, n`` for nothing, ``T, t`` for transpose mode and ``C, c`` for conjugate transpose mode.

:diag: control parameter specifying if the diagonal entries of ``A`` are equal to one, in which case these entries will be ignored. Supported values are ``U, u`` for unit diagonal and ``N, n`` for non-unit diagonal.

:alpha: a scaling factor applied while solving the ``A`` triangular system.

:A: triangular system of size ``k`` equal to:

    - ``A.extent(0)`` when ``side`` is ``L`` or ``l``,
    - ``A.extent(1)`` when ``side`` is ``R`` or ``r``.

   With ``uplo`` set as ``U`` or ``u`` only leading ``k`` by ``k`` upper triangular entries are considered, with ``uplo`` set as ``L`` or ``l`` only the leading ``k`` by ``k`` lower triangular entries are considered. Furthermore if ``diag`` is set to ``U`` or ``u`` then the diagonal entries are also ignored. Note that ``A`` does not have to actually be triangular, only the triangular part of interest (based on the value of ``uplo``) will be used.

:B: a dense matrix representing the right hand side of the system on entry and overwritten with the solution of the system on output.


Type Requirements
-----------------

Example
=======

.. code:: cppkokkos

  #include <Kokkos_Core.hpp>
  #include <Kokkos_Random.hpp>
  #include <KokkosBlas3_trsm.hpp>

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
   
      KokkosBlas::trsm("R","L","T","N",alpha,A,B);
    }
    Kokkos::finalize();
  }
