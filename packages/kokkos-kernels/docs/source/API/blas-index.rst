API: BLAS
#########

.. toctree::
   :maxdepth: 2
   :hidden:

   blas/blas1_abs
   blas/blas1_axpy
   blas/blas1_dot
   blas/blas1_iamax
   blas/blas1_nrm1
   blas/blas1_nrm2
   blas/blas1_rotg
   blas/blas1_rot
   blas/blas1_rotmg
   blas/blas1_rotm
   blas/blas1_scal
   blas/blas1_swap

   blas/blas2_gemv
   blas/blas2_ger
   blas/blas2_syr
   blas/blas2_syr2

   blas/blas3_gemm
   blas/blas3_trmm
   blas/blas3_trsm

Kokkos Kernels strive to make available most BLAS capabilities to users and even more through extensions. Our native implementations using the Kokkos Core library capabilities are providing portability and performance which allows users to explore the performance of early hardware where no other implementation is available. Additional we support a large third party layer that wraps specialized and high performance implementation of the BLAS specification.

Kokkos objects for linear algebra
=================================

Our implementations are heavily leveraging the capabilities of Kokkos::View to simplify the interface of the standard BLAS specification. A Kokkos::View should be thought of as a tensor, depending on its rank it can represent a scalar value, a vector or a matrix. Each rank has an associated extent that specifies the number of values stored by the rank and a stride that records the leading dimension of the rank. Having these encoded in the View means that Kokkos Kernels does not require the user to specify dimension of vector or matrices neither information like leading dimension of matrices. Another information carried by the View is the notion of increment used to step by a given spacing in a vector. That increment is encoded in the layout of the View using a LayoutStride. Below we show the signature of the gemv function of standard BLAS and Kokkos Kernels to illustrate our explanations:

.. code-block:: c

  dgemv(character TRANS,
        integer M,
        integer N,
        double precision ALPHA,
        double precision, dimension(lda,*) A,
        integer LDA,
        double precision, dimension(*) X,
        integer INCX,
        double precision BETA,
        double precision, dimension(*) Y,
        integer INCY)

.. code-block:: c++

  template <class ExecutionSpace, class AViewType, class XViewType, class YViewType>
  void gemv(const ExecutionSpace& space,
            const char trans[],
            typename AViewType::const_value_type& alpha,
            const AViewType& A,
            const XViewType& x,
            typename YViewType::const_value_type& beta,
            const YViewType& y)

Now we provide an example of LayoutStride usage to compute the nrm1 of the odd entries and even entries in a vector, of course the example can be adapted for larger increments.

.. code-block:: c++

  #include "Kokkos_Core.hpp"
  #include "KokkosBlas1_nrm1.hpp"

  int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
      using vector = Kokkos::View<double*, Kokkos::DefaultExecutionSpace>;

      // Create a vector of length 20.
      // Fill it with the value 2.
      // Compute its norm1: 20*2 = 40.
      vector a("a", 20);
      Kokkos::deep_copy(a, 2);
      auto nrm1_a = KokkosBlas::nrm1(a);

      // Create a strided layout of length 10 and stride 2
      // Make a new view using the data from vector a and the new layout
      // Set the values in the even vector to 1
      // Compute the norm1 of the even vector: 10*1 = 10
      Kokkos::LayoutStride layout(10, 2);
      Kokkos::View<double*, Kokkos::LayoutStride, Kokkos::DefaultExecutionSpace> even_a(a.data(), layout);
      Kokkos::deep_copy(even_a, 1);
      auto nrm1_even_a = KokkosBlas::nrm1(even_a);

      // Repeat the operations to create an odd vector using a pointer offset and the layout
      // Set the values of the odd vector to 3
      // Compute the norm1 of the odd vector: 10*3 = 30
      Kokkos::View<double*, Kokkos::LayoutStride, Kokkos::DefaultExecutionSpace> odd_a(a.data()+1, layout);
      Kokkos::deep_copy(odd_a, 3);
      auto nrm1_odd_a = KokkosBlas::nrm1(odd_a);

      // Print all the norms out and the values in the original vector.
      // Note that to do this with device views you will need to create
      // mirror views and copy data from device to host with deep_copy.
      std::cout << "nrm1(a)=" << nrm1_a << ", nrm1(even_a)=" << nrm1_even_a << ", nrm1(odd_a)=" << nrm1_odd_a << std::endl;
      std::cout << "a={ ";
      for(int idx = 0; idx < a.extent_int(0); ++idx) {
        std::cout << a(idx) << " ";
      }
      std::cout << "}" << std::endl;
    }
    Kokkos::finalize();
  }


types and precision
===================

Kokkos Kernels relies on template parameters to specify the type of input and output variables, this allows us to have generic implementation that are not precision specific. This is reflected in the name of our functions that do not include the usual S, D, C and Z letters that identify the parameters type in reference BLAS. While this makes our naming scheme a simpler, the use of template parameters also allows us to naturally support multi-precision operations. For instance the ``gemv`` function above can be called with vectors of double and a matrix of float.

..
   .. note::
      This is a potentially tricky topic and reading what ``std::blas`` does would likely help.
      Three questions come to mind regarding types and precision:

	1. What precision should we use to accumulate values? Using ``std::common_type_t<...>`` can help get the highest precision common type. If needed, after accumulation is performed, the result can be cast to lower precision.
	2. What return type is appropriate especially for functions that return scalar values as that type is currently not a template parameter but rather an implementation detail usually coming from: ``InnerProductSpaceTraits``.
	3. Deciding when a complex type is required for outputs and static assert correctly around that constraint.

Execution spaces, non-blocking behavior and thread safety
=========================================================

All Kokkos Kernels BLAS calls can be made with an execution space instance as the first parameter. This allows users to attach a stream/queue to the space instance which will be used to launch the BLAS kernel in the appropriate stream/queue. For this strategy to be useful, the BLAS functions are implemented in a non-blocking fashion except when a return value is expected, for instance when computing a norm or a dot product, since we need to obtain the return value before the following kernel can start. Finally while we strive to have a thread safe BLAS implementation, when using an execution space instance to set a select a particular stream/queue to execute a kernel, selecting the queue is not a thread safe operation. More detail on that can be obtained from the vendors documentation: cuBLAS.

BLAS support
============

Below are tables summarizing the currently supported function calls and third party libraries in Kokkos Kernels.

BLAS 1
------

.. list-table::
   :widths: 12 26 10 10 10 10 10
   :header-rows: 1

   * - BLAS Call
     - API Call
     - Reference
     - BLAS
     - cuBLAS
     - rocBLAS
     - oneMKL
   * - ROTG
     - :doc:`rotg(a, b, c, s) <blas/blas1_rotg>`
     - X
     - X
     - X
     - X
     - X
   * - ROTMG
     - rotmg(d1, d2, x1, y1, param)
     - X
     - X
     - X
     - X
     - --
   * - ROT
     - :doc:`rot(X, Y, c, s) <blas/blas1_rot>`
     - X
     - X
     - X
     - X
     - --
   * - ROTM
     - :doc:`rotm(X, Y, param) <blas/blas1_rotm>`
     - X
     - X
     - X
     - --
     - --
   * - SWAP
     - :doc:`swap(X, Y) <blas/blas1_swap>`
     - X
     - X
     - X
     - X
     - --
   * - SCAL
     - :doc:`scal(y,a,x) <blas/blas1_scal>`
     - X
     - X
     - X
     - X
     - --
   * - SCAL
     - :doc:`scal(y,a,x) <blas/blas1_scal>`
     - X
     - X
     - X
     - X
     - --
   * - COPY
     - `deep_copy(y,x) <https://kokkos.org/kokkos-core-wiki/API/core/view/deep_copy.html>`_
     - X
     - --
     - --
     - --
     - --
   * - AXPY
     - :doc:`axpy(a,x,y) <blas/blas1_axpy>`
     - X
     - X
     - X
     - --
     - --
   * - DOT*
     - :doc:`dot(x,y) <blas/blas1_dot>`
     - X
     - X
     - X
     - X
     - X
   * - DOTU
     - --
     - --
     - --
     - --
     - --
     - --
   * - DOTC*
     - :doc:`dot(x,y) <blas/blas1_dot>`
     - X
     - X
     - X
     - X
     - X
   * - NRM2
     - :doc:`nrm2(x) <blas/blas1_nrm2>`
     - X
     - X
     - X
     - X
     - X
   * - ASUM
     - :doc:`nrm1(x) <blas/blas1_nrm1>`
     - X
     - X
     - X
     - X
     - X
   * - IAMAX
     - :doc:`iamax(x) <blas/blas1_iamax>`
     - X
     - X
     - X
     - X
     - --

BLAS 2
------

.. list-table::
   :widths: 12 26 10 10 10 10 10
   :header-rows: 1

   * - BLAS Call
     - API Call
     - Reference
     - BLAS
     - cuBLAS
     - rocBLAS
     - oneMKL
   * - GEMV
     - :doc:`gemv(trans,a,A,x,b,y) <blas/blas2_gemv>`
     - X
     - X
     - X
     - X
     - X
   * - GBMV
     - --
     - --
     - --
     - --
     - --
     - --
   * - SYMV
     - --
     - --
     - --
     - --
     - --
     - --
   * - SBMV
     - --
     - --
     - --
     - --
     - --
     - --
   * - SPMV
     - --
     - --
     - --
     - --
     - --
     - --
   * - TRMV
     - :doc:`trmm(side,uplo,trans, diag,alpha,A,B) <blas/blas3_trmm>`
     - X
     - X
     - X
     - --
     - --
   * - TBMV
     - --
     - --
     - --
     - --
     - --
     - --
   * - TPMV
     - --
     - --
     - --
     - --
     - --
     - --
   * - TRSV
     - :doc:`trsm(side,uplo,trans,diag,alpha,A,B) <blas/blas3_trsm>`
     - X
     - X
     - X
     - --
     - --
   * - TBSV
     - --
     - --
     - --
     - --
     - --
     - --
   * - TPSV
     - --
     - --
     - --
     - --
     - --
     - --
   * - GER
     - :doc:`ger(trans,a,x,y,A) <blas/blas2_ger>`
     - X
     - X
     - X
     - X
     - --
   * - GERU
     - :doc:`ger(trans,a,x,y,A) <blas/blas2_ger>`
     - X
     - X
     - X
     - X
     - --
   * - GERC
     - :doc:`ger(trans,a,x,y,A) <blas/blas2_ger>`
     - X
     - X
     - X
     - X
     - --
   * - SYR
     - :doc:`syr(trans,uplo,a,x,A) <blas/blas2_syr>`
     - X
     - X
     - X
     - X
     - --
   * - SPR
     - --
     - --
     - --
     - --
     - --
     - --
   * - SYR2
     - :doc:`syr2(trans,uplo,a,x,y,A) <blas/blas2_syr2>`
     - X
     - X
     - X
     - X
     - --


BLAS 3
------

.. list-table::
   :widths: 12 26 10 10 10 10 10
   :header-rows: 1

   * - BLAS Call
     - API Call
     - Reference
     - BLAS
     - cuBLAS
     - rocBLAS
     - oneMKL
   * - GEMM
     - :doc:`gemm(transA,transB,a,A,B,b,C) <blas/blas3_gemm>`
     - X
     - X
     - X
     - X
     - --
   * - SYMM
     - --
     - --
     - --
     - --
     - --
     - --
   * - SYRK
     - --
     - --
     - --
     - --
     - --
     - --
   * - SYRK2
     - --
     - --
     - --
     - --
     - --
     - --
   * - HEMM
     - --
     - --
     - --
     - --
     - --
     - --
   * - HERK
     - --
     - --
     - --
     - --
     - --
     - --
   * - HERK2
     - --
     - --
     - --
     - --
     - --
     - --
   * - TRMM
     - :doc:`trmm(side, uplo, trans, diag, a, A, B) <blas/blas3_trmm>`
     - X
     - X
     - X
     - --
     - --
   * - TRSM
     - :doc:`trsm(side, uplo, trans, diag, a, A, B) <blas/blas3_trsm>`
     - X
     - X
     - X
     - --
     - --

