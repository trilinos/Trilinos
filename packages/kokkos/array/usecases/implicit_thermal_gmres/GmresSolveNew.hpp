/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <limits>
#include <CRSMatVec.hpp>
#include <impl/Kokkos_Timer.hpp>


template<typename Scalar, class DeviceType, unsigned N = 1>
struct MVScale;

template<typename Scalar>
struct MVScale<Scalar, KOKKOS_MACRO_DEVICE, 1> {
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef Scalar                  value_type ;
  typedef device_type::size_type  size_type ;

  // The multivector to scale elementwise.
  Kokkos::MultiVectorView<Scalar,device_type> Y;
  // The (multiplicative) scaling factor.
  Kokkos::ValueView<Scalar,device_type>     S;

  // Constructor: takes the multivector to scale, and the scaling
  // factor.
  MVScale(
    const Kokkos::MultiVectorView<Scalar,device_type> & argY ,
    const Kokkos::ValueView<Scalar,device_type>       & argS)
    : Y (argY), S (argS) {}

  // Computational kernel for parallel_for.
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator() (size_type iwork) const {
    Y(iwork) *= (*S); 
  }
};


template<typename Scalar, class DeviceType, unsigned N = 1>
struct Dot;

template<typename Scalar>
struct Dot<Scalar, KOKKOS_MACRO_DEVICE, 1> {
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef Scalar                  value_type ;
  typedef device_type::size_type  size_type ;
  typedef Kokkos::MultiVectorView<value_type,device_type> vector_type ;

  vector_type X, Y;

  Dot( const vector_type & argX ,
       const vector_type & argY )
  : X( argX ) , Y( argY ) {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator() (size_type iwork , value_type & update) const {
    update += X(iwork) * Y(iwork); 
  }

  KOKKOS_MACRO_DEVICE_FUNCTION
  static void join (volatile value_type& update, const volatile value_type& source) {
    update += source; 
  }

  KOKKOS_MACRO_DEVICE_FUNCTION
  static void init (value_type& update) {
    update = 0; 
  }
};


//////////////////////////////////////////////////////////////////////

template<class Scalar, class Device, class NextKernel>
class Sqrt {
public:
  typedef Device device_type;
  typedef typename device_type::size_type size_type;
  typedef Kokkos::ValueView<Scalar, Device> value;

  Sqrt (value& out, NextKernel next) 
  : out_ (out), next_ (next) {}

  template< typename ResultScalarType >
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator() (ResultScalarType& input) const {
    ResultScalarType output = sqrt (input);
    *out_ = output;
    next_ (output);
  }

  value result() const {
    return out_;
  }

private:
  value out_;
  NextKernel next_;
};

template<class Scalar, class Device, class NextKernel>
class Inv {
public:
  typedef Device device_type;
  typedef typename device_type::size_type size_type;
  typedef Kokkos::ValueView<Scalar, Device> value;

  Inv (value& out, NextKernel next) 
  : out_ (out), next_ (next) {}

  template< typename ResultScalarType >
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator() (ResultScalarType& input) const {
    ResultScalarType output = ResultScalarType(1) / input;
    *out_ = output;
    next_ (output);
  }

  value result() const {
    return out_;
  }

private:
  value out_;
  NextKernel next_;
};


template<class Scalar, class Device, class NextKernel>
class Store {
public:
  typedef Device device_type;
  typedef typename device_type::size_type size_type;
  typedef Kokkos::ValueView<Scalar, Device> value;
  typedef ??? matrix;

  Store (matrix& A, size_type i, size_type j, NextKernel next) 
    : A_ (A), i_ (i), j_ (j), next_ (next) {}

  template<typename ResultScalarType>
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator() (ResultScalarType& input) const {
    ResultScalarType output = ResultScalarType(1) / input;
    A_(i_, j_) = output;
    next_ (output);
  }

  value result() const {
    return out_;
  }

private:
  matrix A_;
  size_type i_, j_;
  NextKernel next_;
};


// Don't really need this class anymore.
// It's just a trivial wrapper for ValueView.
template<class Scalar, class Device>
class Val {
public:
  typedef Scalar scalar_type;
  typedef Device device_type;
  typedef typename device_type::size_type size_type;
  typedef Kokkos::ValueView<Scalar, Device> value;

  Val (const value& x) :
    : x_ (x) {}

  const scalar_type&
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator* () const {
    return *x_;
  }

private:
  value x_;
};

// Delay accessing the entry of the given matrix until in the kernel.
// We need to do this because the matrix's operator() (i,j) returns a
// Scalar&, not a ValueView.
template<class Scalar, class Device>
class MatVal {
public:
  typedef Scalar scalar_type;
  typedef Device device_type;
  typedef typename device_type::size_type size_type;
  typedef Kokkos::ValueView<Scalar, Device> value;
  typedef ??? matrix;

  Val (matrix& A, size_type i, size_type j) :
    : A_ (A), i_ (i), j_ (j) {}

  const scalar_type&
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator* () const {
    return A_(i_, j_);
  }

private:
  matrix A_;
  size_type i_, j_;
  NextKernel next_;
};



template<class ValueExpression, class VectorExpression>
class Ymax {
public:
  typedef VectorExpression vector_type;
  typedef typename VectorExpression::scalar_type scalar_type;
  typedef VectorExpression::device_type device_type;
  typedef typename device_type::size_type size_type;

  Ymax (const ValueExpression& alpha, 
	const VectorExpression& x, 
	const VectorExpression& y) :
    alpha_ (alpha), x_ (x), y_ (x) {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator() (size_type i) const {
    // ValueExpression defines unary operator*(), which returns "const
    // scalar_type&".  VectorExpression defines operator() (size_type
    // i), which returns either "scalar_type&" or "const
    // scalar_type&".
    y_(i) -= (*alpha_) * x_(i);
  }
  
private:
  ValueExpression alpha_;
  VectorExpression x_;
  VectorExpression y_;
};

//////////////////////////////////////////////////////////////////////

template <class Scalar, class DeviceType>
struct GMRES_Solve;

template<class Scalar>
struct GMRES_Solve<Scalar, KOKKOS_MACRO_DEVICE>
{
  typedef Scalar scalar_type
  typedef KOKKOS_MACRO_DEVICE device_type;
  typedef device_type::size_type size_type;
  typedef Kokkos::MultiVectorView<scalar_type, device_type> vector_type;
  typedef Kokkos::MultiVectorView<size_type, device_type> index_vector_type;

  typedef Kokkos::ValueView<scalar_type, device_type> value_type;
  typedef Kokkos::ValueView<scalar_type, Kokkos::DeviceHost> host_value_type;
  typedef Dot<scalar_type, device_type, 1> dot_type;

  static void 
  run (vector_type & A_value,
       index_vector_type & A_row,
       index_vector_type & A_offsets,
       vector_type & b,
       vector_type & x,
       const int numIters)
  {
    using Kokkos::create_labeled_multivector;
    using Kokkos::create_value;
    using Kokkos::deep_copy;
    using Kokkos::parallel_for;
    using Kokkos::parallel_reduce;
    using std::cout;
    using std::endl;

    // Convenient abbreviations.
    typedef scalar_type S;
    typedef device_type D;

    // Number of rows in the sparse matrix A.  A_row is the "rowptr"
    // array in the CSR (compressed sparse row) data structure, so it
    // has one more entry than the number of rows in A.  (That last
    // entry should hold the number of stored entries in the sparse
    // matrix.)
    const size_t rows = A_row.length() - 1;

    // Current iteration count.  This does not include computing the
    // initial residual vector.
    int iteration = 0;
    
    // Stash useful constants on the device.
    // 
    // TODO (mfh 19 Dec 2011) Explore using constant memory for these
    // constants; this saves (limited) cache space.
    value_type one  = create_value<S, D> ();
    value_type zero = create_value<S, D> ();
    deep_copy (one,  S (1));
    deep_copy (zero, S (0));

    //
    // Compute initial residual r = b - A*x.
    //
    vector_type r = create_labeled_multivector<vector_type>("r", rows);
    { // Compute r := A*x
      CRSMatVec<Scalar,device_type> A_mult_x (A_value, A_row, A_offsets, x, r);
      A_mult_x.apply(); 
    }
    // Compute r := b - one*(A*x).
    parallel_for (rows, Ymax<ValueView<S, D>, MultiVectorView<S, D> > (one, b, r));

    //
    // Compute initial residual norm beta = ||r||_2.  Also compute its
    // inverse, to make scaling easier.  Sqrt's op() takes the dot
    // product result as input, sets its constructor's input to the
    // square root of the input, and then invokes Inv's op().  The
    // latter sets beta_d_inv to the reciprocal of its input (which is
    // the result of the previous op()).
    //
    value beta_d = create_value<S, D>();
    value beta_d_inv = create_value<S, D>();
    parallel_reduce (rows, Dot<S, D> (r, r), Sqrt<S, D, Inv<S, D> > (beta_d, Inv<S, D> (beta_d_inv)));

    // Bring a copy of beta back to host memory and check if zero.
    S beta_h;
    deep_copy (beta_h, beta_d);
    if (beta_h == 0.0) {
      cout << "Initial residual norm was zero" << endl;
      return;
    } 

    // Allocate space Q for the Krylov basis vectors.  Q has as many
    // rows as b and m+1 columns.
    vector_type Q = create_labeled_multivector<vector_type>("Q", rows, num_iters+1);

    // Allocate the m+1 by m upper Hessenberg matrix H, and fill it with zeros.
    vector_type H = 
      create_labeled_multivector<vector_type> ("H", num_iters+1, num_iters);

    // Q(:, 0) = r .* (1/beta) (elementwise).
    vector_type q0 (Q, 0);
    parallel_for (rows, MVScale<ValueView<S, D>, MultiVectorView<S, D> > (q0, beta_d_inv));

    Kokkos::Impl::Timer wall_clock;

    value H_jp1j_inv = create_value<S, D>();
    for (size_t j = 1; j <= num_iters; ++j, ++iteration) {
      //Q(:,j) = A * Q(:, j-1)
      MultiVector qj1 = MultiVector(Q,j-1);
      MultiVector qj = MultiVector(Q,j);
      CRSMatVec<S, D> A_mult_q (A_value, A_row, A_offsets, qj1, qj);
      A_mult_q.apply();

      // Modified Gram-Schmidt (MGS) projection step.
      for (size_t k = 0; k <= j-1; ++k) { 
	MultiVector qk = MultiVector(Q, k);
	// Compute dot(qk,qj) and store the result in H(k,j-1).
        parallel_reduce (rows, Dot<S, D> (qk, qj), Store<S, D> (H, k, j-1));
	// Apply the result stored in H(k,j-1) (from the previous
	// step) to compute qj := qj - H(k,j-1) * qk.
	//
	// Ymax's first input is a "ValueExpression" (a template
	// parameter) whose operator() returns a scalar value, in the
	// kernel (not on the host).
	parallel_for (rows, Ymax<MatVal<S, D>, MultiVectorView<S, D> > (MatVal<S, D> (H, k, j-1), qk, qj));
      }

      // MGS normalization step.
      //
      // Compute the square root of the dot product result, store in
      // H(j+1,j), and assign its reciprocal to H_jp1j_inv.  (Sqrt's
      // constructor takes the next kernel as its input, which is
      // Store.  Store's constructor takes the next kernel as its
      // input, which is Inv.  Inv's constructor takes the target
      // value view as its input.)
      parallel_reduce (rows, Dot<S, D> (qj, qj), Sqrt<S, D, Store<S, D, Inv<S, D> > > (Store<S, D> (H, j+1, j, Inv<S, D> (H_jp1j_inv))));

      // Normalize Q(:,j).
      parallel_for (rows, MVScale<ValueView<S, D>, MultiVectorView<S, D> > (H_jp1j_inv, qj));

#ifdef TEST_KOKKOS_SYNC
      // Simulate a synchronous implementation by first bringing back
      // the dot product result to host memory, then copying it to
      // device memory again.
      scalar_type syncTestValue = 0;
      deep_copy (syncTestValue, h_jp1j_inv);
      deep_copy (h_jp1j_inv, syncTestValue);
#endif // TEST_KOKKOS_SYNC
    }
    device_type::wait_functor_completion();
    return wall_clock.seconds();
  }
};


