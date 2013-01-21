// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#if !defined(Intrepid_MiniTensor_LinearAlgebra_t_cc)
#define Intrepid_MiniTensor_LinearAlgebra_t_cc

namespace Intrepid {

  //
  // R^N 2nd-order tensor inverse
  // Gauss-Jordan elimination. Warning: full pivoting
  // for small tensors. Use Teuchos LAPACK interface for
  // more efficient and robust techniques.
  // \param A nonsingular tensor
  // \return \f$ A^{-1} \f$
  //
  template<typename T>
  Tensor<T>
  inverse(Tensor<T> const & A)
  {
    Index const
    N = A.get_dimension();

    Tensor<T> S = A;
    Tensor<T> B = identity<T>(N);

    typedef std::set<Index> IndexSet;
    typedef std::set<Index>::const_iterator IndexIter;

    IndexSet intact_rows;
    IndexSet intact_cols;

    for (Index k = 0; k < N; ++k) {
      intact_rows.insert(k);
      intact_cols.insert(k);
    }

    // Gauss-Jordan elimination with full pivoting
    for (Index k = 0; k < N; ++k) {

      // Determine full pivot
      T pivot = 0.0;

      IndexIter pivot_row_iter = intact_rows.begin();
      IndexIter pivot_col_iter = intact_cols.begin();

      for (IndexIter rows_iter = intact_rows.begin();
          rows_iter != intact_rows.end(); ++rows_iter) {

        for (IndexIter cols_iter = intact_cols.begin();
            cols_iter != intact_cols.end(); ++cols_iter) {

          Index const row = *rows_iter;
          Index const col = *cols_iter;
          const T s = fabs(S(row, col));

          if (s > pivot) {

            pivot_row_iter = rows_iter;
            pivot_col_iter = cols_iter;

            pivot = s;

          }

        }

      }

      Index const pivot_row = *pivot_row_iter;
      Index const pivot_col = *pivot_col_iter;

      // Gauss-Jordan elimination
      const T t = S(pivot_row, pivot_col);

      if (t == 0.0) {
        std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
        std::cerr << std::endl;
        std::cerr << "Inverse of singular tensor.";
        std::cerr << std::endl;
        exit(1);
      }

      for (Index j = 0; j < N; ++j) {
        S(pivot_row, j) /= t;
        B(pivot_row, j) /= t;
      }

      for (Index i = 0; i < N; ++i) {
        if (i == pivot_row) continue;

        const T c = S(i, pivot_col);

        for (Index j = 0; j < N; ++j) {
          S(i, j) -= c * S(pivot_row, j);
          B(i, j) -= c * B(pivot_row, j);
        }
      }

      // Eliminate current row and col from intact rows and cols
      intact_rows.erase(pivot_row_iter);
      intact_cols.erase(pivot_col_iter);

    }

    Tensor<T>
    X = t_dot(S, B);

    return X;
  }

  //
  // R^N Subtensor
  // \param A tensor
  // \param i index
  // \param j index
  // \return Subtensor with i-row and j-col deleted.
  //
  template<typename T>
  Tensor<T>
  subtensor(Tensor<T> const & A, Index const i, Index const j)
  {
    Index const
    N = A.get_dimension();

    assert(i < N);
    assert(j < N);

    Tensor<T> B(N - 1);

    Index p = 0;
    Index q = 0;
    for (Index m = 0; m < N; ++m) {
      if (m == i) continue;
      for (Index n = 0; n < N; ++n) {
        if (n == j) continue;
        B(p, q) = A(m, n);
        ++q;
      }
      ++p;
    }

    return B;
  }

  //
  // Exponential map
  //
  template<typename T>
  Tensor<T>
  exp(Tensor<T> const & A)
  {
    return exp_pade(A);
  }

  //
  // R^N exponential map by Taylor series, radius of convergence is infinity
  // \param A tensor
  // \return \f$ \exp A \f$
  //
  template<typename T>
  Tensor<T>
  exp_taylor(Tensor<T> const & A)
  {
    Index const
    max_iter = 128;

    const T
    tol = machine_epsilon<T>();

    Index const
    N = A.get_dimension();

    Tensor<T>
    term = identity<T>(N);

    // Relative error taken wrt to the first term, which is I and norm = 1
    T
    relative_error = 1.0;

    Tensor<T>
    B = term;

    Index
    k = 0;

    while (relative_error > tol && k < max_iter) {
      term = T(1.0 / (k + 1.0)) * term * A;
      B = B + term;
      relative_error = norm_1(term);
      ++k;
    }

    return B;
  }

  namespace {

    //
    // Scaling parameter theta for scaling and squaring exponential.
    //
    template<typename T>
    T
    scaling_squaring_theta(Index const order)
    {
      assert(order > 0 && order < 22);

      T const theta[] =
      {
          0.0e-0, 3.7e-8, 5.3e-4, 1.5e-2, 8.5e-2, 2.5e-1, 5.4e-1, 9.5e-1,
          1.5e-0, 2.1e-0, 2.8e-0, 3.6e-0, 4.5e-0, 5.4e-0, 6.3e-0, 7.3e-0,
          8.4e-0, 9,4e-0, 1.1e+1, 1.2e+1, 1.3e+1, 1.4e+1
      };

      return theta[order];
    }

    //
    // Polynomial coefficients for Padé approximants.
    //
    template<typename T>
    T
    polynomial_coefficient(Index const order, Index const index)
    {
      assert(index <= order);

      T
      c = 0.0;

      switch (order) {

      default:
        std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
        std::cerr << std::endl;
        std::cerr << "Wrong order in Padé polynomial coefficient: ";
        std::cerr << order << std::endl;
        exit(1);
        break;

      case 3:
        {
          T const
          b[] = {120, 60, 12, 1};

          c = b[index];
        }
        break;

      case 5:
        {
          T const
          b[] = {30240, 15120, 3360, 420, 30, 1};

          c = b[index];
        }
        break;

      case 7:
        {
          T const
          b[] = {17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1};

          c = b[index];
        }
        break;

      case 9:
        {
          T const
          b[] = {17643225600, 8821612800, 2075673600, 302702400, 30270240,
              2162160, 110880, 3960, 90, 1};

          c = b[index];
        }
        break;

      case 13:
        {
          T const
          b[] = {64764752532480000, 32382376266240000, 7771770303897600,
              1187353796428800, 129060195264000, 10559470521600,
              670442572800, 33522128640, 1323241920, 40840800,
              960960, 16380, 182, 1};

          c = b[index];
        }
        break;

      }

      return c;
    }

    //
    // Padé approximant polynomial odd and even terms.
    //
    template<typename T>
    std::pair<Tensor<T>, Tensor<T> >
    pade_polynomial_terms(Tensor<T> const & A, Index const order)
    {
      Index const
      N = A.get_dimension();

      Tensor<T>
      B = identity<T>(N);

      Tensor<T>
      U = polynomial_coefficient<Real>(order, 1) * B;

      Tensor<T>
      V = polynomial_coefficient<Real>(order, 0) * B;

      Tensor<T> const
      A2 = A * A;

      for (Index i = 3; i <= order; i += 2) {

        B = B * A2;

        Tensor<T> const
        O = polynomial_coefficient<Real>(order, i) * B;

        Tensor<T> const
        E = polynomial_coefficient<Real>(order, i - 1) * B;

        U += O;

        V += E;

      }

      U = A * U;

      return std::make_pair(U, V);
    }

    //
    // Compute an integer power of a tensor by binary manipulation.
    //
    template<typename T>
    Tensor<T>
    binary_powering(Tensor<T> const & A, Index const exponent)
    {
      if (exponent == 0) return A;

      Index const
      rightmost_bit = 1;

      Index const
      number_digits = std::numeric_limits<Index>::digits;

      Index const
      leftmost_bit = rightmost_bit << (number_digits - 1);

      Index
      t = 0;

      for (Index j = 0; j < number_digits; ++j) {

        if (((exponent << j) & leftmost_bit) != 0) {

          t = number_digits - j - 1;
          break;

        }

      }

      Tensor<T>
      P = A;

      Index
      i = 0;

      Index
      m = exponent;

      while ((m & rightmost_bit) == 0) {
        P = P * P;
        ++i;
        m = m >> 1;
      }

      Tensor<T>
      X = P;

      for (Index j = i + 1; j <= t; ++j) {
        P = P * P;

        if (((exponent >> j) & rightmost_bit) != 0) {
          X = X * P;
        }
      }

      return X;
    }

  } // anonymous namespace

  //
  // Exponential map by squaring and scaling and Padé approximants.
  // See algorithm 10.20 in Functions of Matrices, N.J. Higham, SIAM, 2008.
  // \param A tensor
  // \return \f$ \exp A \f$
  //
  template<typename T>
  Tensor<T>
  exp_pade(Tensor<T> const & A)
  {
    Index const
    N = A.get_dimension();

    Index const
    orders[] = {3, 5, 7, 9, 13};

    Index const
    number_orders = 5;

    Index const
    highest_order = orders[number_orders - 1];

    Tensor<T>
    B;

    Real const
    norm = Sacado::ScalarValue<T>::eval((norm_1(A)));

    for (Index i = 0; i < number_orders; ++i) {

      Index const
      order = orders[i];

      Real const
      theta = scaling_squaring_theta<Real>(order);

      if (order < highest_order && norm < theta) {

        Tensor<T>
        U;

        Tensor<T>
        V;

        boost::tie(U, V) = pade_polynomial_terms(A, order);

        B = inverse(V - U) * (U + V);

        break;

      } else if (order == highest_order) {

        Real const
        theta_highest = scaling_squaring_theta<Real>(order);

        Index const
        power_two = Index(std::ceil(log2(norm / theta_highest)));

        Real
        scale = 1.0;

        for (Index i = 0; i < power_two; ++i) {
          scale /= 2.0;
        }

        Tensor<T> const
        I = identity<T>(N);

        Tensor<T> const
        A1 = scale * A;

        Tensor<T> const
        A2 = A1 * A1;

        Tensor<T> const
        A4 = A2 * A2;

        Tensor<T> const
        A6 = A2 * A4;

        Real const b0  = polynomial_coefficient<Real>(order, 0);
        Real const b1  = polynomial_coefficient<Real>(order, 1);
        Real const b2  = polynomial_coefficient<Real>(order, 2);
        Real const b3  = polynomial_coefficient<Real>(order, 3);
        Real const b4  = polynomial_coefficient<Real>(order, 4);
        Real const b5  = polynomial_coefficient<Real>(order, 5);
        Real const b6  = polynomial_coefficient<Real>(order, 6);
        Real const b7  = polynomial_coefficient<Real>(order, 7);
        Real const b8  = polynomial_coefficient<Real>(order, 8);
        Real const b9  = polynomial_coefficient<Real>(order, 9);
        Real const b10 = polynomial_coefficient<Real>(order, 10);
        Real const b11 = polynomial_coefficient<Real>(order, 11);
        Real const b12 = polynomial_coefficient<Real>(order, 12);
        Real const b13 = polynomial_coefficient<Real>(order, 13);

        Tensor<T> const
        U = A1 * (
            (A6 * (b13 * A6 + b11 * A4 + b9 * A2) +
             b7 * A6 + b5 * A4 + b3 * A2 + b1 * I));

        Tensor<T> const
        V = A6 * (b12 * A6 + b10 * A4 + b8 * A2) +
          b6 * A6 + b4 * A4 + b2 * A2 + b0 * I;

        Tensor<T> const
        R = inverse(V - U) * (U + V);

        Index const
        exponent = (1U << power_two);

        B = binary_powering(R, exponent);

      }

    }

    return B;
  }

  //
  // R^N logarithmic map by Taylor series, converges for \f$ |A-I| < 1 \f$
  // \param A tensor
  // \return \f$ \log A \f$
  //
  template<typename T>
  Tensor<T>
  log(Tensor<T> const & A)
  {
    // Check whether skew-symmetric holds

    Index const
    max_iter = 128;

    const T
    tol = machine_epsilon<T>();

    const T
    norm_arg = norm_1(A);

    Index const
    N = A.get_dimension();

    const Tensor<T>
    Am1 = A - identity<T>(N);

    Tensor<T>
    term = Am1;

    T
    norm_term = norm_1(term);

    T
    relative_error = norm_term / norm_arg;

    Tensor<T>
    B = term;

    Index
    k = 1;

    while (relative_error > tol && k <= max_iter) {
      term = T(- (k / (k + 1.0))) * term * Am1;
      B = B + term;
      norm_term = norm_1(term);
      relative_error = norm_term / norm_arg;
      ++k;
    }

    return B;
  }

  //
  // R^N logarithmic map of a rotation. Not implemented yet.
  // \param R with \f$ R \in SO(N) \f$
  // \return \f$ r = \log R \f$ with \f$ r \in so(N) \f$
  //
  template<typename T>
  Tensor<T>
  log_rotation(Tensor<T> const & R)
  {
    Index const
    N = R.get_dimension();

    //firewalls, make sure R \in SO(N)
    assert(norm(dot_t(R,R) - eye<T>(N)) < 100.0 * machine_epsilon<T>());
    assert(abs(det(R) - 1.0) < 100.0 * machine_epsilon<T>());

    // acos requires input between -1 and +1
    T
    cosine = 0.5*(trace(R) - 1.0);

    if (cosine < -1.0) {
      cosine = -1.0;
    } else if(cosine > 1.0) {
      cosine = 1.0;
    }

    T
    theta = acos(cosine);

    Tensor<T>
    r(N);

    switch (N) {

    default:
      std::cerr << "Logarithm of SO(N) N != 2,3 not implemented." << std::endl;
      exit(1);
      break;

    case 3:
      if (theta == 0.0) {

        r = zero<T>(3);

      } else if (fabs(cosine + 1.0) < 10.0*machine_epsilon<T>())  {

        r = log_rotation_pi(R);

      } else {

        r = static_cast<T>(theta / (2.0 * sin(theta))) * (R - transpose(R));

      }
      break;

    case 2:
      r(0,0) = 0.0;
      r(0,1) = -theta;
      r(1,0) = theta;
      r(1,1) = 0.0;
      break;

    }

    return r;
  }

  // R^N Logarithmic map of a 180-degree rotation. Not implemented.
  // \param R with \f$ R \in SO(N) \f$
  // \return \f$ r = \log R \f$ with \f$ r \in so(N) \f$
  //
  template<typename T>
  Tensor<T>
  log_rotation_pi(Tensor<T> const & R)
  {
    Index const
    N = R.get_dimension();

    T
    cosine = 0.5*(trace(R) - 1.0);
    // set firewall to make sure the rotation is indeed 180 degrees
    assert(fabs(cosine + 1.0) < 10.0 * machine_epsilon<T>());

    Tensor<T>
    r(N);

    switch (N) {

    default:
      std::cerr << "Logarithm of SO(N) N != 2,3 not implemented." << std::endl;
      exit(1);
      break;

    case 3:
      {
        // obtain U from R = LU
        r = gaussian_elimination((R - identity<T>(3)));

        // backward substitution (for rotation exp(R) only)
        const T tol = 10.0 * machine_epsilon<T>();

        Vector<T> normal(3);

        if (fabs(r(2,2)) < tol){
          normal(2) = 1.0;
        } else {
          normal(2) = 0.0;
        }

        if (fabs(r(1,1)) < tol){
          normal(1) = 1.0;
        } else {
          normal(1) = -normal(2)*r(1,2)/r(1,1);
        }

        if (fabs(r(0,0)) < tol){
          normal(0) = 1.0;
        } else {
          normal(0) = -normal(1)*r(0,1) - normal(2)*r(0,2)/r(0,0);
        }

        normal = normal / norm(normal);

        r.clear();
        r(0,1) = -normal(2);
        r(0,2) =  normal(1);
        r(1,0) =  normal(2);
        r(1,2) = -normal(0);
        r(2,0) = -normal(1);
        r(2,1) =  normal(0);

        const T pi = acos(-1.0);

        r = pi * r;
      }
      break;

    case 2:
      {
        T theta = acos(-1.0);

        if (R(0,0) > 0.0) {
          theta = - theta;
        }

        r(0,0) = 0.0;
        r(0,1) = -theta;
        r(1,0) = theta;
        r(1,1) = 0.0;
      }
      break;

    }

    return r;
  }

  // Gaussian Elimination with partial pivot
  // \param matrix \f$ A \f$
  // \return \f$ U \f$ where \f$ A = LU \f$
  //
  template<typename T>
  Tensor<T>
  gaussian_elimination(Tensor<T> const & A)
  {
    Index const
    N = A.get_dimension();

    Tensor<T>
    U = A;

    const T
    tol = 10.0 * machine_epsilon<T>();

    Index i = 0;
    Index j = 0;
    Index i_max = 0;

    while ((i <  N) && (j < N)) {
      // find pivot in column j, starting in row i
      i_max = i;
      for (Index k = i + 1; k < N; ++k) {
        if (fabs(U(k,j) > fabs(U(i_max,j)))) {
          i_max = k;
        }
      }

      // Check if A(i_max,j) equal to or very close to 0
      if (fabs(U(i_max,j)) > tol){
        // swap rows i and i_max and divide each entry in row i
        // by U(i,j)
        for (Index k = 0; k < N; ++k) {
          std::swap(U(i,k), U(i_max,k));
        }

        for (Index k = 0; k < N; ++k) {
          U(i,k) = U(i,k) / U(i,j);
        }

        for (Index l = i + 1; l < N; ++l) {
          for (Index k = 0; k < N; ++k) {
            U(l,k) = U(l,k) - U(l,i) * U(i,k) / U(i,i);
          }
        }
        ++i;
      }
      ++j;
    }

    return U;
  }

  // Apply Givens-Jacobi rotation on the left in place.
  // \param c and s for a rotation G in form [c, s; -s, c]
  // \param A
  //
  template<typename T>
  void
  givens_left(T const & c, T const & s, Index i, Index k, Tensor<T> & A)
  {
    Index const
    N = A.get_dimension();

    for (Index j = 0; j < N; ++j) {
      T const t1 = A(i,j);
      T const t2 = A(k,j);
      A(i,j) = c * t1 - s * t2;
      A(k,j) = s * t1 + c * t2;
    }
    return;
  }

  // Apply Givens-Jacobi rotation on the right in place.
  // \param A
  // \param c and s for a rotation G in form [c, s; -s, c]
  //
  template<typename T>
  void
  givens_right(T const & c, T const & s, Index i, Index k, Tensor<T> & A)
  {
    Index const
    N = A.get_dimension();

    for (Index j = 0; j < N; ++j) {
      T const t1 = A(j,i);
      T const t2 = A(j,k);
      A(j,i) = c * t1 - s * t2;
      A(j,k) = s * t1 + c * t2;
    }
    return;
  }

  //
  // R^N exponential map of a skew-symmetric tensor. Not implemented.
  // \param r \f$ r \in so(N) \f$
  // \return \f$ R = \exp R \f$ with \f$ R \in SO(N) \f$
  //
  template<typename T>
  Tensor<T>
  exp_skew_symmetric(Tensor<T> const & r)
  {
    // Check whether skew-symmetry holds
    assert(norm(r+transpose(r)) < machine_epsilon<T>());

    Index const
    N = r.get_dimension();

    Tensor<T>
    R = identity<T>(N);

    T
    theta = 0.0;

    switch (N) {

    default:
      std::cerr << "Exp of so(N) N != 2,3 not implemented." << std::endl;
      exit(1);
      break;

    case 3:
      theta = sqrt(r(2,1)*r(2,1)+r(0,2)*r(0,2)+r(1,0)*r(1,0));

      //Check whether norm == 0. If so, return identity.
      if (theta >= machine_epsilon<T>()) {
        R += sin(theta)/theta*r + (1.0-cos(theta))/(theta*theta)*r*r;
      }
      break;

    case 2:
      theta = r(1,0);

      T
      c = cos(theta);

      T
      s = sin(theta);

      R(0,0) = c;
      R(0,1) = -s;
      R(1,0) = s;
      R(1,1) = c;

      break;

    }

    return R;
  }

  //
  // R^N off-diagonal norm. Useful for SVD and other algorithms
  // that rely on Jacobi-type procedures.
  // \param A
  // \return \f$ \sqrt(\sum_i \sum_{j, j\neq i} a_{ij}^2) \f$
  //
  template<typename T>
  T
  norm_off_diagonal(Tensor<T> const & A)
  {
    Index const
    N = A.get_dimension();

    T s = 0.0;

    switch (N) {

    default:
      for (Index i = 0; i < N; ++i) {
        for (Index j = 0; j < N; ++j) {
          if (i != j) s += A(i,j)*A(i,j);
        }
      }
      break;

    case 3:
      s = A(0,1)*A(0,1) + A(0,2)*A(0,2) + A(1,2)*A(1,2) +
          A(1,0)*A(1,0) + A(2,0)*A(2,0) + A(2,1)*A(2,1);
      break;

    case 2:
      s = A(0,1)*A(0,1) + A(1,0)*A(1,0);
      break;

    }

    return sqrt(s);
  }

  //
  // R^N arg max abs. Useful for inverse and other algorithms
  // that rely on Jacobi-type procedures.
  // \param A
  // \return \f$ (p,q) = arg max_{i,j} |a_{ij}| \f$
  //
  template<typename T>
  std::pair<Index, Index>
  arg_max_abs(Tensor<T> const & A)
  {
    Index p = 0;
    Index q = 0;

    T s = fabs(A(p,q));

    Index const
    N = A.get_dimension();

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        if (fabs(A(i,j)) > s) {
          p = i;
          q = j;
          s = fabs(A(i,j));
        }
      }
    }

    return std::make_pair(p,q);
  }

  //
  // R^N arg max off-diagonal. Useful for SVD and other algorithms
  // that rely on Jacobi-type procedures.
  // \param A
  // \return \f$ (p,q) = arg max_{i \neq j} |a_{ij}| \f$
  //
  template<typename T>
  std::pair<Index, Index>
  arg_max_off_diagonal(Tensor<T> const & A)
  {
    Index p = 0;
    Index q = 1;

    T s = fabs(A(p,q));

    Index const
    N = A.get_dimension();

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        if (i != j && fabs(A(i,j)) > s) {
          p = i;
          q = j;
          s = fabs(A(i,j));
        }
      }
    }

    return std::make_pair(p,q);
  }

  namespace {

    //
    // Singular value decomposition (SVD) for 2x2
    // bidiagonal matrix. Used for general 2x2 SVD.
    // Adapted from LAPAPCK's DLASV2, Netlib's dlasv2.c
    // and LBNL computational crystalography toolbox
    // \param f, g, h where A = [f, g; 0, h]
    // \return \f$ A = USV^T\f$
    //
    template<typename T>
    boost::tuple<Tensor<T>, Tensor<T>, Tensor<T> >
    svd_bidiagonal(T f, T g, T h)
    {
      T fa = std::abs(f);
      T ga = std::abs(g);
      T ha = std::abs(h);

      T s0 = 0.0;
      T s1 = 0.0;

      T cu = 1.0;
      T su = 0.0;
      T cv = 1.0;
      T sv = 0.0;

      bool swap_diag = (ha > fa);

      if (swap_diag) {
        std::swap(fa, ha);
        std::swap(f, h);
      }

      if (ga == 0.0) {
        s1 = ha;
        s0 = fa;
      } else if (ga > fa && fa / ga < machine_epsilon<T>()) {
        // case of very large ga
        s0 = ga;
        s1 = ha > 1.0 ?
            T(fa / (ga / ha)) :
            T((fa / ga) * ha);
        cu = 1.0;
        su = h / g;
        cv = f / g;
        sv = 1.0;
      } else {
        // normal case
        T d = fa - ha;
        T l = d != fa ?
            T(d / fa) :
            T(1.0); // l \in [0,1]
        T m = g / f; // m \in (-1/macheps, 1/macheps)
        T t = 2.0 - l; // t \in [1,2]
        T mm = m * m;
        T tt = t * t;
        T s = sqrt(tt + mm); // s \in [1,1 + 1/macheps]
        T r = l != 0.0 ?
            T(sqrt(l * l + mm)) :
            T(fabs(m)); // r \in [0,1 + 1/macheps]
        T a = 0.5 * (s + r); // a \in [1,1 + |m|]
        s1 = ha / a;
        s0 = fa * a;

        // Compute singular vectors
        T tau; // second assignment to T in DLASV2
        if (mm != 0.0) {
          tau = (m / (s + t) + m / (r + l)) * (1.0 + a);
        } else {
          // note that m is very tiny
          tau = l == 0.0 ?
              T(copysign(T(2.0), f) * copysign(T(1.0), g)) :
              T(g / copysign(d, f) + m / t);
        }
        T lv = sqrt(tau * tau + 4.0); // second assignment to L in DLASV2
        cv = 2.0 / lv;
        sv = tau / lv;
        cu = (cv + sv * m) / a;
        su = (h / f) * sv / a;
      }

      // Fix signs of singular values in accordance to sign of singular vectors
      s0 = copysign(s0, f);
      s1 = copysign(s1, h);

      if (swap_diag) {
        std::swap(cu, sv);
        std::swap(su, cv);
      }

      Tensor<T> U(cu, -su, su, cu);

      Tensor<T> S(s0, 0.0, 0.0, s1);

      Tensor<T> V(cv, -sv, sv, cv);

      return boost::make_tuple(U, S, V);
    }

    //
    // R^2 singular value decomposition (SVD)
    // \param A tensor
    // \return \f$ A = USV^T\f$
    //
    template<typename T>
    boost::tuple<Tensor<T>, Tensor<T>, Tensor<T> >
    svd_2x2(Tensor<T> const & A)
    {
      assert(A.get_dimension() == 2);

      // First compute a givens rotation to eliminate 1,0 entry in tensor
      T c = 1.0;
      T s = 0.0;
      boost::tie(c, s) = givens(A(0,0), A(1,0));

      Tensor<T>
      R(c, -s, s, c);

      Tensor<T>
      B = R * A;

      // B is bidiagonal. Use specialized algorithm to compute its SVD
      Tensor<T>
      X(2), S(2), V(2);

      boost::tie(X, S, V) = svd_bidiagonal(B(0,0), B(0,1), B(1,1));

      // Complete general 2x2 SVD with givens rotation calculated above
      Tensor<T>
      U = transpose(R) * X;

      return boost::make_tuple(U, S, V);
    }

    //
    // R^N singular value decomposition (SVD)
    // \param A tensor
    // \return \f$ A = USV^T\f$
    //
    template<typename T>
    boost::tuple<Tensor<T>, Tensor<T>, Tensor<T> >
    svd_NxN(Tensor<T> const & A)
    {
      Tensor<T>
      S = A;

      Index const
      N = A.get_dimension();

      Tensor<T>
      U = identity<T>(N);

      Tensor<T>
      V = identity<T>(N);

      T
      off = norm_off_diagonal(S);

      const T
      tol = machine_epsilon<T>() * norm(A);

      Index const
      max_iter = 128;

      Index
      num_iter = 0;

      while (off > tol && num_iter < max_iter) {

        // Find largest off-diagonal entry
        Index
        p = 0;

        Index
        q = 0;

        boost::tie(p,q) = arg_max_off_diagonal(S);

        if (p > q) {
          std::swap(p, q);
        }

        // Obtain left and right Givens rotations by using 2x2 SVD
        Tensor <T>
        Spq(S(p,p), S(p,q), S(q,p), S(q,q));

        Tensor <T>
        L(2), D(2), R(2);

        boost::tie(L, D, R) = svd_2x2(Spq);

        T const &
        cl = L(0,0);

        T const &
        sl = L(0,1);

        T const &
        cr = R(0,0);

        T const &
        sr = (sgn(R(0,1)) == sgn(R(1,0))) ? T(-R(0,1)) : T(R(0,1));

        // Apply both Givens rotations to matrices
        // that are converging to singular values and singular vectors
        givens_left(cl, sl, p, q, S);
        givens_right(cr, sr, p, q, S);

        givens_right(cl, sl, p, q, U);
        givens_left(cr, sr, p, q, V);

        off = norm_off_diagonal(S);
        num_iter++;
      }

      if (num_iter == max_iter) {
        std::cerr << "WARNING: SVD iteration did not converge." << std::endl;
      }

      // Fix signs for entries in the diagonal matrix S
      // that are negative
      for (Index i = 0; i < N; ++i) {
        if (S(i,i) < 0.0) {
          S(i,i) = -S(i,i);
          for (Index j = 0; j < N; ++j) {
            U(j,i) = -U(j,i);
          }
        }
      }

      Vector<T> s(N);
      Tensor<T> P(N);

      boost::tie(s, P) = sort_permutation(diag(S));
      S = diag(s);
      U = U * P;
      V = V * P;

      return boost::make_tuple(U, diag(diag(S)), transpose(V));
    }

  } // anonymous namespace

  //
  // R^N singular value decomposition (SVD)
  // \param A tensor
  // \return \f$ A = USV^T\f$
  //
  template<typename T>
  boost::tuple<Tensor<T>, Tensor<T>, Tensor<T> >
  svd(Tensor<T> const & A)
  {
    Index const
    N = A.get_dimension();

    Tensor<T>
    U(N), S(N), V(N);

    switch (N) {

    default:
      boost::tie(U, S, V) = svd_NxN(A);
      break;

    case 2:
      boost::tie(U, S, V) = svd_2x2(A);
      break;

    }

    return boost::make_tuple(U, S, V);
  }

  //
  // Project to O(N) (Orthogonal Group) using a Newton-type algorithm.
  // See Higham's Functions of Matrices p210 [2008]
  // \param A tensor (often a deformation-gradient-like tensor)
  // \return \f$ R = \argmin_Q \|A - Q\|\f$
  // This algorithm projects a given tensor in GL(N) to O(N).
  // The rotation/reflection obtained through this projection is
  // the orthogonal component of the real polar decomposition
  //
  template<typename T>
  Tensor<T>
  polar_rotation(Tensor<T> const & A)
  {
    Index const
    N = A.get_dimension();

    bool
    scale = true;

    const T
    tol_scale = 0.01;

    const T
    tol_conv = sqrt(N) * machine_epsilon<T>();

    Tensor<T>
    X = A;

    T
    gamma = 2.0;

    Index const
    max_iter = 128;

    Index
    num_iter = 0;

    while (num_iter < max_iter) {

      Tensor<T>
      Y = inverse(X);

      T
      mu = 1.0;

      if (scale == true) {
        mu = (norm_1(Y) * norm_infinity(Y)) / (norm_1(X) * norm_infinity(X));
        mu = sqrt(sqrt(mu));
      }

      Tensor<T>
      Z = 0.5 * (mu * X + transpose(Y) / mu);

      Tensor<T>
      D = Z - X;

      T
      delta = norm(D) / norm(Z);

      if (scale == true && delta < tol_scale) {
        scale = false;
      }

      bool
      end_iter =
          norm(D) <= sqrt(tol_conv) ||
          (delta > 0.5 * gamma && scale == false);

      X = Z;
      gamma = delta;

      if (end_iter == true) {
        break;
      }

      num_iter++;

    }

    if (num_iter == max_iter) {
      std::cerr << "WARNING: Polar iteration did not converge." << std::endl;
    }

    return X;
  }

  //
  // R^N Left polar decomposition
  // \param A tensor (often a deformation-gradient-like tensor)
  // \return \f$ VR = A \f$ with \f$ R \in SO(N) \f$ and \f$ V \in SPD(N) \f$
  //
  template<typename T>
  std::pair<Tensor<T>, Tensor<T> >
  polar_left(Tensor<T> const & A)
  {
    Tensor<T>
    R = polar_rotation(A);

    Tensor<T>
    V = symm(A * transpose(R));

    return std::make_pair(V, R);
  }

  //
  // R^N Right polar decomposition
  // \param A tensor (often a deformation-gradient-like tensor)
  // \return \f$ RU = A \f$ with \f$ R \in SO(N) \f$ and \f$ U \in SPD(N) \f$
  //
  template<typename T>
  std::pair<Tensor<T>, Tensor<T> >
  polar_right(Tensor<T> const & A)
  {
    Tensor<T>
    R = polar_rotation(A);

    Tensor<T>
    U = symm(transpose(R) * A);

    return std::make_pair(R, U);
  }

  //
  // R^3 left polar decomposition with eigenvalue decomposition
  // \param F tensor (often a deformation-gradient-like tensor)
  // \return \f$ VR = F \f$ with \f$ R \in SO(3) \f$ and V SPD(3)
  //
  template<typename T>
  std::pair<Tensor<T>, Tensor<T> >
  polar_left_eig(Tensor<T> const & F)
  {
    assert(F.get_dimension() == 3);

    // set up return tensors
    Tensor<T>
    R(3);

    Tensor<T>
    V(3);

    // temporary tensor used to compute R
    Tensor<T>
    Vinv(3);

    // compute spd tensor
    Tensor<T>
    b = F * transpose(F);

    // get eigenvalues/eigenvectors
    Tensor<T>
    eVal(3);

    Tensor<T>
    eVec(3);

    boost::tie(eVec, eVal) = eig_spd(b);

    // compute sqrt() and inv(sqrt()) of eigenvalues
    Tensor<T>
    x = zero<T>(3);

    x(0,0) = sqrt(eVal(0,0));
    x(1,1) = sqrt(eVal(1,1));
    x(2,2) = sqrt(eVal(2,2));

    Tensor<T>
    xi = zero<T>(3);

    xi(0,0) = 1.0 / x(0,0);
    xi(1,1) = 1.0 / x(1,1);
    xi(2,2) = 1.0 / x(2,2);

    // compute V, Vinv, and R
    V    = eVec * x * transpose(eVec);
    Vinv = eVec * xi * transpose(eVec);
    R    = Vinv * F;

    return std::make_pair(V, R);
  }

  //
  // R^3 right polar decomposition with eigenvalue decomposition
  // \param F tensor (often a deformation-gradient-like tensor)
  // \return \f$ RU = F \f$ with \f$ R \in SO(3) \f$ and U SPD(3)
  //
  template<typename T>
  std::pair<Tensor<T>, Tensor<T> >
  polar_right_eig(Tensor<T> const & F)
  {
    assert(F.get_dimension() == 3);

    Tensor<T>
    R(3);

    Tensor<T>
    U(3);

    // temporary tensor used to compute R
    Tensor<T>
    Uinv(3);

    // compute spd tensor
    Tensor<T>
    C = transpose(F) * F;

    // get eigenvalues/eigenvectors
    Tensor<T>
    eVal(3);

    Tensor<T>
    eVec(3);

    boost::tie(eVec, eVal) = eig_spd(C);

    // compute sqrt() and inv(sqrt()) of eigenvalues
    Tensor<T>
    x = zero<T>(3);

    x(0,0) = sqrt(eVal(0,0));
    x(1,1) = sqrt(eVal(1,1));
    x(2,2) = sqrt(eVal(2,2));

    Tensor<T>
    xi = zero<T>(3);

    xi(0,0) = 1.0 / x(0,0);
    xi(1,1) = 1.0 / x(1,1);
    xi(2,2) = 1.0 / x(2,2);

    // compute U, Uinv, and R
    U    = eVec * x * transpose(eVec);
    Uinv = eVec * xi * transpose(eVec);
    R    = F * Uinv;

    return std::make_pair(R, U);
  }

  //
  // R^N left polar decomposition with matrix logarithm for V
  // \param F tensor (often a deformation-gradient-like tensor)
  // \return \f$ VR = F \f$ with \f$ R \in SO(N) \f$ and V SPD(N), and log V
  //
  template<typename T>
  boost::tuple<Tensor<T>, Tensor<T>, Tensor<T> >
  polar_left_logV(Tensor<T> const & F)
  {
    Index const
    N = F.get_dimension();

    Tensor<T>
    X(N), S(N), Y(N);

    boost::tie(X, S, Y) = svd(F);

    Tensor<T>
    R = X * transpose(Y);

    Tensor<T>
    V = X * S * transpose(X);

    Tensor<T>
    s = S;

    for (Index i = 0; i < N; ++i) {
      s(i,i) = std::log(s(i,i));
    }

    Tensor<T>
    v = X * s * transpose(X);

    return boost::make_tuple(V, R, v);
  }

  //
  // R^N left polar decomposition with matrix logarithm for V
  // \param F tensor (often a deformation-gradient-like tensor)
  // \return \f$ VR = F \f$ with \f$ R \in SO(N) \f$ and V SPD(N), and log V
  //
  template<typename T>
  boost::tuple<Tensor<T>, Tensor<T>, Tensor<T> >
  polar_left_logV_lame(Tensor<T> const & F)
  {
    Index const
    N = F.get_dimension();

    // set up return tensors
    Tensor<T> R(N), V(N), v(N), Vinv(N);

    // compute spd tensor
    Tensor<T> b = F*transpose(F);

    // get eigenvalues/eigenvectors
    Tensor<T> eVal(N);
    Tensor<T> eVec(N);
    boost::tie(eVec,eVal) = eig_spd_cos(b);

    // compute sqrt() and inv(sqrt()) of eigenvalues
    Tensor<T> x = zero<T>(3);
    x(0,0) = sqrt(eVal(0,0));
    x(1,1) = sqrt(eVal(1,1));
    x(2,2) = sqrt(eVal(2,2));
    Tensor<T> xi = zero<T>(3);
    xi(0,0) = 1.0/x(0,0);
    xi(1,1) = 1.0/x(1,1);
    xi(2,2) = 1.0/x(2,2);
    Tensor<T> lnx = zero<T>(3);
    lnx(0,0) = std::log(x(0,0));
    lnx(1,1) = std::log(x(1,1));
    lnx(2,2) = std::log(x(2,2));

    // compute V, Vinv, log(V)=v, and R
    V    = eVec*x*transpose(eVec);
    Vinv = eVec*xi*transpose(eVec);
    v    = eVec*lnx*transpose(eVec);
    R    = Vinv*F;

    return boost::make_tuple(V,R,v);
  }

  //
  // R^N logarithmic map using BCH expansion (4 terms)
  // \param x tensor
  // \param y tensor
  // \return Baker-Campbell-Hausdorff series up to 4 terms
  //
  template<typename T>
  Tensor<T>
  bch(Tensor<T> const & x, Tensor<T> const & y)
  {
    return
        // first order term
        x + y
        +
        // second order term
        0.5*(x*y - y*x)
        +
        // third order term
        1.0/12.0 *
          (x*x*y - 2.0*x*y*x + x*y*y + y*x*x - 2.0*y*x*y + y*y*x)
        +
        // fourth order term
        1.0/24.0 *
        (x*x*y*y - 2.0*x*y*x*y + 2.0*y*x*y*x - y*y*x*x);
  }

  //
  // Symmetric Schur algorithm for R^2.
  // \param \f$ A = [f, g; g, h] \in S(2) \f$
  // \return \f$ c, s \rightarrow [c, -s; s, c]\f diagonalizes A$
  //
  template<typename T>
  std::pair<T, T>
  schur_sym(const T f, const T g, const T h)
  {
    T c = 1.0;
    T s = 0.0;

    if (g != 0.0) {
      T t = (h - f) / (2.0 * g);

      if (t >= 0.0) {
        t = 1.0 / (sqrt(1.0 + t * t) + t);
      } else {
        t = -1.0 / (sqrt(1.0 + t * t) - t);
      }
      c = 1.0 / sqrt(1.0 + t * t);
      s = t * c;
    }

    return std::make_pair(c, s);
  }

  //
  // Givens rotation. [c, -s; s, c] [a; b] = [r; 0]
  // \param a, b
  // \return c, s
  //
  template<typename T>
  std::pair<T, T>
  givens(T const & a, T const & b)
  {
    T c = 1.0;
    T s = 0.0;

    if (b != 0.0) {
      if (fabs(b) > fabs(a)) {
        const T t = - a / b;
        s = 1.0 / sqrt(1.0 + t * t);
        c = t * s;
      } else {
        const T t = - b / a;
        c = 1.0 / sqrt(1.0 + t * t);
        s = t * c;
      }
    }

    return std::make_pair(c, s);
  }

  namespace {

    //
    // R^N eigenvalue decomposition for symmetric 2nd-order tensor
    // \param A tensor
    // \return V eigenvectors, D eigenvalues in diagonal Matlab-style
    //
    template<typename T>
    std::pair<Tensor<T>, Tensor<T> >
    eig_sym_NxN(Tensor<T> const & A)
    {
      Tensor<T>
      D = symm(A);

      Index const
      N = A.get_dimension();

      Tensor<T>
      V = identity<T>(N);

      T
      off = norm_off_diagonal(D);

      const T
      tol = machine_epsilon<T>() * norm(A);

      Index const
      max_iter = 128;

      Index
      num_iter = 0;

      while (off > tol && num_iter < max_iter) {

        // Find largest off-diagonal entry
        Index
        p = 0;

        Index
        q = 0;

        boost::tie(p,q) = arg_max_off_diagonal(D);
        if (p > q) {
          std::swap(p,q);
        }

        // Obtain Givens rotations by using 2x2 symmetric Schur algorithm
        T const &
        f = D(p,p);

        T const &
        g = D(p,q);

        T const &
        h = D(q,q);

        T
        c, s;

        boost::tie(c, s) = schur_sym(f, g, h);

        // Apply Givens rotation to matrices
        // that are converging to eigenvalues and eigenvectors
        givens_left(c, s, p, q, D);
        givens_right(c, s, p, q, D);

        givens_right(c, s, p, q, V);

        off = norm_off_diagonal(D);
        num_iter++;
      }

      if (num_iter == max_iter) {
        std::cerr << "WARNING: EIG iteration did not converge." << std::endl;
      }

      Vector<T> d(N);
      Tensor<T> P(N);

      boost::tie(d, P) = sort_permutation(diag(D));
      D = diag(d);
      V = V * P;

      return std::make_pair(V, D);
    }

    //
    // R^2 eigenvalue decomposition for symmetric 2nd-order tensor
    // \param A tensor
    // \return V eigenvectors, D eigenvalues in diagonal Matlab-style
    //
    template<typename T>
    std::pair<Tensor<T>, Tensor<T> >
    eig_sym_2x2(Tensor<T> const & A)
    {
      assert(A.get_dimension() == 2);

      const T f = A(0,0);
      const T g = 0.5 * (A(0,1) + A(1,0));
      const T h = A(1,1);

      //
      // Eigenvalues, based on LAPACK's dlae2
      //
      const T sum = f + h;
      const T dif = fabs(f - h);
      const T g2 = fabs(g + g);

      T fhmax = f;
      T fhmin = h;

      const bool swap_diag = fabs(h) > fabs(f);

      if (swap_diag == true) {
        std::swap(fhmax, fhmin);
      }

      T r = 0.0;
      if (dif > g2) {
        const T t = g2 / dif;
        r = dif * sqrt(1.0 + t * t);
      } else if (dif < g2) {
        const T t = dif / g2;
        r = g2 * sqrt(1.0 + t * t);
      } else {
        // dif == g2, including zero
        r = g2 * sqrt(2.0);
      }

      T s0 = 0.0;
      T s1 = 0.0;

      if (sum != 0.0) {
        s0 = 0.5 * (sum + copysign(r, sum));
        // Order of execution important.
        // To get fully accurate smaller eigenvalue,
        // next line needs to be executed in higher precision.
        s1 = (fhmax / s0) * fhmin - (g / s0) * g;
      } else {
        // s0 == s1, including zero
        s0 = 0.5 * r;
        s1 = -0.5 * r;
      }

      Tensor<T>
      D(s0, 0.0, 0.0, s1);

      //
      // Eigenvectors
      //
      T
      c, s;

      boost::tie(c, s) = schur_sym(f, g, h);

      Tensor<T>
      V(c, -s, s, c);

      if (swap_diag == true) {
        // swap eigenvectors if eigenvalues were swapped
        std::swap(V(0,0), V(0,1));
        std::swap(V(1,0), V(1,1));
      }

      return std::make_pair(V, D);
    }

  } // anonymous namespace

  //
  // R^N eigenvalue decomposition for symmetric 2nd-order tensor
  // \param A tensor
  // \return V eigenvectors, D eigenvalues in diagonal Matlab-style
  //
  template<typename T>
  std::pair<Tensor<T>, Tensor<T> >
  eig_sym(Tensor<T> const & A)
  {
    Index const
    N = A.get_dimension();

    Tensor<T>
    V(N), D(N);

    switch (N) {

    default:
      boost::tie(V, D) = eig_sym_NxN(A);
      break;

    case 2:
      boost::tie(V, D) = eig_sym_2x2(A);
      break;

    }

    return std::make_pair(V, D);
  }

  //
  // R^N eigenvalue decomposition for SPD 2nd-order tensor
  // \param A tensor
  // \return V eigenvectors, D eigenvalues in diagonal Matlab-style
  //
  template<typename T>
  std::pair<Tensor<T>, Tensor<T> >
  eig_spd(Tensor<T> const & A)
  {
    return eig_sym(A);
  }

  //
  // R^3 eigenvalue decomposition for SPD 2nd-order tensor
  // \param A tensor
  // \return V eigenvectors, D eigenvalues in diagonal Matlab-style
  //
  template<typename T>
  std::pair<Tensor<T>, Tensor<T> >
  eig_spd_cos(Tensor<T> const & A)
  {
    assert(A.get_dimension() == 3);

    // This algorithm comes from the journal article
    // Scherzinger and Dohrmann, CMAME 197 (2008) 4007-4015

    // this algorithm will return the eigenvalues in D
    // and the eigenvectors in V
    Tensor<T>
    D = zero<T>(3);

    Tensor<T>
    V = zero<T>(3);

    // not sure if this is necessary...
    T
    pi = acos(-1);

    // convenience operators
    const Tensor<T>
    I = identity<T>(3);

    int
    ii[3][2] = { { 1, 2 }, { 2, 0 }, { 0, 1 } };

    Tensor<T>
    rm = zero<T>(3);

    // scale the matrix to reduce the characteristic equation
    T
    trA = (1.0/3.0) * I1(A);

    Tensor<T>
    Ap(A - trA*I);

    // compute other invariants
    T
    J2 = I2(Ap);

    T
    J3 = det(Ap);

    // deal with volumetric tensors
    if (-J2 <= 1.e-30)
    {
      D(0,0) = trA;
      D(1,1) = trA;
      D(2,2) = trA;

      V(0,0) = 1.0;
      V(1,0) = 0.0;
      V(2,0) = 0.0;

      V(0,1) = 0.0;
      V(1,1) = 1.0;
      V(2,1) = 0.0;

      V(0,2) = 0.0;
      V(1,2) = 0.0;
      V(2,2) = 1.0;
    }
    else
    {
      // first things first, find the most dominant e-value
      // Need to solve cos(3 theta)=rhs for theta
      T
      t1 = 3.0 / -J2;

      T
      rhs = (J3 / 2.0) * T(sqrt(t1 * t1 * t1));

      T
      theta = pi / 2.0 * (1.0 - (rhs < 0 ? -1.0 : 1.0));

      if (fabs(rhs) <= 1.0) theta = acos(rhs);

      T
      thetad3 = theta / 3.0;

      if (thetad3 > pi / 6.0) thetad3 += 2.0 * pi / 3.0;

      // most dominant e-value
      D(2,2) = 2.0 * cos(thetad3) * sqrt(-J2 / 3.0);

      // now reduce the system
      Tensor<T>
      R = Ap - D(2,2) * I;

      // QR factorization with column pivoting
      Vector<T> a(3);
      a(0) = R(0,0)*R(0,0) + R(1,0)*R(1,0) + R(2,0)*R(2,0);
      a(1) = R(0,1)*R(0,1) + R(1,1)*R(1,1) + R(2,1)*R(2,1);
      a(2) = R(0,2)*R(0,2) + R(1,2)*R(1,2) + R(2,2)*R(2,2);

      // find the most dominant column
      int k = 0;
      T max = a(0);
      if (a(1) > max)
      {
        k = 1;
        max = a(1);
      }
      if (a(2) > max)
      {
        k = 2;
      }

      // normalize the most dominant column to get s1
      a(k) = sqrt(a(k));
      for (int i(0); i < 3; ++i)
        R(i,k) /= a(k);

      // dot products of dominant column with other two columns
      T d0 = 0.0;
      T d1 = 0.0;
      for (int i(0); i < 3; ++i)
      {
        d0 += R(i,k) * R(i,ii[k][0]);
        d1 += R(i,k) * R(i,ii[k][1]);
      }

      // projection
      for (int i(0); i < 3; ++i)
      {
        R(i,ii[k][0]) -= d0 * R(i,k);
        R(i,ii[k][1]) -= d1 * R(i,k);
      }

      // now finding next most dominant column
      a.clear();
      for (int i(0); i < 3; ++i)
      {
        a(0) += R(i,ii[k][0]) * R(i,ii[k][0]);
        a(1) += R(i,ii[k][1]) * R(i,ii[k][1]);
      }

      int p = 0;
      if (fabs(a(1)) > fabs(a(0))) p = 1;

      // normalize next most dominant column to get s2
      a(p) = sqrt(a(p));
      int k2 = ii[k][p];

      for (int i(0); i < 3; ++i)
        R(i,k2) /= a(p);

      // set first eigenvector as cross product of s1 and s2
      V(0,2) = R(1,k) * R(2,k2) - R(2,k) * R(1,k2);
      V(1,2) = R(2,k) * R(0,k2) - R(0,k) * R(2,k2);
      V(2,2) = R(0,k) * R(1,k2) - R(1,k) * R(0,k2);

      // normalize
      T
      mag = sqrt(V(0,2) * V(0,2) + V(1,2) * V(1,2) + V(2,2) * V(2,2));

      V(0,2) /= mag;
      V(1,2) /= mag;
      V(2,2) /= mag;

      // now for the other two eigenvalues, extract vectors
      Vector<T>
      rk(R(0,k), R(1,k), R(2,k));

      Vector<T>
      rk2(R(0,k2), R(1,k2), R(2,k2));

      // compute projections
      Vector<T>
      ak = Ap * rk;

      Vector<T>
      ak2 = Ap * rk2;

      // set up reduced remainder matrix
      rm(0,0) = dot(rk,ak);
      rm(0,1) = dot(rk,ak2);
      rm(1,1) = dot(rk2,ak2);

      // compute eigenvalues 2 and 3
      T
      b = 0.5 * (rm(0,0) - rm(1,1));

      T
      fac = (b < 0 ? -1.0 : 1.0);

      T
      arg = b * b + rm(0,1) * rm(0,1);

      if (arg == 0)
        D(0,0) = rm(1,1) + b;
      else
        D(0,0) = rm(1,1) + b - fac * sqrt(b * b + rm(0,1) * rm(0,1));

      D(1,1) = rm(0,0) + rm(1,1) - D(0,0);

      // update reduced remainder matrix
      rm(0,0) -= D(0,0);
      rm(1,0) = rm(0,1);
      rm(1,1) -= D(0,0);

      // again, find most dominant column
      a.clear();
      a(0) = rm(0,0) * rm(0,0) + rm(0,1) * rm(0,1);
      a(1) = rm(0,1) * rm(0,1) + rm(1,1) * rm(1,1);

      int k3 = 0;
      if (a(1) > a(0)) k3 = 1;
      if (a(k3) == 0.0)
      {
        rm(0,k3) = 1.0;
        rm(1,k3) = 0.0;
      }

      // set 2nd eigenvector via cross product
      V(0,0) = rm(0,k3) * rk2(0) - rm(1,k3) * rk(0);
      V(1,0) = rm(0,k3) * rk2(1) - rm(1,k3) * rk(1);
      V(2,0) = rm(0,k3) * rk2(2) - rm(1,k3) * rk(2);

      // normalize
      mag = sqrt(V(0,0) * V(0,0) + V(1,0) * V(1,0) + V(2,0) * V(2,0));
      V(0,0) /= mag;
      V(1,0) /= mag;
      V(2,0) /= mag;

      // set last eigenvector as cross product of other two
      V(0,1) = V(1,0) * V(2,2) - V(2,0) * V(1,2);
      V(1,1) = V(2,0) * V(0,2) - V(0,0) * V(2,2);
      V(2,1) = V(0,0) * V(1,2) - V(1,0) * V(0,2);

      // normalize
      mag = sqrt(V(0,1) * V(0,1) + V(1,1) * V(1,1) + V(2,1) * V(2,1));
      V(0,1) /= mag;
      V(1,1) /= mag;
      V(2,1) /= mag;

      // add back in the offset
      for (int i(0); i < 3; ++i)
        D(i,i) += trA;
    }

    return std::make_pair(V, D);
  }

  //
  // Cholesky decomposition, rank-1 update algorithm
  // (Matrix Computations 3rd ed., Golub & Van Loan, p145)
  // \param A assumed symmetric tensor
  // \return G Cholesky factor A = GG^T
  // \return completed (bool) algorithm ran to completion
  //
  template<typename T>
  std::pair<Tensor<T>, bool >
  cholesky(Tensor<T> const & A)
  {
    Tensor<T>
    G = symm(A);

    Index const
    N = A.get_dimension();

    for (Index k = 0; k < N; ++k) {

      // Zeros above the diagonal
      for (Index j = k + 1; j < N; ++j) {
        G(k,j) = 0.0;
      }

      T
      s = G(k,k);

      if (s <= 0.0) {
        return std::make_pair(G, false);
      }

      s = sqrt(s);

      for (Index j = k + 1; j < N; ++j) {
        G(j,k) /= s;
      }

      G(k,k) = s;

      for (Index j = k + 1; j < N; ++j) {
        for (Index i = j; i < N; ++i) {
          G(i,j) -= G(i,k) * G(j,k);
        }
      }

    }

    return std::make_pair(G, true);
  }

} // namespace Intrepid

#endif // Intrepid_MiniTensor_LinearAlgebra_t_cc
