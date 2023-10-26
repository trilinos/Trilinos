#ifndef STK_MIDDLE_GRID_MATRIX_H
#define STK_MIDDLE_GRID_MATRIX_H

#include <cassert>
#include <ostream>
#include <stdexcept>
#include <vector>

#include <iostream>

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

// defines a dynamically sized matrix type
// Data is stores in row major order

template <typename T>
class Matrix
{
  public:
    explicit Matrix(const int m = 0, const int n = 0)
      : m_m(m)
      , m_n(n)
      , m_data(m * n)
    {}

    explicit Matrix(const int m, const int n, const std::vector<T>& data)
      : m_m(m)
      , m_n(n)
      , m_data(data)
    {
      assert(data.size() >= static_cast<unsigned>(m * n));
      m_data.resize(m * n);
    }

    T& operator()(const int i, const int j) { return m_data[get_idx(i, j)]; }

    const T& operator()(const int i, const int j) const { return m_data[get_idx(i, j)]; }

    int extent(const int dim) const
    {
      if (dim == 0)
        return m_m;
      else if (dim == 1)
        return m_n;
      else
        throw std::invalid_argument("matrix only has 2 dimensions");
    }

    void fill(const T& val)
    {
      for (auto& v : m_data)
        v = val;
    }

    int extent0() const { return m_m; }

    int extent1() const { return m_n; }

    // resizes matrix.  No guarantees about what the matrix contains after
    // this is done
    void resize(const int m, const int n)
    {
      m_data.resize(m * n);
      m_m = m;
      m_n = n;
    }

    // returns a pointer to the (contiguous) underlying data.
    // Prefer operator(), this is here for interfacing with BLAS/LAPACK
    T* data() { return m_data.data(); }

    const T* data() const { return m_data.data(); }

  private:
    // computes linear index
    int get_idx(const int i, const int j) const
    {
      assert(i >= 0 && i < m_m);
      assert(j >= 0 && j < m_n);
      return i * m_n + j;
    }

    int m_m;
    int m_n;
    std::vector<T> m_data;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& a)
{
  os << "Matrix: " << a.extent0() << " x " << a.extent1() << std::endl;
  for (int i = 0; i < a.extent0(); ++i)
  {
    os << "[";
    for (int j = 0; j < a.extent1(); ++j)
    {
      os << a(i, j);
      if (j < a.extent1() - 1)
        os << ", ";
    }

    os << "]" << std::endl;
  }

  return os;
}

template <typename T>
Matrix<T> transpose(const Matrix<T>& a)
{
  Matrix<T> b(a.extent1(), a.extent0());

  for (int i = 0; i < a.extent0(); ++i)
    for (int j = 0; j < a.extent1(); ++j)
      b(j, i) = a(i, j);

  return b;
}

enum class BlasTrans
{
  NoTrans = 0,
  Trans,
  CTrans
};

enum class LapackTrans : char
{
  NoTrans = 'N',
  Trans   = 'T'
};

void matvec(const double alpha, const Matrix<double>& a, const double* x, const double beta, double* y,
            BlasTrans trans = BlasTrans::NoTrans);

// does matrix-vector multiplication for matrix in column major format
void matvec(BlasTrans trans, int m, int n, double alpha, const double* a, const double* x, double beta, double* y);

// note: B must be max(A.extent0, A.extent1) x nrhs, because it is used
//       to return the result
// recommend work be the same size as A
// A gets overwritten by this function
void solve_least_squares(Matrix<double>& a, Matrix<double>& b, Matrix<double>& work,
                         LapackTrans trans = LapackTrans::NoTrans);

// calls lapack dgeqrf to compute QR factorization. A is overwritten with R in the upper
// triangle and intermediate results needed to compute Q in the rest.
void compute_qr_factorization(Matrix<double>& a, Matrix<double>& work, double* tau);

// calls lapack dgeqp3 to compute rank revealing QR
void compute_rank_revealing_qr(Matrix<double>& a, Matrix<double>& work, double* tau, int* jpvt);

void solve_upper_triangular(Matrix<double>& a, double* rhs);

// given the result of compute_qr_factorization stored in A and tau, solves a linear
// system.
void solve_qr_factorization(Matrix<double>& a, Matrix<double>& work, double* tau, double* rhs);

// solves the qrp (rank-revealing QR) factorization
void solve_qrp_factorization(Matrix<double>& a, Matrix<double>& work, double* tau, int* jpvt, double* rhs);

// solves the qrp (rank-revealing QR) factorization, using only the first numSingularValues rows/columns of the R factor.
// this is roughly equivalent to solving a truncated SVD using the first numSingularValues
void solve_qrp_factorization(Matrix<double>& a, Matrix<double>& work, double* tau, int* jpvt, double* rhs, int numSingularVals);

void solve_linear_system(Matrix<double>& a, int* ipiv, double* rhs, int nrhs=1);

// computes the (full) SVD factorization
// A is m x n
// s is min(m, n)
// u is m x m
// vt is n x n
// work is scratch space for Lapack
// lworks is length of work
// Note that u and vt are column-major
void compute_svd_factorization(Matrix<double>& a, double* s, double* u, double* vt, double* work, int lwork);


} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif
