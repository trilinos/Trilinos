/* This is a reference implementation of block CRS (BCRS) operations for a
   structured mesh: matrix-vector product, block-tridiagonal factorization and
   solves along lines. Others may want to copy this file and add or modify
   implementations to performance test new KokkosKernels code.
     StructuredBlock represents a 3D mesh having ni, nj, nk cells in each
   dimension. Variable ordering is such that the k index is the fastest and the
   i index is slowest. Smoothing lines are built in the k direction.
     BlockCrsMatrix is a simple block CRS data structure.
     BlockTridiagMatrices holds the block tridiagonal matrices.
     This program is intended to be simple to build on testbeds. First, build
   Kokkos as a standalone library using Kokkos's generate_makefile.sh. Then
   build with something like these lines, where KH is the location of the
   standalone Kokkos installation:
     cpu: (KH=/home/ambradl/lib/kokkos/cpu; g++-4.7 -O3 -I$KH/include --std=c++11 -pedantic -Wall -fopenmp bcrs.cpp -L$KH/lib -lkokkos -ldl)
     gpu: (KH=/home/ambradl/lib/kokkos/gpu; $KH/bin/nvcc_wrapper -O3 -I$KH/include --std=c++11 -pedantic -Wall bcrs.cpp -L$KH/lib -lkokkos)

   An example run is
     ./a.out -ni 32 -nj 32 -nk 128 -bs 5 -c
   This runs a sequence of unit tests, then runs a problem having a 32x32x128
   structured block with the lines oriented along the third dimension (line
   length = 128). The block size is 5. -c adds a somewhat expensive check of the
   answer. It's good to run with -c once in a while, but the cheap unit tests
   that always run before the big problem already provide good coverage.

   todo
   - might set up a map of A block to T block so extract's index
     finding is faster
   - for <= 8x8 blocks, can do 2 factorizations per team
   - team abstraction/wrapper over Kokkos::TeamPolicy::member_type so i can do
     multiple useful teams per kokkos team cleanly.
       struct Team {
          int team_size, team_rank;
          void* wrk;
          void team_barrier();
       };
   - don't assume all lines are the same length
   - ray t.'s suggestion: block cyclic reduction. could do this down
     to the level necessary to make enough tridiags to solve.
 */

#include <cassert>
#include <cstring>
#include <vector>
#include <iostream>
#include <sstream>
#include <memory>
#include <limits>
#include <stdexcept>
#include <sys/time.h>

#define pr(m) do {                              \
    std::stringstream _ss_;                     \
    _ss_ << m << std::endl;                     \
    std::cerr << _ss_.str();                    \
  } while (0)
#define prc(m) pr(#m << " | " << (m))
#define puf(m) "(" << #m << " " << (m) << ")"
#define pu(m) << " " << puf(m)
template <typename T>
static void prarr (const std::string& name, const T* const v, const size_t n) {
  std::cerr << name << " = [";
  for (size_t i = 0; i < n; ++i) std::cerr << " " << v[i];
  std::cerr << "];\n";
}

#define TEST_ASSERT(m, success) \
  do { if ( ! (m)) { success = false; std::cerr << "FAILED: " << #m << "\n"; } } while (0)

static double urand () { return rand() / ((double) RAND_MAX + 1.0); }

#include <Kokkos_Core.hpp>
#define KIF KOKKOS_INLINE_FUNCTION
#define KFIF KOKKOS_FORCEINLINE_FUNCTION
namespace ko = Kokkos;

#if 0
KIF static unsigned int nextpow2 (unsigned int v) {
  // From https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}
#endif

static constexpr int max_team_size = 256;

static unsigned int calc_team_size (unsigned int v, const unsigned int unit = 32) {
#if 0
  if (v < unit) v = unit;
  v = nextpow2(v);
  v = std::min<unsigned int>(max_team_size, v);
#else
  v = std::min<unsigned int>(max_team_size, unit*((v + unit - 1)/unit));
#endif
  return v;
}

namespace impl {
template <typename T> KIF const T& min (const T& a, const T& b) { return a < b ? a : b; }
}

namespace vdec {
template <typename MemoryTraitsType, ko::MemoryTraitsFlags flag>
using MemoryTraits = ko::MemoryTraits<MemoryTraitsType::Unmanaged |
                                      MemoryTraitsType::RandomAccess |
                                      MemoryTraitsType::Atomic |
                                      flag>;
template <typename View>
using Unmanaged = ko::View<typename View::data_type, typename View::array_layout,
                           typename View::device_type,
                           MemoryTraits<typename View::memory_traits, ko::Unmanaged> >;
template <typename View>
using Const = ko::View<typename View::const_data_type, typename View::array_layout,
                       typename View::device_type, typename View::memory_traits>;
template <typename View>
using ConstUnmanaged = Const<Unmanaged<View> >;

template <typename ViewType>
KIF Unmanaged<ViewType> get_unmanaged (ViewType& vt) { return Unmanaged<ViewType>(vt); }
template <typename ViewType>
KIF Unmanaged<ViewType> get_const (ViewType& vt) { return Const<ViewType>(vt); }
} // namespace vdec

class SimpleTimer {
  std::string name_;
  timeval t1_;
  static double calc_et (const timeval& t1, const timeval& t2) {
    static const double us = 1.0e6;
    return (t2.tv_sec * us + t2.tv_usec - t1.tv_sec * us - t1.tv_usec) / us;
  }

public:
  SimpleTimer (const std::string& name)
    : name_(name)
  { gettimeofday(&t1_, 0); }

  double elapsed_time () const {
    timeval t2;
    gettimeofday(&t2, 0);
    return calc_et(t1_, t2);
  }

  ~SimpleTimer () {
    ko::fence();
    double et = elapsed_time();
    std::cout << "    + Timer: " << name_ << " = " << et << "\n";
  }
};

template <typename T> KIF constexpr T square (const T& x) { return x*x; }

typedef int Int;
typedef int Size;
typedef double Scalar;

template <typename T, typename Int>
inline void copy (T* const d, const T* const s, const Int& n) {
  memcpy(d, s, n*sizeof(T));
}

template <typename ExeSpace = ko::DefaultExecutionSpace,
          typename IndexT = Int, typename SizeT = Size>
struct CrsGraph {
  typedef IndexT Index;
  typedef SizeT Size;
  typedef ko::View<Size*, ExeSpace> Rowptr;
  typedef ko::View<Index*, ExeSpace> Colidx;

  Rowptr rowptr;
  Colidx colidx;

  CrsGraph () : rowptr("rowptr", 1), colidx("colidx", 0) {}
  KIF bool empty () const { return rowptr.dimension_0() == 0 || colidx.dimension_0() == 0; }
  KIF Index nrow () const { return empty() ? 0 : static_cast<Index>(rowptr.dimension_0()) - 1; }
  KIF Size nnz () const { return empty() ? 0 : colidx.dimension_0(); }
};

namespace impl {
#ifdef KOKKOS_HAVE_CUDA
template <typename ExeSpace> using is_cuda = std::is_same<ExeSpace, ko::Cuda>;
#else
template <typename ExeSpace> struct is_cuda { enum { value = false }; };
#endif

#ifdef MIMIC_SPARC
enum { mimic_sparc = true };
#else
enum { mimic_sparc = false };
#endif

template <typename ExeSpace, typename Scalar> struct BcrsValuesType {
  typedef typename ko::View<Scalar***, ExeSpace> type;
  enum { is_LayoutStride = false };
};
#ifdef KOKKOS_HAVE_CUDA
template <typename Scalar> struct BcrsValuesType<ko::Cuda, Scalar> {
  typedef typename ko::View<Scalar***, ko::LayoutStride, ko::Cuda> type;
  enum { is_LayoutStride = true };
};
#endif
} // namespace impl

template <typename ExeSpace, typename Scalar> struct BcrsValuesType {
  typedef typename impl::BcrsValuesType<ExeSpace, Scalar>::type type;
  enum { is_LayoutStride = impl::BcrsValuesType<ExeSpace, Scalar>::is_LayoutStride };
};
template <typename Scalar> struct BcrsValuesType<ko::DefaultHostExecutionSpace, Scalar> {
  typedef typename impl::BcrsValuesType<ko::DefaultExecutionSpace, Scalar>::type::HostMirror type;
  enum { is_LayoutStride = impl::BcrsValuesType<ko::DefaultExecutionSpace, Scalar>::is_LayoutStride };
};

// Block size bs is a RT parameter, but some impl details make it CT.
template <typename ExeSpaceT = ko::DefaultExecutionSpace,
          typename ScalarT = Scalar, typename IndexT = Int, typename SizeT = Size>
class BcrsMatrix {
public:
  typedef ExeSpaceT ExeSpace;
  typedef IndexT Index;
  typedef SizeT Size;
  typedef ScalarT Scalar;
  typedef CrsGraph<ExeSpace, Index, Size> CrsGraphType;
  typedef typename BcrsValuesType<ExeSpace, Scalar>::type Values;

private:
  enum { is_LayoutStride = BcrsValuesType<ExeSpace, Scalar>::is_LayoutStride };
  const std::shared_ptr<const CrsGraphType> g_;
  const Int bs_;
  Values v_;

  template <bool is_LayoutStride>
  void init_v (const typename std::enable_if<is_LayoutStride, int>::type* = 0) {
    // Good layout for the MVP by row on the GPU.
    const int order[] = {1, 2, 0}, dims[] = {g_->nnz(), bs_, bs_};
    v_ = Values(ko::ViewAllocateWithoutInitializing("values"),
                ko::LayoutStride::order_dimensions(3, order, dims));    
  }
  template <bool is_LayoutStride>
  void init_v (const typename std::enable_if< ! is_LayoutStride, int>::type* = 0) {
    v_ = Values(ko::ViewAllocateWithoutInitializing("values"), g_->nnz(), bs_, bs_);
  }

public:
  BcrsMatrix (const std::shared_ptr<const CrsGraphType>& g, const Int block_size)
    : g_(g), bs_(block_size)
  { init_v<is_LayoutStride>(); }

  const Int& bs () const { return bs_; }
  const CrsGraphType& g () const { return *g_; }
  const Values& v () const { return v_; }
};

// Representation of a structured block mesh. The fastest index is k.
struct StructuredBlock {
  const Int ni, nj, nk;
  StructuredBlock (const Int ni, const Int nj, const Int nk)
    : ni(ni), nj(nj), nk(nk), njnk_(nj*nk) {}
  KIF Int size () const { return ni*njnk_; }
  KIF Int ijk2id (const Int i, const Int j, const Int k) const { return (i*nj + j)*nk + k; }
  KIF void id2ijk (const Int id, Int& i, Int& j, Int& k) const {
    i = id / njnk_;
    k = id % njnk_;
    j = k / nk;
    k = k % nk;
  }
private:
  const Int njnk_;
};

struct StencilShape { enum Enum { cross }; };

// Given a structured block and a stencil (at present, just a 3D 1-hop cross),
// construct a corresponding CRS graph.
std::shared_ptr<CrsGraph<ko::DefaultHostExecutionSpace> >
make_graph_for_structured_block (const StructuredBlock& sb, const StencilShape::Enum shape) {
  typedef CrsGraph<ko::DefaultHostExecutionSpace> CrsGraphT;
  auto g = std::make_shared<CrsGraphT>();
  ko::resize(g->rowptr, sb.size() + 1);
  g->rowptr[0] = 0;
  std::vector<typename CrsGraphT::Index> colidx;
  switch (shape) {
  case StencilShape::cross:
    for (Int c = 0; c < sb.size(); ++c) {
      Int i, j, k, n = 0;
      sb.id2ijk(c, i, j, k);
      colidx.push_back(c); ++n;
      if (i > 0)       { colidx.push_back(sb.ijk2id(i-1, j, k  )); ++n; }
      if (i+1 < sb.ni) { colidx.push_back(sb.ijk2id(i+1, j, k  )); ++n; }
      if (j > 0)       { colidx.push_back(sb.ijk2id(i, j-1, k  )); ++n; }
      if (j+1 < sb.nj) { colidx.push_back(sb.ijk2id(i, j+1, k  )); ++n; }
      if (k > 0)       { colidx.push_back(sb.ijk2id(i, j,   k-1)); ++n; }
      if (k+1 < sb.nk) { colidx.push_back(sb.ijk2id(i, j,   k+1)); ++n; }
      g->rowptr[c+1] = g->rowptr[c] + n;
    }
    break;
  }
  assert(g->rowptr[sb.size()] == static_cast<Int>(colidx.size()));
  for (Int c = 0; c < sb.size(); ++c)
    std::sort(colidx.begin() + g->rowptr[c], colidx.begin() + g->rowptr[c+1]);
  ko::resize(g->colidx, colidx.size());
  copy(g->colidx.data(), colidx.data(), colidx.size());
  return g;
}

// Fill a matrix with random numbers. The on-diag blocks are s.p.d.
void fill_matrix (BcrsMatrix<ko::DefaultHostExecutionSpace>& m) {
#ifdef FAST_MATRIX_FILL
#pragma message "BUILD WITH FAST_MATRIX_FILL"
  const auto bs = m.bs();
  Scalar diag_block[bs][bs], offdiag_block[bs][bs];
  
  // precompute 
  {
    Scalar tmp[bs*bs];
    // If this is a diagonal block, make spd.
    for (Int i = 0; i < bs*bs; ++i)
      tmp[i] = 2*(urand() - 0.5);
    for (Int i = 0; i < bs; ++i)
      for (Int j = i; j < bs; ++j) {
        diag_block[i][j] = 0;
        for (Int k = 0; k < bs; ++k)
          diag_block[i][j] += tmp[i*bs + k] * tmp[j*bs + k];
        if (i != j) diag_block[j][i] = diag_block[i][j];
        else {
          // Improve conditioning.
          diag_block[i][j] *= 0.5*bs;
        }
      }
  }
  {
    for (Int i = 0; i < bs; ++i)
      for (Int j = 0; j < bs; ++j) 
        // Down-weight off-diag blocks to improve conditioning.
        offdiag_block[i][j] = 0.1 * 2*(urand() - 0.5);
  }
  const auto& g = m.g();
  for (Int r = 0; r < g.nrow(); ++r)
    for (Int ci = g.rowptr[r]; ci < g.rowptr[r+1]; ++ci) {
      auto b = ko::subview(m.v(), ci, ko::ALL(), ko::ALL());
      if (g.colidx[ci] == r) 
        for (Int i = 0; i < bs; ++i)
          for (Int j = 0; j < bs; ++j) 
            b(i,j) = diag_block[i][j];
      else 
        for (Int i = 0; i < bs; ++i)
          for (Int j = 0; j < bs; ++j) 
            b(i,j) = offdiag_block[i][j];
    }  
#else
  const auto bs = m.bs();
  std::vector<Scalar> tmp(square(bs));
  const auto& g = m.g();
  for (Int r = 0; r < g.nrow(); ++r)
    for (Int ci = g.rowptr[r]; ci < g.rowptr[r+1]; ++ci) {
      auto b = ko::subview(m.v(), ci, ko::ALL(), ko::ALL());
      if (g.colidx[ci] == r) {
        // If this is a diagonal block, make spd.
        for (Int i = 0; i < bs*bs; ++i)
          tmp[i] = 2*(urand() - 0.5);
        for (Int i = 0; i < bs; ++i)
          for (Int j = i; j < bs; ++j) {
            b(i,j) = 0;
            for (Int k = 0; k < bs; ++k)
              b(i,j) += tmp[i*bs + k] * tmp[j*bs + k];
            if (i != j) b(j,i) = b(i,j);
            else {
              // Improve conditioning.
              b(i,j) *= 0.5*bs;
            }
          }
      } else {
        for (Int i = 0; i < bs; ++i)
          for (Int j = 0; j < bs; ++j) {
            // Down-weight off-diag blocks to improve conditioning.
            b(i,j) = 0.1 * 2*(urand() - 0.5);
          }
      }
    }
#endif
}

std::shared_ptr<BcrsMatrix<ko::DefaultExecutionSpace> >
copy_to_device (const BcrsMatrix<ko::DefaultHostExecutionSpace>& mh) {
  auto g = std::make_shared<CrsGraph<ko::DefaultExecutionSpace> >();
  ko::resize(g->rowptr, mh.g().rowptr.size());
  ko::deep_copy(g->rowptr, mh.g().rowptr);
  ko::resize(g->colidx, mh.g().colidx.size());
  ko::deep_copy(g->colidx, mh.g().colidx);
  auto md = std::make_shared<BcrsMatrix<ko::DefaultExecutionSpace> >(g, mh.bs());
  ko::deep_copy(md->v(), mh.v());
  return md;
}

template <typename Matrix>
KIF void cholfac_serial (Matrix a) {
  const Int n = static_cast<Int>(a.dimension_0());
  for (Int k = 0; k < n; ++k) {
    a(k,k) = std::sqrt(a(k,k));
    const auto scale = 1 / a(k,k);
    for (Int j = k+1; j < n; ++j)
      a(k,j) *= scale;
    for (Int i = k+1; i < n; ++i) {
      const auto& aik = a(k,i);
      for (Int j = i; j < n; ++j)
        a(i,j) -= aik * a(k,j);
    }
  }
}

template <typename Matrix, typename MVecB, typename MVecX>
KIF void cholslv_serial (const Matrix a, const MVecB b, MVecX x) {
  const Int n = static_cast<Int>(a.dimension_0()), nrhs = static_cast<Int>(b.dimension_1());
  if (x.data() != b.data())
    for (Int irhs = 0; irhs < nrhs; ++irhs)
      for (Int i = 0; i < n; ++i)
        x(i,irhs) = b(i,irhs);
  for (Int j = 0; j < n; ++j) {
    const auto& ajj = a(j,j);
    for (Int irhs = 0; irhs < nrhs; ++irhs) {
      x(j,irhs) /= ajj;
      const auto& xj = x(j,irhs);
      for (Int i = j+1; i < n; ++i)
        x(i,irhs) -= a(j,i) * xj;
    }
  }
  for (Int irhs = 0; irhs < nrhs; ++irhs)
    for (Int i = n-1; i >= 0; --i) {
      for (Int j = i+1; j < n; ++j)
        x(i,irhs) -= a(i,j) * x(j,irhs);
      x(i,irhs) /= a(i,i);
    }
}

// LU, no pivot, in place.
template <typename Matrix>
KIF void lunpfac_serial (Matrix a) {
  const Int n = static_cast<Int>(a.dimension_0());
  for (Int k = 0; k < n; ++k) {
    const auto scale = 1 / a(k,k);
    for (Int i = k+1; i < n; ++i) {
      a(i,k) *= scale;
      const auto& aik = a(i,k);
      for (Int j = k+1; j < n; ++j)
        a(i,j) -= aik * a(k,j);
    }
  }
}

// LU, no pivot, in place.
template <typename ExeSpace, typename Matrix>
KIF void lunpfac (const typename ko::TeamPolicy<ExeSpace>::member_type& mbr, Matrix a) {
  const Int n = static_cast<Int>(a.dimension_0());
  const Int ncol_batch = impl::min<Int>(n, mbr.team_size() / n);
  const bool ut = mbr.team_rank() < n*ncol_batch;
  const Int i = mbr.team_rank() % n;
  const Int col_start = mbr.team_rank() / n;
  for (Int k = 0; k < n; ++k) {
    if (ut && i > k && col_start == 0)
      a(i,k) /= a(k,k);
    mbr.team_barrier();
    if (ut && i > k) {
      const auto& aik = a(i,k);
      for (Int j = col_start; j < n; j += ncol_batch)
        if (j > k)
          a(i,j) -= aik * a(k,j);
    }
    mbr.team_barrier();
  }
}

// Solve L x = b, L(i,i) = 1. b can an alias of x.
template <typename Matrix, typename MVecB, typename MVecX>
KIF void lunpfwd_serial (const Matrix a, const MVecB b, MVecX x) {
  const Int n = static_cast<Int>(a.dimension_0()), nrhs = static_cast<Int>(b.dimension_1());
  for (Int i = 0; i < n; ++i)
    for (Int irhs = 0; irhs < nrhs; ++irhs) {
      auto accum = b(i,irhs);
      for (Int j = 0; j < i; ++j)
        accum -= a(i,j) * x(j,irhs);
      x(i,irhs) = accum;
    }  
}

static KIF Int calc_stride (const Int& n) { return n; /*nextpow2(n);*/ }

// Solve L x = b, L(i,i) = 1. b can an alias of x.
template <typename ExeSpace, typename Matrix, typename MVecB, typename MVecX>
KIF void lunpfwd (const typename ko::TeamPolicy<ExeSpace>::member_type& mbr,
                  const Matrix a, const MVecB b, MVecX x, const Int team_size,
                  const Int thread_id_base) {
  const Int n = static_cast<Int>(a.dimension_0()), nrhs = static_cast<Int>(b.dimension_1());
  const Int n_p2 = calc_stride(n);
  const Int nrhs_batch = impl::min<Int>(nrhs, team_size / n_p2);
  const Int nrhs_per_thread = (nrhs + nrhs_batch - 1) / nrhs_batch;
  const Int team_rank = mbr.team_rank() - thread_id_base;
  const Int row = team_rank % n_p2;
  const Int k = team_rank / n_p2;
  const bool use_thread = row < n && k < nrhs_batch;
  Int irhs = k;
  for (Int ctr = 0; ctr < nrhs_per_thread; ++ctr) {
    const bool ut = use_thread && irhs < nrhs;
    if (ut)
      x(row,irhs) = b(row,irhs);
    mbr.team_barrier();
    for (Int col = 0; col < n-1; ++col) {
      if (ut && row > col)
        x(row,irhs) -= x(col,irhs) * a(row,col);
      mbr.team_barrier();
    }
    irhs += nrhs_batch;
  }
}

// Solve U x = b. b can an alias of x.
template <typename Matrix, typename MVecB, typename MVecX>
KIF void lunpbwd_serial (const Matrix a, const MVecB b, MVecX x) {
  const Int n = static_cast<Int>(a.dimension_0()), nrhs = static_cast<Int>(b.dimension_1());
  for (Int irhs = 0; irhs < nrhs; ++irhs)
    for (Int i = n-1; i >= 0; --i) {
      auto accum = b(i,irhs);
      for (Int j = i+1; j < n; ++j)
        accum -= a(i,j) * x(j,irhs);
      x(i,irhs) = accum / a(i,i);
    }
}

// Solve U x = b. b can an alias of x.
template <typename ExeSpace, typename Matrix, typename MVecB, typename MVecX>
KIF void lunpbwd (const typename ko::TeamPolicy<ExeSpace>::member_type& mbr,
                  const Matrix a, const MVecB b, MVecX x, const Int team_size,
                  const Int thread_id_base) {
  const Int n = static_cast<Int>(a.dimension_0()), nrhs = static_cast<Int>(b.dimension_1());
  const Int n_p2 = calc_stride(n);
  const Int nrhs_batch = impl::min<Int>(nrhs, team_size / n_p2);
  const Int nrhs_per_thread = (nrhs + nrhs_batch - 1) / nrhs_batch;
  const Int team_rank = mbr.team_rank() - thread_id_base;
  const Int row = team_rank % n_p2;
  const Int k = team_rank / n_p2;
  const bool use_thread = row < n && k < nrhs_batch;
  Int irhs = k;
  for (Int ctr = 0; ctr < nrhs_per_thread; ++ctr) {
    const bool ut = use_thread && irhs < nrhs;
    if (ut)
      x(row,irhs) = b(row,irhs);
    mbr.team_barrier();
    for (Int col = n-1; col >= 0; --col) {
      if (ut && row == col)
        x(row,irhs) /= a(col,col);
      mbr.team_barrier();
      if (ut && row < col)
        x(row,irhs) -= x(col,irhs) * a(row,col);
    }
    irhs += nrhs_batch;
  }
}

// Solve using the LU factorization. b can an alias of x.
template <typename Matrix, typename MVecB, typename MVecX>
KIF void lunpslv_serial (const Matrix a, const MVecB b, MVecX x) {
  lunpfwd_serial(a, b, x);
  lunpbwd_serial(a, x, x);
}

// Solve U' X' = B'. B can be an alias of X.
template <typename Matrix, typename MVecB, typename MVecX>
KIF void lunpbwdtt_serial (const Matrix U, const MVecB B, MVecX X) {
  const Int n = static_cast<Int>(U.dimension_0()), nrhs = static_cast<Int>(B.dimension_0());
  for (Int k = 0; k < nrhs; ++k) {
    auto b = ko::subview(B, k, ko::ALL());
    auto x = ko::subview(X, k, ko::ALL());
    {
      auto u = ko::subview(U, 0, ko::ALL());
      x(0) = b(0) / U(0,0);
      const auto x0 = x(0);
      for (Int i = 1; i < n; ++i)
        x(i) = b(i) - u(i)*x0;
    }
    for (Int j = 1; j < n; ++j) {
      auto u = ko::subview(U, j, ko::ALL());
      x(j) /= U(j,j);
      const auto xj = x(j);
      for (Int i = j+1; i < n; ++i)
        x(i) -= u(i)*xj;
    }
  }
}

// Solve U' X' = B'. B can be an alias of X. X and B are square.
template <typename ExeSpace, typename Matrix, typename MVecB, typename MVecX>
KIF void lunpbwdtt (const typename ko::TeamPolicy<ExeSpace>::member_type& mbr,
                    const Matrix U, const MVecB B, MVecX X) {
  const Int n = static_cast<Int>(U.dimension_0());
  const Int ncol_batch = impl::min<Int>(n, mbr.team_size() / n);
  const Int ncol_per_thread = (n + ncol_batch - 1) / ncol_batch;
  const Int k = mbr.team_rank() % n;
  const Int col_start = mbr.team_rank() / n;
  const bool ut = mbr.team_rank() < n*ncol_batch;
  // There are n separate trisolves. k indexes the trisolve. The solve
  // pattern is that of a CSC L matrix.
  if (ut && col_start == 0)
    X(k,0) = B(k,0) / U(0,0);
  mbr.team_barrier();
  if (ut) {
    const auto& Xk0 = X(k,0);
    Int i = col_start;
    for (Int ctr = 0; ctr < ncol_per_thread; ++ctr) {
      if (i > 0 && i < n)
        X(k,i) = B(k,i) - U(0,i)*Xk0;
      i += ncol_batch;
    }
  }
  mbr.team_barrier();
  for (Int j = 1; j < n; ++j) {
    if (ut && col_start == 0)
      X(k,j) /= U(j,j);
    mbr.team_barrier();
    if (ut) {
      const auto& Xkj = X(k,j);
      Int i = col_start;
      for (Int ctr = 0; ctr < ncol_per_thread; ++ctr) {
        if (ut && i > j && i < n)
          X(k,i) -= U(j,i)*Xkj;
        i += ncol_batch;
      }
    }
    mbr.team_barrier();
  }
}

// a := a - b*c for all nxn matrices.
template <typename MatrixA, typename ConstMatrixB, typename ConstMatrixC>
KIF void ngemm (MatrixA a, ConstMatrixB b, ConstMatrixC c) {
  const Int n = static_cast<Int>(a.dimension_0());
  for (Int i = 0; i < n; ++i) {
    auto ai = ko::subview(a, i, ko::ALL());
    auto bi = ko::subview(b, i, ko::ALL());
    for (Int k = 0; k < n; ++k) {
      const auto& bik = bi(k);
      auto ck = ko::subview(c, k, ko::ALL());
      for (Int j = 0; j < n; ++j)
        ai(j) -= bik * ck(j);
    }
  }
}

// d := a - b*c, for b nxn.
template <typename ConstMatrixA, typename ConstMatrixB, typename ConstMatrixC, typename MatrixD>
KIF void ngemm (ConstMatrixA a, ConstMatrixB b, ConstMatrixC c, MatrixD d) {
  const Int n = static_cast<Int>(a.dimension_0()), nrhs = static_cast<Int>(a.dimension_1());
  for (Int i = 0; i < n; ++i) {
    for (Int j = 0; j < nrhs; ++j)
      d(i,j) = a(i,j);
    for (Int k = 0; k < n; ++k) {
      const auto& bik = b(i,k);
      for (Int j = 0; j < nrhs; ++j)
        d(i,j) -= bik * c(k,j);
    }
  }
}

// d := a - b*c, for b nxn.
template <typename ExeSpace, typename ConstMatrixA, typename ConstMatrixB, typename ConstMatrixC, typename MatrixD>
KIF void ngemm (const typename ko::TeamPolicy<ExeSpace>::member_type& mbr,
                ConstMatrixA a, ConstMatrixB b, ConstMatrixC c, MatrixD d) {
  const Int n = static_cast<Int>(a.dimension_0()), nrhs = static_cast<Int>(a.dimension_1());
  const Int ncol_batch = impl::min<Int>(nrhs, mbr.team_size() / n);
  if (mbr.team_rank() < ncol_batch * n) {
    const Int i = mbr.team_rank() % n;
    Int j = mbr.team_rank() / n;
    while (j < nrhs) {
      auto accum = a(i,j);
      for (Int k = 0; k < n; ++k)
        accum -= b(i,k) * c(k,j);
      d(i,j) = accum;
      j += ncol_batch;
    }
  }
  mbr.team_barrier();
}

// Fuse lunpfwd and ngemm.
template <typename ExeSpace, typename Matrix, typename MVecX>
KIF void lunpfwd_ngemm (const typename ko::TeamPolicy<ExeSpace>::member_type& mbr,
                        const Matrix L, MVecX x, const Int team_size,
                        const Int thread_id_base) {
  Int n, nrhs, nrhs_batch, nrhs_per_thread, row, irhs, tdidx, row_minus_base;
  bool usethread;
  {
    n = static_cast<Int>(L.dimension_1());
    nrhs = static_cast<Int>(x.dimension_1());
    const Int n2 = n << 1;
    const Int n_p2 = calc_stride(n2);
    nrhs_batch = impl::min<Int>(nrhs, team_size / n_p2);
    nrhs_per_thread = (nrhs + nrhs_batch - 1) / nrhs_batch;
    const Int team_rank = mbr.team_rank() - thread_id_base;
    row = team_rank % n_p2;
    irhs = team_rank / n_p2;
    tdidx = row < n ? 0 : 1;
    row_minus_base = row - (row < n ? 0 : n);
    usethread = row < n2 && irhs < nrhs_batch;
  }
  for (Int ctr = 0; ctr < nrhs_per_thread; ++ctr) {
    const bool ut = usethread && irhs < nrhs;
    for (Int col = 0; col < n; ++col) {
      if (ut && row > col)
        x(row,irhs) -= x(col,irhs) * L(tdidx,row_minus_base,col);
      mbr.team_barrier();
    }
    irhs += nrhs_batch;
  }
}

// Fuse lunpbwd and ngemm.
template <typename ExeSpace, typename Matrix, typename MVecX>
KIF void lunpbwd_ngemm (const typename ko::TeamPolicy<ExeSpace>::member_type& mbr,
                        const Matrix U, MVecX x, const Int team_size,
                        const Int thread_id_base) {
  Int n, nrhs, nrhs_batch, nrhs_per_thread, row, irhs, tdidx, row_minus_base;
  bool usethread;
  {
    n = static_cast<Int>(U.dimension_1());
    nrhs = static_cast<Int>(x.dimension_1());
    const Int n2 = n << 1;
    const Int n_p2 = calc_stride(n2);
    nrhs_batch = impl::min<Int>(nrhs, team_size / n_p2);
    nrhs_per_thread = (nrhs + nrhs_batch - 1) / nrhs_batch;
    const Int team_rank = mbr.team_rank() - thread_id_base;
    row = team_rank % n_p2;
    irhs = team_rank / n_p2;
    tdidx = row < n ? 2 : 3;
    row_minus_base = row - (row < n ? 0 : n);
    usethread = row < n2 && irhs < nrhs_batch;
  }
  for (Int ctr = 0; ctr < nrhs_per_thread; ++ctr) {
    const bool ut = usethread && irhs < nrhs;
    for (Int col = n-1; col >= 0; --col) {
      const auto n_plus_col = n + col;
      typename Matrix::non_const_value_type xoU = 0;
      if (ut && row <= n_plus_col) {
        xoU = x(n_plus_col,irhs) / U(3,col,col);
        if (row < n_plus_col)
          x(row,irhs) -= xoU * U(tdidx,row_minus_base,col);
      }
      mbr.team_barrier();
      if (ut && row == n_plus_col)
        x(row,irhs) = xoU;
    }
    irhs += nrhs_batch;
  }
}

/* Factorization of a 2x2 set of blocks from a block tridiag matrix. Given
      E = [(al,au) b; c d],
   form
      L = [al 0; L21 L22]
      U = [au U12; 0 U22],
   yielding the identities
      b = al U12
      c = L21 au => au' L21' = c'
      d = L21 U12 + L22 U22.
   Algorithm:
      Factor (EL,EU) <- E in place:
      1a. Solve al b = b.
      1b. Solve au' c' = c'.
      2.  d := d - c*b.
      3.  Factor d.
 */
template <typename Matrix>
KIF void blktridiag_fac2x2_serial (const Matrix a, Matrix b, Matrix c, Matrix d) {
  typedef int Int;
  const Int n = static_cast<Int>(a.extent_int(0));
  lunpfwd_serial(a, b, b);
  lunpbwdtt_serial(a, c, c);
  ngemm(d, c, b);
  lunpfac_serial(d);
}

template <typename ExeSpace, typename Matrix>
KIF void blktridiag_fac2x2 (const typename ko::TeamPolicy<ExeSpace>::member_type& mbr,
                            const Matrix a, Matrix b, Matrix c, Matrix d) {
  typedef int Int;
  const Int n = static_cast<Int>(a.extent_int(0));
  lunpfwd<ExeSpace>(mbr, a, b, b, mbr.team_size(), 0);
  lunpbwdtt<ExeSpace>(mbr, a, c, c);
  ngemm<ExeSpace>(mbr, d, c, b, d);
  lunpfac<ExeSpace>(mbr, d);
}

namespace impl {
template <typename ExeSpace, typename Scalar> struct BtmValuesType {
  typedef typename ko::View<Scalar***, ExeSpace> type;
  enum { is_LayoutStride = false };
};
#ifdef KOKKOS_HAVE_CUDA
template <typename Scalar> struct BtmValuesType<ko::Cuda, Scalar> {
  typedef typename ko::View<Scalar***, ko::LayoutStride, ko::Cuda> type;
  enum { is_LayoutStride = true };
};
#endif
} // namespace impl

template <typename ExeSpace, typename Scalar> struct BtmValuesType {
  typedef typename impl::BtmValuesType<ExeSpace, Scalar>::type type;
  enum { is_LayoutStride = impl::BtmValuesType<ExeSpace, Scalar>::is_LayoutStride };
};
template <typename Scalar> struct BtmValuesType<ko::DefaultHostExecutionSpace, Scalar> {
  typedef typename impl::BtmValuesType<ko::DefaultExecutionSpace, Scalar>::type::HostMirror type;
    enum { is_LayoutStride = impl::BtmValuesType<ko::DefaultExecutionSpace, Scalar>::is_LayoutStride };
};

// A collection of block tridiag matrices.
template <typename ExeSpaceT = ko::DefaultExecutionSpace,
          typename IndexT = Int, typename ScalarT = Scalar>
struct BlockTridiagMatrices {
  typedef ExeSpaceT ExeSpace;
  typedef IndexT Index;
  typedef ScalarT Scalar;
  typedef typename BtmValuesType<ExeSpace, Scalar>::type Blocks;
  typedef typename BtmValuesType<ExeSpace, const Scalar>::type ConstBlocks;
  typedef typename ko::View<Index*, ExeSpace> IndexList;
  typedef typename ko::View<const Index*, ExeSpace> ConstIndexList;

  // Initialize with n tridiags, each having nblkrows, and each block having
  // size bs x bs.
  BlockTridiagMatrices (const Index& ntridiags, const Index& nblkrows, const Index& bs) {
    blocks_ = init_blocks<is_LayoutStride, Blocks>(ntridiags*(3*nblkrows - 2), bs);
    use_blocks_ = blocks_;
    nblkrows_matrix_ = make_structured<IndexList>(ntridiags, nblkrows, tridiag_ptr_, blkrow_);
  }

  Index ntridiags () const { return blkrow_.dimension_0() - 1; }
  Index blocksize () const { return blocks_.dimension_1(); }
  // Number of block rows in the whole matrix.
  Index nblkrows_matrix () const { return nblkrows_matrix_; }

  // Bare getters for optional use in kernel ctors.
  Blocks get_blocks () { return blocks_; }
  ConstBlocks get_blocks () const { return blocks_; }
  IndexList get_tridiag_ptr () { return tridiag_ptr_; }
  ConstIndexList get_tridiag_ptr () const { return tridiag_ptr_; }
  IndexList get_blkrow () { return blkrow_; }
  ConstIndexList get_blkrow () const { return blkrow_; }

  template <typename ES = ExeSpace, typename Ind = Index, typename Scal = Scalar>
  std::shared_ptr<BlockTridiagMatrices<ES, Ind, Scal> > clone () const {
    typedef typename ko::View<Ind*, ES> IndexListOther;
    IndexListOther tridiag_ptr(ko::ViewAllocateWithoutInitializing("BlockTridiagMatrices::tridiag_ptr"),
                                                                   tridiag_ptr_.dimension_0());
    ko::deep_copy(tridiag_ptr, tridiag_ptr_);
    IndexListOther blkrow(ko::ViewAllocateWithoutInitializing("BlockTridiagMatrices::blkrow"),
                                                              blkrow_.dimension_0());
    ko::deep_copy(blkrow, blkrow_);
    typedef typename BtmValuesType<ES, Scal>::type BlocksOther;
    auto blocks = init_blocks<is_LayoutStride, BlocksOther>(blocks_.dimension_0(), blocks_.dimension_1());
    ko::deep_copy(blocks, blocks_);
    auto btm = std::make_shared<BlockTridiagMatrices<ES, Ind, Scal> >(tridiag_ptr, blkrow, blocks, nblkrows_matrix_);
    return btm;
  }

  // For clone's use.
  BlockTridiagMatrices (const IndexList& tridiag_ptr, const IndexList& blkrow, const Blocks& blocks,
                        const Index nblkrows_matrix)
    : tridiag_ptr_(tridiag_ptr), blkrow_(blkrow), blocks_(blocks), use_blocks_(blocks_),
      nblkrows_matrix_(nblkrows_matrix)
  {}

private:
  enum { is_LayoutStride = BtmValuesType<ExeSpace, Scalar>::is_LayoutStride };
  // tridiag_ptr_[i] is the block idx of the start of the i'th tridiag.
  IndexList tridiag_ptr_;
  // blkrow_[i] is the block row at the start of the i'th tridiag.
  IndexList blkrow_;
  // Blocks packed by tridiag.
  Blocks blocks_;
  vdec::Unmanaged<Blocks> use_blocks_;
  Index nblkrows_matrix_;

  // Initialize for a structured mesh.
  template <bool is_LayoutStride, typename BlocksT>
  static BlocksT init_blocks (const Index& nblks, const Index& bs,
                              const typename std::enable_if<is_LayoutStride, int>::type* = 0) {
    // Decent layout for solve, team per tridiag.
    const int order[] = {1, 2, 0};
    const int dims[] = {nblks, bs, bs};
    return BlocksT(ko::ViewAllocateWithoutInitializing("BlockTridiagMatrices::blocks"),
                   ko::LayoutStride::order_dimensions(3, order, dims));
  }
  template <bool is_LayoutStride, typename BlocksT>
  static BlocksT init_blocks (const Index& nblks, const Index& bs,
                              const typename std::enable_if< ! is_LayoutStride, int>::type* = 0) {
    return BlocksT(ko::ViewAllocateWithoutInitializing("BlockTridiagMatrices::blocks"), nblks, bs, bs);
  }
  template <typename IndexListT>
  static Index make_structured (const Index& ntridiags, const Index& nblkrows, IndexListT& tridiag_ptr,
                                IndexListT& blkrow) {
    ko::resize(tridiag_ptr, ntridiags + 1);
    ko::resize(blkrow, ntridiags + 1);
    auto tph = ko::create_mirror_view(tridiag_ptr);
    auto brh = ko::create_mirror_view(blkrow);
    const auto nblks = 3*nblkrows - 2;
    tph[0] = 0;
    brh[0] = 0;
    for (Index i = 0; i < ntridiags; ++i) {
      tph[i+1] = tph[i] + nblks;
      brh[i+1] = brh[i] + nblkrows;
    }
    ko::deep_copy(tridiag_ptr, tph);
    ko::deep_copy(blkrow, brh);
    return brh[ntridiags];
  }
};

//< From SPARC. I wrote this code in SPARC, but it has nothing to do with SPARC
// in particular; it's basic numerical linear algebra serial kernels. The
// motivation is to make this reference impl self-contained w.r.t. performance
// testing.
namespace sparc {
#ifdef MIMIC_SPARC
#pragma message "BUILD WITH MIMIC_SPARC"

typedef unsigned int UInt;
typedef int Index;
typedef double Real;

namespace impl {
// z = y + sign A x. A is row major.
template <int sign>
static inline void gemv (const Real* const A, const UInt n, const Real* x, const Real* y, Real* z,
                         const UInt nrhs=1, const UInt ldx=0) {
  for (UInt irhs = 0; ; ) {
    for (UInt r = 0; r < n; ++r) {
      Real accum = 0;
      const Real* const Ar = A + r*n;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#     pragma ivdep
#endif
      for (UInt c = 0; c < n; ++c)
        accum += Ar[c] * x[c];
      z[r] = y[r] + sign * accum;
    }
    ++irhs;
    if (irhs == nrhs) break;
    x += ldx;
    y += ldx;
    z += ldx;
  }
}

// A -= B C. All are row major.
static inline void ngemm (Real* const A, const UInt n, const Real* const B, const Real* const C) {
  for (UInt i = 0; i < n; ++i) {
    Real* const Ai = A + i*n;
    const Real* const Bi = B + i*n;
    for (UInt k = 0; k < n; ++k) {
      const Real Bik = Bi[k];
      const Real* const Ck = C + k*n;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#     pragma ivdep
#endif
      for (UInt j = 0; j < n; ++j)
        Ai[j] -= Bik * Ck[j];
    }
  }
}

// A = L U with L(i,i) = 1. Danger! No pivoting.
static inline void factor_nopivot (Real* const a, const UInt n) {
  // Outer-product LU.
  for (UInt i = 0; i < n; ++i) {
    Real* const ai = a + i*n;
    Real aii = ai[i];
    aii = 1.0 / aii;
    ai[i] = aii;
    for (UInt r = i+1; r < n; ++r) {
      Real* const ar = a + r*n;
      ar[i] *= aii;
      const Real ari = ar[i];
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#     pragma ivdep
#endif
      for (UInt c = i+1; c < n; ++c)
        ar[c] -= ari*ai[c];
    }
  }
}

// x = L \ b with L(i,i) = 1.
static inline void fwdsub (const Real* const A, const UInt n, const Real* b, Real* x,
                           const UInt nrhs=1, const UInt ldx=0) {
  for (UInt irhs = 0; ; ) {
    x[0] = b[0];
    for (UInt r = 1; r < n; ++r) {
      const Real* const Ar = A + r*n;
      Real a = 0;
      for (UInt c = 0; c < r; ++c)
        a += Ar[c] * x[c];
      x[r] = b[r] - a;
    }
    ++irhs;
    if (irhs == nrhs) break;
    b += ldx;
    x += ldx;
  }
}

// x = U \ b.
static inline void bcksub (const Real* const A, const UInt n, const Real* b, Real* x,
                           const UInt nrhs=1, const UInt ldx=0) {
  for (UInt irhs = 0; ; ) {
    for (int r = static_cast<int>(n)-1; r >= 0; --r) {
      const Real* const Ar = A + r*n;
      Real a = 0;
      for (UInt c = r+1; c < n; ++c)
        a += Ar[c] * x[c];
      x[r] = (b[r] - a) * Ar[r];
    }
    ++irhs;
    if (irhs == nrhs) break;
    b += ldx;
    x += ldx;
  }
}

// In A = [P U; L D], all blocks row major, with P = Lp Up from factor_nopivot,
// set U = Lp \ U and L = L / Up.
static inline void blktridiag_factor_offdiag (const Real* const prev_block, Real* const upper_block,
                                              Real* const lower_block, const UInt n) {
  for (UInt i = 0; i < n; ++i) {
    const Real* const pbi = prev_block + i*n;
    for (UInt j = 0; j < n; ++j) {
      Real* const lbj = lower_block + j*n;
      Real a = 0;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#     pragma ivdep
#endif
      for (UInt k = 0; k < i; ++k)
        a += lbj[k] * prev_block[k*n + i];
      lbj[i] = (lbj[i] - a) * pbi[i];
    }
  }
  for (UInt i = 0; i < n; ++i) {
    const Real* const pbi = prev_block + i*n;
    Real* const ubi = upper_block + i*n;
    for (UInt k = 0; k < i; ++k) {
      const Real pbk = pbi[k];
      const Real* const ubr = upper_block + k*n;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#     pragma ivdep
#endif
      for (UInt j = 0; j < n; ++j) {
        ubi[j] -= pbk * ubr[j];
      }
    }
  }
}

// Template on block size BSZ. If BSZ == 0, do runtime sizing. If BSZ != 0, then
// the block computation is compiled for the specific size BSZ.

template <UInt BSZ, int sign>
static inline void gemv (const Real* const A, const UInt n, const Real* const x, const Real* const y, Real* const z,
                         const UInt nrhs=1, const UInt ldx=0) {
  if (BSZ == 0)
    gemv<sign>(A, n, x, y, z, nrhs, ldx);
  else
    gemv<sign>(A, BSZ, x, y, z, nrhs, ldx);
}

template <UInt BSZ>
static inline void ngemm (Real* const A, const UInt n, const Real* const B, const Real* const C) {
  if (BSZ == 0)
    ngemm(A, n, B, C);
  else
    ngemm(A, BSZ, B, C);
}

template <UInt BSZ>
static inline void factor_nopivot (Real* const a, const UInt n) {
  if (BSZ == 0)
    factor_nopivot(a, n);
  else
    factor_nopivot(a, BSZ);
}

template <UInt BSZ>
static inline void fwdsub (const Real* const A, const UInt n, const Real* const b, Real* const x,
                           const UInt nrhs=1, const UInt ldx=0) {
  if (BSZ == 0)
    fwdsub(A, n, b, x, nrhs, ldx);
  else
    fwdsub(A, BSZ, b, x, nrhs, ldx);
}

template <UInt BSZ>
static inline void bcksub (const Real* const A, const UInt n, const Real* const b, Real* const x,
                           const UInt nrhs=1, const UInt ldx=0) {
  if (BSZ == 0)
    bcksub(A, n, b, x, nrhs, ldx);
  else
    bcksub(A, BSZ, b, x, nrhs, ldx);
}

template <UInt BSZ>
static inline void blktridiag_factor_offdiag (const Real* const prev_block, Real* const upper_block,
                                              Real* const lower_block, const UInt n) {
  if (BSZ == 0)
    blktridiag_factor_offdiag(prev_block, upper_block, lower_block, n);
  else
    blktridiag_factor_offdiag(prev_block, upper_block, lower_block, BSZ);
}

template <UInt BSZ, typename BcrsMatrixT, typename ConstVector, typename Vector>
static void mvp (const BcrsMatrixT& A, const ConstVector& xv, Vector& yv) {
  const auto* const rowptr = A.g().rowptr.data();
  const auto* const colidx = A.g().colidx.data();
  const auto* const values = A.v().data();
  const auto* const x = xv.data();
  auto* const y = yv.data();
  const Index blockrow_end = A.g().rowptr.dimension_0() - 1;
  const Index bsz = A.bs(), bsz2 = square(bsz);
  const Index nrhs = xv.dimension_1(), ldx = xv.dimension_0();
#ifdef _OPENMP
# pragma omp parallel for
#endif
  for (Index blockrow = 0; blockrow < blockrow_end; ++blockrow) {
    auto* const y_lcl = y + bsz*blockrow;
    const Index j_lim = rowptr[blockrow+1];
    bool first = true;
    for (Index j = rowptr[blockrow]; j < j_lim; ++j) {
      const auto* const b = values + j*bsz2;
      const auto blockcol = colidx[j];
      const auto x_start = blockcol*bsz;
      auto* const x_lcl = x + x_start;
      for (Index irhs = 0; irhs < nrhs; ++irhs) {
        const Real* const xk = x_lcl + irhs*ldx;
        Real* const yk = y_lcl + irhs*ldx;
        for (Index lclrow = 0; lclrow < bsz; ++lclrow) {
          Real accum = 0;
          const Real* const br = b + bsz*lclrow;
          for (Index i = 0; i < bsz; ++i)
            accum += br[i] * xk[i];
          if (first) yk[lclrow] = accum;
          else yk[lclrow] += accum;
        }
      }
    }
  }  
}

template <UInt BSZ, typename BlockTridiagMatricesType>
static void factor (BlockTridiagMatricesType& btm) {
  Real* const tds = btm.get_blocks().data();
  const auto* const tridiag_ptr = btm.get_tridiag_ptr().data();
  const Index bsz = btm.blocksize(), bsz2 = square(bsz);
  const Index td_end = btm.get_tridiag_ptr().dimension_0() - 1;
#ifdef _OPENMP
# pragma omp parallel for
#endif
  for (Index td_idx = 0; td_idx < td_end; ++td_idx) {
    Real* const td = tds + bsz2*tridiag_ptr[td_idx];
    factor_nopivot<BSZ>(td, bsz);
    Index b_start = 0;
    const auto nblk = tridiag_ptr[td_idx + 1] - tridiag_ptr[td_idx];
    if (nblk == 1) continue;
    const Index n2x2 = (nblk + 2)/3 - 1;
    for (Index i = 0; i < n2x2; ++i) {
      Real* const B = td + bsz2*(b_start + 2);
      Real* const C = td + bsz2*(b_start + 1);
      blktridiag_factor_offdiag<BSZ>(td + bsz2*b_start, B, C, bsz);
      Real* const D = td + bsz2*(b_start + 3);
      ngemm<BSZ>(D, bsz, C, B);
      factor_nopivot<BSZ>(D, bsz);
      b_start += 3;
    }
  }
}

template <UInt BSZ, typename BlockTridiagMatricesType, typename ConstVector, typename Vector>
static void solve (const BlockTridiagMatricesType& btm, const ConstVector& bv, Vector& xv) {
  const Real* const tds = btm.get_blocks().data();
  const auto* const blkrow = btm.get_blkrow().data();
  const auto* const tridiag_ptr = btm.get_tridiag_ptr().data();
  const auto* const b = bv.data();
  auto* const x = xv.data();
  const Index bsz = btm.blocksize(), bsz2 = square(bsz);
  const Index td_end = btm.get_tridiag_ptr().dimension_0() - 1;
  const Index nrhs = bv.dimension_1(), ldx = bv.dimension_0();
#ifdef _OPENMP
# pragma omp parallel
#endif
  {
#ifdef _OPENMP
#   pragma omp for
#endif
    for (Index td_idx = 0; td_idx < td_end; ++td_idx) {
      const Real* const td = tds + bsz2*tridiag_ptr[td_idx];
      Index br = blkrow[td_idx], b_start = 0;
      fwdsub<BSZ>(td, bsz, b + bsz*br, x + bsz*br, nrhs, ldx);
      const auto nblk = tridiag_ptr[td_idx + 1] - tridiag_ptr[td_idx];
      if (nblk == 1) continue;
      const Index n2x2 = (nblk + 2)/3 - 1;
      ++br;
      for (Index i = 0; i < n2x2; ++i) {
        const auto os = bsz*br;
        const Real* const bp = b + os;
        Real* const xp = x + os;
        gemv<BSZ, -1>(td + (b_start+1)*bsz2, bsz, x + bsz*(br-1), bp, xp, nrhs, ldx);
        fwdsub<BSZ>(td + (b_start+3)*bsz2, bsz, xp, xp, nrhs, ldx);
        b_start += 3;
        ++br;
      }
    }
#ifdef _OPENMP
#   pragma omp for
#endif
    for (Index td_idx = 0; td_idx < td_end; ++td_idx) {
      const Real* const td = tds + bsz2*tridiag_ptr[td_idx];
      const auto nblk = tridiag_ptr[td_idx + 1] - tridiag_ptr[td_idx];
      Index br = blkrow[td_idx+1] - 1, b_start = nblk - 1;
      bcksub<BSZ>(td + bsz2*b_start, bsz, x + bsz*br, x + bsz*br, nrhs, ldx);
      if (nblk == 1) continue;
      const Index n2x2 = (nblk + 2)/3 - 1;
      --br;
      for (Index i = 0; i < n2x2; ++i) {
        Real* const xp = x + bsz*br;
        gemv<BSZ, -1>(td + (b_start-1)*bsz2, bsz, x + bsz*(br+1), xp, xp, nrhs, ldx);
        bcksub<BSZ>(td + (b_start-3)*bsz2, bsz, xp, xp, nrhs, ldx);
        b_start -= 3;
        --br;
      }
    }
  }
}
} // namespace impl
#endif // MIMIC_SPARC

template <typename BcrsMatrixT, typename ConstVector, typename Vector>
static void mvp (const BcrsMatrixT& A, const ConstVector& x, Vector& y) {
#ifdef MIMIC_SPARC
  switch (A.bs()) {
  case 5:  impl::mvp<5> (A, x, y); break;
  case 9:  impl::mvp<9> (A, x, y); break;
  case 15: impl::mvp<15>(A, x, y); break;
  default: impl::mvp<0> (A, x, y); break;
  }
#endif
}

template <typename BlockTridiagMatricesType>
static void factor (BlockTridiagMatricesType& btm) {
#ifdef MIMIC_SPARC
  switch (btm.blocksize()) {
  case 5:  impl::factor<5> (btm); break;
  case 9:  impl::factor<9> (btm); break;
  case 15: impl::factor<15>(btm); break;
  default: impl::factor<0> (btm); break;
  }
#endif
}

template <typename BlockTridiagMatricesType, typename ConstVector, typename Vector>
static void solve (const BlockTridiagMatricesType& btm, const ConstVector& b, Vector& x) {
#ifdef MIMIC_SPARC
  switch (btm.blocksize()) {
  case 5:  impl::solve<5> (btm, b, x); break;
  case 9:  impl::solve<9> (btm, b, x); break;
  case 15: impl::solve<15>(btm, b, x); break;
  default: impl::solve<0> (btm, b, x); break;
  }
#endif
}
} // namespace sparc
//> From SPARC.

template <typename VA, typename VB>
double reldif (VA& a, VB& b) {
  // Bring the vectors to the host. This is just a correctness checker.
  auto x = ko::create_mirror_view(a);
  ko::deep_copy(x, a);
  auto y = ko::create_mirror_view(b);
  ko::deep_copy(y, b);
  double num = 0, den = 0;
  const Int n = static_cast<Int>(x.dimension_0()), nrhs = static_cast<Int>(x.dimension_1());
  for (Int i = 0; i < n; ++i)
    for (Int j = 0; j < nrhs; ++j) {
      num += square(x(i,j) - y(i,j));
      den += square(x(i,j));
    }
  return std::sqrt(num/den);
}

template <typename ExeSpace = ko::DefaultExecutionSpace, typename ...Parms>
struct BcrsMvp {
public:
  typedef ::BcrsMatrix<ExeSpace, Parms...> BcrsMatrix;
  typedef typename BcrsMatrix::Index Index;
  typedef typename BcrsMatrix::Size Size;
  typedef typename BcrsMatrix::Scalar Scalar;
  typedef typename BcrsMatrix::Values Values;
  typedef typename BcrsMatrix::CrsGraphType::Rowptr Rowptr;
  typedef typename BcrsMatrix::CrsGraphType::Colidx Colidx;
  typedef ko::View<Scalar**, Kokkos::LayoutLeft, ExeSpace> Vector;
  typedef ko::View<const Scalar**, Kokkos::LayoutLeft, ExeSpace> ConstVector;
};

template <typename ExeSpace = ko::DefaultExecutionSpace, typename ...Parms>
class BcrsMvpByBlockRow {
public:
  typedef BcrsMvp<ExeSpace, Parms...> Traits;
  typedef typename Traits::BcrsMatrix BcrsMatrixT;
  typedef typename Traits::Index Index;
  typedef typename Traits::Size Size;
  typedef typename Traits::Scalar Scalar;
  typedef typename Traits::Values Values;
  typedef typename Traits::Rowptr Rowptr;
  typedef typename Traits::Colidx Colidx;
  typedef typename Traits::Vector Vector;
  typedef typename Traits::ConstVector ConstVector;

private:
  vdec::Unmanaged<Rowptr> rowptr;
  vdec::Unmanaged<Colidx> colidx;
  vdec::Unmanaged<Values> values;
  vdec::Unmanaged<ConstVector> x;
  vdec::Unmanaged<Vector> y;

public:
  BcrsMvpByBlockRow (const BcrsMatrixT& A)
    : rowptr(A.g().rowptr), colidx(A.g().colidx), values(A.v())
  {}

  KIF void operator() (const Index& blockrow) const {
    const Index bsz = values.dimension_1();
    const Index nrhs = x.dimension_1();
    const auto y_start = blockrow*bsz;
    const auto y_lcl = ko::subview(y, ko::make_pair(y_start, y_start + bsz), ko::ALL());
    for (Index irhs = 0; irhs < nrhs; ++irhs)
      for (Index i = 0; i < bsz; ++i)
        y_lcl(i, irhs) = 0;
    const Index j_lim = rowptr(blockrow+1);
    for (Index j = rowptr(blockrow); j < j_lim; ++j) {
      const Index blockcol = colidx(j);
      const auto b = ko::subview(values, j, ko::ALL(), ko::ALL());
      const Index x_start = blockcol*bsz;
      const auto x_lcl = ko::subview(x, ko::make_pair(x_start, x_start + bsz), ko::ALL());
      for (Index lclrow = 0; lclrow < bsz; ++lclrow) {
        for (Index irhs = 0; irhs < nrhs; ++irhs) {
          Scalar accum = 0;
          for (Index i = 0; i < bsz; ++i)
            accum += b(lclrow, i) * x_lcl(i, irhs);
          y_lcl(lclrow, irhs) += accum;
        }
      }
    }
  }

  void run (const ConstVector& x_, const Vector& y_) {
    x = x_;
    y = y_;
    ko::parallel_for(rowptr.dimension_0() - 1, *this);
  }
};

template <typename ExeSpace = ko::DefaultExecutionSpace, typename ...Parms>
class BcrsMvpByRow {
public:
  typedef BcrsMvp<ExeSpace, Parms...> Traits;
  typedef typename Traits::BcrsMatrix BcrsMatrixT;
  typedef typename Traits::Index Index;
  typedef typename Traits::Size Size;
  typedef typename Traits::Scalar Scalar;
  typedef typename Traits::Values Values;
  typedef typename Traits::Rowptr Rowptr;
  typedef typename Traits::Colidx Colidx;
  typedef typename Traits::Vector Vector;
  typedef typename Traits::ConstVector ConstVector;

private:
  vdec::Unmanaged<Rowptr> rowptr;
  vdec::Unmanaged<Colidx> colidx;
  vdec::Unmanaged<Values> values;
  vdec::Unmanaged<ConstVector> x;
  vdec::Unmanaged<Vector> y;

public:
  BcrsMvpByRow (const BcrsMatrixT& A)
    : rowptr(A.g().rowptr), colidx(A.g().colidx), values(A.v())
  {}

  // A thread maps to a point row of the matrix.
  KIF void operator() (const Index& k) const {
    const Index irhs = k / y.dimension_0();
    const Index row = k % y.dimension_0();
    const Index bsz = values.dimension_1();
    const Index blockrow = row / bsz;
    const Index lclrow = row % bsz;
    Scalar accum = 0;
    const Index j_lim = rowptr(blockrow+1);
    for (Index j = rowptr(blockrow); j < j_lim; ++j) {
      const Index blockcol = colidx(j);
      const auto b = ko::subview(values, j, ko::ALL(), ko::ALL());
      const Index x_start = blockcol*bsz;
      const auto x_lcl = ko::subview(x, ko::make_pair(x_start, x_start + bsz), irhs);
      for (Index i = 0; i < bsz; ++i)
        accum += b(lclrow, i) * x_lcl(i);
    }
    y(row, irhs) = accum;
  }

  void run (const ConstVector& x_, const Vector& y_) {
    x = x_;
    y = y_;
    ko::parallel_for(x.size(), *this);
  }
};

struct MvpImpl { enum Enum { by_block_row, by_row, sparc }; };

template <typename ExeSpace, typename ...Parms>
void mvp (const BcrsMatrix<ExeSpace, Parms...>& A,
          const typename BcrsMvp<ExeSpace, Parms...>::ConstVector& x,
          const typename BcrsMvp<ExeSpace, Parms...>::Vector& y,
          const MvpImpl::Enum impl = (impl::is_cuda<ExeSpace>::value ? MvpImpl::by_row :
                                      impl::mimic_sparc ? MvpImpl::sparc :
                                      MvpImpl::by_block_row)) {
  switch (impl) {
  case MvpImpl::by_block_row:
    BcrsMvpByBlockRow<ExeSpace, Parms...>(A).run(x, y);
    break;
  case MvpImpl::sparc:
    sparc::mvp(A, x, y);
    break;
  case MvpImpl::by_row:
  default:
    BcrsMvpByRow<ExeSpace, Parms...>(A).run(x, y);
  }
}

// Extract the k-line block tridiag matrices from a BcrsMatrix and put the data
// into a BlockTridiagMatrices. Tridiags are stored in block-column-major order,
// not (as one might think) block-row-major. That is so that we can fuse the
// trisolve and gemv in the solve phase.
template <typename BcrsMatrixT, typename BlockTridiagMatricesT,
          typename ExeSpace = ko::DefaultExecutionSpace>
class ExtractBlockTridiag {
  const StructuredBlock sb;
  const Int nlines;
  vdec::ConstUnmanaged<typename BcrsMatrixT::CrsGraphType::Rowptr> rowptr;
  vdec::ConstUnmanaged<typename BcrsMatrixT::CrsGraphType::Colidx> colidx;
  vdec::ConstUnmanaged<typename BcrsMatrixT::Values> v;
  std::shared_ptr<BlockTridiagMatricesT> btm;
  vdec::Unmanaged<typename BlockTridiagMatricesT::ConstIndexList> td_ptr;
  vdec::Unmanaged<typename BlockTridiagMatricesT::Blocks> td_blocks;

  KIF Size copy (const Int& tridiag_idx, const Int& blocks_idx, Size v_idx, const Size& id,
                 const int& direction) const {
    const auto d = ko::subview(td_blocks, td_ptr[tridiag_idx] + blocks_idx, ko::ALL(), ko::ALL());
    const Size idd = id + direction;
    while (colidx(v_idx) != idd) ++v_idx;
    const auto s = ko::subview(v, v_idx, ko::ALL(), ko::ALL());
    for (int i = 0; i < s.extent_int(0); ++i)
      for (int j = 0; j < s.extent_int(1); ++j)
        d(i,j) = s(i,j);
    return v_idx;    
  }

public:
  struct ByLine_tag {};
  struct ByBlock_tag {};

  ExtractBlockTridiag (const StructuredBlock& sb_, const BcrsMatrixT& mat)
    : sb(sb_), nlines(sb.ni*sb.nj), rowptr(mat.g().rowptr), colidx(mat.g().colidx), v(mat.v())
  {
    btm = std::make_shared<BlockTridiagMatricesT>(nlines, sb.nk, mat.bs());
    td_ptr = btm->get_tridiag_ptr();
    td_blocks = btm->get_blocks();
  }

  KIF void operator() (const ByLine_tag&, const Int tridiag_idx) const {
    const Int i = tridiag_idx / sb.nj, j = tridiag_idx % sb.nj;
    typename BlockTridiagMatricesT::Index blocks_idx = 0;
    const auto id = sb.ijk2id(i,j,0);
    auto v_idx = copy(tridiag_idx, blocks_idx, rowptr(id), id, 0);
    if (sb.nk < 2) return;
    copy(tridiag_idx, blocks_idx+2, v_idx+1, id, 1);
    blocks_idx += 3;
    if (sb.nk > 2)
      for (Int k = 1; k < sb.nk - 1; ++k) {
        const auto id = sb.ijk2id(i,j,k);
        v_idx = copy(tridiag_idx, blocks_idx-2, rowptr(id), id, -1);
        v_idx = copy(tridiag_idx, blocks_idx, v_idx+1, id,  0);
        copy(tridiag_idx, blocks_idx+2, v_idx+1, id,  1);
        blocks_idx += 3;
      }
    {
      const auto id = sb.ijk2id(i,j,sb.nk-1);
      v_idx = copy(tridiag_idx, blocks_idx-2, rowptr(id), id, -1);
      copy(tridiag_idx, blocks_idx, v_idx, id,  0);
    }
  }

  KIF void operator() (const ByBlock_tag&, const Int& tid) const {
    const Int n = td_blocks.dimension_1(), n2 = n*n;
    // Indices into tridiags.
    Int block_idx = tid / n2;
    const auto nblks_per_tridiag = td_ptr[1];
    const Int tridiag_idx = block_idx / nblks_per_tridiag;
    const Int tridiag_block_idx = block_idx % nblks_per_tridiag;
    const Int i = tridiag_idx / sb.nj;
    const Int j = tridiag_idx % sb.nj;
    const Int k = (tridiag_block_idx + 1) / 3;
    // Indices into v.
    const int direction = static_cast<int>((tridiag_block_idx + 1) % 3) - 1;
    const auto id = sb.ijk2id(i,j,k);
    const auto idd = id + direction;
    Int v_idx = rowptr(id);
    while (colidx(v_idx) != idd) ++v_idx;
    // Assign row and col values for this thread to handle.
    const Int bi = tid % n2;
    const Int row = bi % n;
    const Int col = bi / n;
    // Now we're ready to copy. The +direction is because the tridiag is stored
    // block-column-major. Hence we want to swap the off-diag blocks, which is
    // accomplished by adding +/-1 to them. (This does +0 if an on-diag block.)
    td_blocks(block_idx + direction, row, col) = v(v_idx, row, col);
  }

  // This can be called multiple times. It makes sense to do so if (i) the
  // values in the BCRS matrix change and (ii) you want the same
  // BlockTridiagMatrices structure to hold the new values.
  std::shared_ptr<BlockTridiagMatricesT> run () {
    if (true || impl::is_cuda<ExeSpace>::value) {
      auto policy = ko::RangePolicy<ExeSpace, ByBlock_tag>(0, td_blocks.dimension_0()*square(td_blocks.dimension_1()));
      ko::parallel_for(policy, *this);
    } else {
      ko::parallel_for(ko::RangePolicy<ExeSpace, ByLine_tag>(0, nlines), *this);
    }
    return btm;
  }
};

// Convenience function.
template <typename BcrsMatrixT>
std::shared_ptr<BlockTridiagMatrices<typename BcrsMatrixT::ExeSpace, typename BcrsMatrixT::Index,
                                     typename BcrsMatrixT::Scalar> >
extract_block_tridiag (const StructuredBlock& sb, const BcrsMatrixT& mat) {
  typedef BlockTridiagMatrices<typename BcrsMatrixT::ExeSpace, typename BcrsMatrixT::Index, 
                               typename BcrsMatrixT::Scalar> BlockTridiagMatricesT;
  return ExtractBlockTridiag<BcrsMatrixT, BlockTridiagMatricesT>(sb, mat).run();
}

// Check that the extraction is correct.
template <typename BcrsMatrixT, typename BlockTridiagMatricesT>
bool check_block_tridiag (const StructuredBlock& sb, const BcrsMatrixT& mat,
                          BlockTridiagMatricesT& btm) {
  // Check this on the host in serial.
  auto rowptr = ko::create_mirror_view(mat.g().rowptr);
  ko::deep_copy(rowptr, mat.g().rowptr);
  auto colidx = ko::create_mirror_view(mat.g().colidx);
  ko::deep_copy(colidx, mat.g().colidx);
  auto values = ko::create_mirror_view(mat.v());
  ko::deep_copy(values, mat.v());
  auto blocks = ko::create_mirror_view(btm.get_blocks());
  ko::deep_copy(blocks, btm.get_blocks());
  auto tridiag_ptr = ko::create_mirror_view(btm.get_tridiag_ptr());
  ko::deep_copy(tridiag_ptr, btm.get_tridiag_ptr());
  const auto nblks = tridiag_ptr[1];
  auto is_same = [&] (const Size& tridiag_idx, const Size& b_idx, const Size& v_idx) -> bool {
    auto b = ko::subview(blocks, nblks*tridiag_idx + b_idx, ko::ALL(), ko::ALL());
    auto v = ko::subview(values, v_idx, ko::ALL(), ko::ALL());
    for (int i = 0; i < v.extent_int(0); ++i)
      for (int j = 0; j < v.extent_int(1); ++j)
        if (b(i,j) != v(i,j)) return false;
    return true;
  };
  // Very straightforward loops for clarity, since this is a correctness
  // check. There is some complexity, expressed in the b_idx de/increments,
  // because the tridiag is stored block-col-major but the matrix is
  // block-row-major (CRS).
  bool same = true;
  for (Int i = 0, td_idx = 0; i < sb.ni; ++i)
    for (Int j = 0; j < sb.nj; ++j, ++td_idx) {
      Size b_idx = 0;
      for (Int k = 0; k < sb.nk; ++k) {
        const auto row = sb.ijk2id(i,j,k);
        Size v_idx = rowptr(row);
        const Size v_idx_lim = rowptr(row + 1);
        if (k == 0) {
          for ( ; v_idx < v_idx_lim; ++v_idx) {
            const auto col = colidx(v_idx);
            if (col == row || col == row+1) {
              same = is_same(td_idx, b_idx, v_idx) && same;
              b_idx += 2;
            }
          }
          b_idx = 1;
        } else if (k+1 == sb.nk) {
          for ( ; v_idx < v_idx_lim; ++v_idx) {
            const auto col = colidx(v_idx);
            if (col == row-1 || col == row) {
              same = is_same(td_idx, b_idx, v_idx) && same;
              b_idx += 2;
            }
          }
        } else {
          for ( ; v_idx < v_idx_lim; ++v_idx) {
            const auto col = colidx(v_idx);
            if (col == row-1 || col == row || col == row+1) {
              same = is_same(td_idx, b_idx, v_idx) && same;
              b_idx += 2;
            }
          }
          b_idx -= 3;
        }
      }
    }
  return same;
}

template <typename ExeSpace, typename VectorT>
class FillVector {
  typedef VectorT Vector;
  vdec::Unmanaged<Vector> v;
public:
  FillVector (const Vector& v_) : v(v_) {}
  KIF void operator() (const Int i) const {
    for (int j = 0; j < v.extent_int(1); ++j)
      v(i,j) = static_cast<double>((i+j) % 7) - 3;
  }
  void run () { ko::parallel_for(v.dimension_0(), *this); }
};

template <typename ExeSpace, typename ...Parms>
void fill_vector (const ko::View<ExeSpace, Parms...>& v) {
  FillVector<ExeSpace, ko::View<ExeSpace, Parms...> >(v).run();
}

// Factor the block tridiag matrices. In this impl, each matrix is factored in
// serial.
template <typename ExeSpace, typename ...Parms>
class FactorSerialByLine {
  typedef BlockTridiagMatrices<ExeSpace, Parms...> BlockTridiagMatricesT;
  vdec::Unmanaged<typename BlockTridiagMatricesT::ConstIndexList> tridiag_ptr;
  vdec::Unmanaged<typename BlockTridiagMatricesT::Blocks> blocks;

public:
  FactorSerialByLine (BlockTridiagMatricesT& btm)
    : tridiag_ptr(btm.get_tridiag_ptr()), blocks(btm.get_blocks())
  {}

  KIF void operator() (const Int td_idx) const {
    typedef typename BlockTridiagMatricesT::Index Index;
    const auto td = ko::subview(blocks, ko::make_pair(tridiag_ptr[td_idx], tridiag_ptr[td_idx+1]),
                                ko::ALL(), ko::ALL());
    Index b_start = 0;
    lunpfac_serial(ko::subview(td, b_start, ko::ALL(), ko::ALL()));
    const auto nblk = td.dimension_0();
    if (nblk == 1) return;
    const Index n2x2 = (nblk + 2)/3 - 1;
    for (Index i = 0; i < n2x2; ++i) {
      blktridiag_fac2x2_serial(ko::subview(td, b_start,     ko::ALL(), ko::ALL()),
                               ko::subview(td, b_start + 2, ko::ALL(), ko::ALL()),
                               ko::subview(td, b_start + 1, ko::ALL(), ko::ALL()),
                               ko::subview(td, b_start + 3, ko::ALL(), ko::ALL()));
      b_start += 3;
    }
  }

  void factor () { ko::parallel_for(tridiag_ptr.dimension_0() - 1, *this); }
};

template <typename ExeSpace, typename ...Parms>
class FactorTeamByLine {
  typedef BlockTridiagMatrices<ExeSpace, Parms...> BlockTridiagMatricesT;
  typedef typename ko::TeamPolicy<ExeSpace>::member_type Member;

  vdec::Unmanaged<typename BlockTridiagMatricesT::ConstIndexList> tridiag_ptr;
  vdec::Unmanaged<typename BlockTridiagMatricesT::Blocks> blocks;
  Int team_size;

public:
  FactorTeamByLine (BlockTridiagMatricesT& btm)
    : tridiag_ptr(btm.get_tridiag_ptr()), blocks(btm.get_blocks())
  {
    if (static_cast<int>(blocks.dimension_1()) > max_team_size) {
      std::stringstream ss;
      ss << "GPU FactorTeamByLine can't handle block size > " << max_team_size;
      throw std::runtime_error(ss.str());
    }
  }

  KIF void operator() (const Member& mbr) const {    
    typedef typename BlockTridiagMatricesT::Index Index;
    const auto td_idx = mbr.league_rank();
    const auto td = ko::subview(blocks, ko::make_pair(tridiag_ptr[td_idx], tridiag_ptr[td_idx+1]),
                                ko::ALL(), ko::ALL());
    Index b_start = 0;
    lunpfac<ExeSpace>(mbr, ko::subview(td, b_start, ko::ALL(), ko::ALL()));
    const auto nblk = tridiag_ptr[td_idx+1] - tridiag_ptr[td_idx];
    if (nblk == 1) return;
    const Index n2x2 = (nblk + 2)/3 - 1;
    for (Index i = 0; i < n2x2; ++i) {
      blktridiag_fac2x2<ExeSpace>(mbr,
                                  ko::subview(td, b_start,     ko::ALL(), ko::ALL()),
                                  ko::subview(td, b_start + 2, ko::ALL(), ko::ALL()),
                                  ko::subview(td, b_start + 1, ko::ALL(), ko::ALL()),
                                  ko::subview(td, b_start + 3, ko::ALL(), ko::ALL()));
      b_start += 3;
    }
  }

  void factor () {
    team_size = calc_team_size(square(blocks.dimension_1()));
    ko::parallel_for(ko::TeamPolicy<ExeSpace>(tridiag_ptr.dimension_0() - 1, team_size), *this);
  }
};

struct TridiagImpl { enum Enum { serial, team, sparc }; };

// Convenience function.
template <typename ExeSpace, typename ...Parms>
void factor (BlockTridiagMatrices<ExeSpace, Parms...>& btm,
             const TridiagImpl::Enum impl = (impl::is_cuda<ExeSpace>::value ? TridiagImpl::team :
                                            impl::mimic_sparc ? TridiagImpl::sparc :
                                            TridiagImpl::serial)) {
  switch (impl) {
  case TridiagImpl::team:
    FactorTeamByLine<ExeSpace, Parms...>(btm).factor();
    break;
  case TridiagImpl::sparc:
    sparc::factor(btm);
    break;
  case TridiagImpl::serial:
  default:
    FactorSerialByLine<ExeSpace, Parms...>(btm).factor();
  }
}

// Solve A x = b, where A is given by a BlockTridiagMatrices.
template <typename ExeSpace = ko::DefaultExecutionSpace, typename ...Parms>
struct Solve {
  typedef ::BlockTridiagMatrices<ExeSpace, Parms...> BlockTridiagMatrices;
  typedef typename BlockTridiagMatrices::Scalar Scalar;
  typedef typename BlockTridiagMatrices::Index Index;
  typedef ko::View<Scalar**, Kokkos::LayoutLeft, ExeSpace> Vector;
  typedef ko::View<const Scalar**, Kokkos::LayoutLeft, ExeSpace> ConstVector;
};

template <typename ExeSpace = ko::DefaultExecutionSpace, typename ...Parms>
class CheckSolve {
  typedef Solve<ExeSpace, Parms...> Traits;
  typedef typename Traits::BlockTridiagMatrices BlockTridiagMatricesT;

public:
  typedef typename Traits::Scalar Scalar;
  typedef typename Traits::Index Index;
  typedef typename Traits::Vector Vector;
  typedef typename Traits::ConstVector ConstVector;

private:
  vdec::Unmanaged<typename BlockTridiagMatricesT::ConstIndexList> tridiag_ptr, blkrow;
  vdec::Unmanaged<typename BlockTridiagMatricesT::ConstBlocks> blocks;
  vdec::Unmanaged<ConstVector> b;
  vdec::Unmanaged<Vector> x;

public:
  CheckSolve (const BlockTridiagMatricesT& btm)
    : tridiag_ptr(btm.get_tridiag_ptr()), blkrow(btm.get_blkrow()), blocks(btm.get_blocks())
  {}

  KIF void operator() (const Int td_idx) const {
    const auto td = ko::subview(blocks, ko::make_pair(tridiag_ptr[td_idx], tridiag_ptr[td_idx+1]),
                                ko::ALL(), ko::ALL());
    const Index bsz = static_cast<Index>(td.dimension_1());
    const Index nblk = static_cast<Index>(td.dimension_0());
    const Index nrhs = static_cast<Index>(b.dimension_1());
    const Index r0 = bsz*blkrow[td_idx];
    const Index n2x2 = static_cast<Index>((nblk + 2)/3 - 1);
    for (Index bi = 0, tdidx = 0; bi < n2x2; ++bi) {
      for (Index cnt = bi == 0 ? 0 : 1; cnt < 4; ++cnt, ++tdidx) {
        const auto br = bi + (cnt % 2);
        const auto bc = bi + (cnt / 2);
        const auto blk = ko::subview(td, tdidx, ko::ALL(), ko::ALL());
        const auto b_bc = ko::subview(b, ko::make_pair(r0 + bsz*bc, r0 + bsz*(bc + 1)), ko::ALL());
        const auto x_br = ko::subview(x, ko::make_pair(r0 + bsz*br, r0 + bsz*(br + 1)), ko::ALL());
        for (Index r = 0; r < bsz; ++r)
          for (Index i = 0; i < nrhs; ++i) {
            Scalar accum = 0;
            for (Index c = 0; c < bsz; ++c)
              accum += blk(r,c) * b_bc(c,i);
            x_br(r,i) += accum;
          }
      }
    }
  }

  bool check (const ConstVector& b_, const Vector& x_, const double tol) {
    b = x_;
    const auto x_mgd = Vector("mvp", b.dimension_0(), b.dimension_1());
    x = x_mgd;
    ko::parallel_for(tridiag_ptr.dimension_0() - 1, *this);
    const bool pass = reldif(x, b_) <= tol;
    return pass;
  }
};

template <typename ExeSpace, typename ...Parms>
bool check_solve (const BlockTridiagMatrices<ExeSpace, Parms...>& btm,
                  const typename Solve<ExeSpace, Parms...>::ConstVector& b,
                  const typename Solve<ExeSpace, Parms...>::Vector& x,
                  const double tol) {
  return CheckSolve<ExeSpace, Parms...>(btm).check(b, x, tol);
}

// This impl is decent, although not great, on CPU/KNL, but it is
// catastrophically bad on GPU. The obvious thing to impl next for ko::Cuda is a
// team per tridiag. On CPU/KNL, one could get nearer-optimal performance by
// templating operations on block size, but teams are probably unnecessary and
// possibly counter-productive.
template <typename ExeSpace = ko::DefaultExecutionSpace, typename ...Parms>
class SolveSerialByLine {
  typedef Solve<ExeSpace, Parms...> Traits;
  typedef typename Traits::BlockTridiagMatrices BlockTridiagMatricesT;

public:
  typedef typename Traits::Scalar Scalar;
  typedef typename Traits::Index Index;
  typedef typename Traits::Vector Vector;
  typedef typename Traits::ConstVector ConstVector;

private:
  vdec::Unmanaged<typename BlockTridiagMatricesT::ConstIndexList> tridiag_ptr, blkrow;
  vdec::Unmanaged<typename BlockTridiagMatricesT::ConstBlocks> blocks;
  vdec::Unmanaged<ConstVector> b;
  vdec::Unmanaged<Vector> x;

public:
  struct L_tag {};
  struct U_tag {};

  SolveSerialByLine (const BlockTridiagMatricesT& btm)
    : tridiag_ptr(btm.get_tridiag_ptr()), blkrow(btm.get_blkrow()), blocks(btm.get_blocks())
  {}

  KIF void operator() (const L_tag&, const Int td_idx) const {
    const auto td = ko::subview(blocks, ko::make_pair(tridiag_ptr[td_idx], tridiag_ptr[td_idx+1]),
                                ko::ALL(), ko::ALL());
    Index b_start = 0, br = blkrow[td_idx];
    const auto bsz = td.dimension_1();
    lunpfwd_serial(ko::subview(td, b_start, ko::ALL(), ko::ALL()),
                   ko::subview(b, ko::make_pair(bsz*br, bsz*(br + 1)), ko::ALL()),
                   ko::subview(x, ko::make_pair(bsz*br, bsz*(br + 1)), ko::ALL()));
    const Index nblk = td.dimension_0();
    if (nblk == 1) return;
    const Index n2x2 = static_cast<Index>((nblk + 2)/3 - 1);
    ++br;
    for (Index i = 0; i < n2x2; ++i) {
      const auto b_br = ko::subview(b, ko::make_pair(bsz*br, bsz*(br + 1)), ko::ALL());
      const auto x_br = ko::subview(x, ko::make_pair(bsz*br, bsz*(br + 1)), ko::ALL());
      ngemm(b_br, ko::subview(td, b_start + 1, ko::ALL(), ko::ALL()),
            ko::subview(x, ko::make_pair(bsz*(br-1), bsz*br), ko::ALL()), x_br);
      lunpfwd_serial(ko::subview(td, b_start + 3, ko::ALL(), ko::ALL()), x_br, x_br);
      b_start += 3;
      ++br;
    }
  }

  KIF void operator() (const U_tag&, const Int td_idx) const {
    const auto td = ko::subview(blocks, ko::make_pair(tridiag_ptr[td_idx], tridiag_ptr[td_idx+1]),
                                ko::ALL(), ko::ALL());
    Index b_start = td.dimension_0() - 1, br = blkrow[td_idx+1] - 1;
    const auto bsz = td.dimension_1();
    lunpbwd_serial(ko::subview(td, b_start, ko::ALL(), ko::ALL()),
                   ko::subview(b, ko::make_pair(bsz*br, bsz*(br + 1)), ko::ALL()),
                   ko::subview(x, ko::make_pair(bsz*br, bsz*(br + 1)), ko::ALL()));
    const Index nblk = static_cast<Index>(td.dimension_0());
    if (nblk == 1) return;
    const Index n2x2 = (nblk + 2)/3 - 1;
    --br;
    for (Index i = 0; i < n2x2; ++i) {
      const auto b_br = ko::subview(b, ko::make_pair(bsz*br, bsz*(br + 1)), ko::ALL());
      const auto x_br = ko::subview(x, ko::make_pair(bsz*br, bsz*(br + 1)), ko::ALL());
      ngemm(b_br, ko::subview(td, b_start - 1, ko::ALL(), ko::ALL()),
            ko::subview(x, ko::make_pair(bsz*(br + 1), bsz*(br + 2)), ko::ALL()), x_br);
      lunpbwd_serial(ko::subview(td, b_start - 3, ko::ALL(), ko::ALL()), x_br, x_br);
      b_start -= 3;
      --br;
    }
  }

  void solve (const ConstVector& b_, const Vector& x_) {
    // Solve L x = b.
    b = b_;
    x = x_;
    ko::parallel_for(ko::RangePolicy<ExeSpace, L_tag>(0, tridiag_ptr.dimension_0() - 1), *this);
    ko::fence();
    // Solve U x = x.
    b = x;
    ko::parallel_for(ko::RangePolicy<ExeSpace, U_tag>(0, tridiag_ptr.dimension_0() - 1), *this);
  }
};

// "Fused" refers to the fact that the trisolve and the ngemm are fused so that
// more threads in a team can be used. (This is motivated by the fact that to
// get 100% occupancy, team size must be >= 128.)
template <typename ExeSpace = ko::DefaultExecutionSpace, typename ...Parms>
class SolveFusedTeamByLine {
  typedef Solve<ExeSpace, Parms...> Traits;
  typedef typename Traits::BlockTridiagMatrices BlockTridiagMatricesT;
  typedef typename ko::TeamPolicy<ExeSpace>::member_type Member;

public:
  typedef typename Traits::Scalar Scalar;
  typedef typename Traits::Index Index;
  typedef typename Traits::Vector Vector;
  typedef typename Traits::ConstVector ConstVector;

private:
  vdec::Unmanaged<typename BlockTridiagMatricesT::ConstIndexList> tridiag_ptr, blkrow;
  vdec::Unmanaged<typename BlockTridiagMatricesT::ConstBlocks> blocks;
  vdec::Unmanaged<ConstVector> b;
  vdec::Unmanaged<Vector> x;
  Int ko_team_size, team_size, tridiags_per_team;

  KIF void calc_thread_indices (const Member& mbr, Index& team_tridiag_idx, Index& td_idx,
                                Index& thread_id_base) const {
    team_tridiag_idx = mbr.team_rank() / team_size;
    td_idx = tridiags_per_team * mbr.league_rank() + team_tridiag_idx;
    thread_id_base = team_tridiag_idx * team_size;
  }

public:
  struct L_tag {};
  struct U_tag {};

  SolveFusedTeamByLine (const BlockTridiagMatricesT& btm)
    : tridiag_ptr(btm.get_tridiag_ptr()), blkrow(btm.get_blkrow()), blocks(btm.get_blocks())
  {
    if (blocks.dimension_1() > max_team_size/2) {
      std::stringstream ss;
      ss << "GPU SolveFusedTeamByLine can't handle block size > " << max_team_size/2;
      throw std::runtime_error(ss.str());
    }
  }

  KIF void operator() (const L_tag&, const Member& mbr) const {
    Index team_tridiag_idx, td_idx, thread_id_base;
    calc_thread_indices(mbr, team_tridiag_idx, td_idx, thread_id_base);
    //todo Determine whether this is truly allowed. It happens to pass the
    // tests, but maybe in general the fact that the __syncthreads() won't be
    // done in all threads is bad.
    if (td_idx >= static_cast<Index>(blkrow.dimension_0() - 1) ||
        team_tridiag_idx >= tridiags_per_team)
      return;
    const auto td = ko::subview(blocks, ko::make_pair(tridiag_ptr[td_idx], tridiag_ptr[td_idx+1]),
                                ko::ALL(), ko::ALL());
    Index b_start = 0, br = blkrow[td_idx];
    const auto bsz = td.dimension_1();
    const Index nblk = td.dimension_0();
    if (nblk == 1) {
      const auto x_br = ko::subview(x, ko::make_pair(bsz*br, bsz*(br + 1)), ko::ALL());
      lunpfwd<ExeSpace>(mbr, ko::subview(td, b_start, ko::ALL(), ko::ALL()), x_br, x_br,
                        team_size, thread_id_base);
      return;
    }
    const Index n2x2 = static_cast<Index>((nblk + 2)/3 - 1);
    for (Index i = 0; i < n2x2; ++i) {
      lunpfwd_ngemm<ExeSpace>(mbr,
                              ko::subview(td, ko::make_pair(b_start, b_start + 3), ko::ALL(), ko::ALL()),
                              ko::subview(x, ko::make_pair(bsz*br, bsz*(br + 2)), ko::ALL()),
                              team_size, thread_id_base);
      b_start += 3;
      ++br;
    }
    {
      const auto x_br = ko::subview(x, ko::make_pair(bsz*br, bsz*(br + 1)), ko::ALL());
      lunpfwd<ExeSpace>(mbr, ko::subview(td, b_start, ko::ALL(), ko::ALL()), x_br, x_br,
                        team_size, thread_id_base);
    }
  }

  KIF void operator() (const U_tag&, const Member& mbr) const {
    Index team_tridiag_idx, td_idx, thread_id_base;
    calc_thread_indices(mbr, team_tridiag_idx, td_idx, thread_id_base);
    if (td_idx >= static_cast<Index>(blkrow.dimension_0() - 1) ||
        team_tridiag_idx >= tridiags_per_team)
      return;
    const auto td = ko::subview(blocks, ko::make_pair(tridiag_ptr[td_idx], tridiag_ptr[td_idx+1]),
                                ko::ALL(), ko::ALL());
    Index b_start = td.dimension_0() - 1, br = blkrow[td_idx+1] - 1;
    const auto bsz = td.dimension_1();
    const Index nblk = static_cast<Index>(td.dimension_0());
    if (nblk == 1) {
      lunpbwd<ExeSpace>(mbr, ko::subview(td, b_start, ko::ALL(), ko::ALL()),
                        ko::subview(b, ko::make_pair(bsz*br, bsz*(br + 1)), ko::ALL()),
                        ko::subview(x, ko::make_pair(bsz*br, bsz*(br + 1)), ko::ALL()),
                        team_size, thread_id_base);
      return;
    }
    const Index n2x2 = (nblk + 2)/3 - 1;
    for (Index i = 0; i < n2x2; ++i) {
      lunpbwd_ngemm<ExeSpace>(mbr,
                              ko::subview(td, ko::make_pair(b_start - 3, b_start + 1), ko::ALL(), ko::ALL()),
                              ko::subview(x, ko::make_pair(bsz*(br - 1), bsz*(br + 1)), ko::ALL()),
                              team_size, thread_id_base);
      b_start -= 3;
      --br;
    }
    {
      const auto x_br = ko::subview(x, ko::make_pair(bsz*br, bsz*(br + 1)), ko::ALL());
      lunpbwd<ExeSpace>(mbr, ko::subview(td, b_start, ko::ALL(), ko::ALL()), x_br, x_br,
                        team_size, thread_id_base);
    }
  }

  std::size_t team_shmem_size (const int& team_size) const { return 0; }

  void solve (const ConstVector& b_, const Vector& x_) {
    static const Int desired_team_size = 128;
    tridiags_per_team = 1;
    const auto bsz = blocks.dimension_1();
    ko_team_size = calc_team_size(2 * bsz * b_.dimension_1(), 1);
    if (ko_team_size > desired_team_size && bsz <= desired_team_size >> 1)
      ko_team_size = desired_team_size;
#if 1 // Set to 0 to get 1 tridiag per team.
    if (ko_team_size <= desired_team_size >> 1) {
      tridiags_per_team = desired_team_size / ko_team_size;
      ko_team_size = desired_team_size;
    }
#endif
    team_size = ko_team_size / tridiags_per_team;
    const auto ntridiags = blkrow.dimension_0() - 1;
    const Int league_size = (ntridiags + tridiags_per_team - 1) / tridiags_per_team;
    // Solve L x = b.
    x = x_;
    ko::deep_copy(x, b_);
    ko::parallel_for(ko::TeamPolicy<ExeSpace, L_tag>(league_size, ko_team_size), *this);
    ko::fence();
    // Solve U x = x.
    b = x;
    ko::parallel_for(ko::TeamPolicy<ExeSpace, U_tag>(league_size, ko_team_size), *this);
  }
};

// Convenience function.
template <typename ExeSpace, typename ...Parms>
void solve (const BlockTridiagMatrices<ExeSpace, Parms...>& btm,
            const typename Solve<ExeSpace, Parms...>::ConstVector& b,
            const typename Solve<ExeSpace, Parms...>::Vector& x,
            const TridiagImpl::Enum impl = (impl::is_cuda<ExeSpace>::value ? TridiagImpl::team :
                                            impl::mimic_sparc ? TridiagImpl::sparc :
                                            TridiagImpl::serial)) {
  switch (impl) {
  case TridiagImpl::team:
    SolveFusedTeamByLine<ExeSpace, Parms...>(btm).solve(b, x);
    break;
  case TridiagImpl::sparc:
    sparc::solve(btm, b, x);
    break;
  case TridiagImpl::serial:
  default:
    SolveSerialByLine<ExeSpace, Parms...>(btm).solve(b, x);
  }
}

// Check the factorization given in a_fac by solving an equation. a holds the
// original tridiag matrices prior to factorization in place.
template <typename ExeSpace, typename ...Parms>
bool check_factor_and_solve (const BlockTridiagMatrices<ExeSpace, Parms...>& a,
                             const BlockTridiagMatrices<ExeSpace, Parms...>& a_fac,
                             const Int nrhs) {
  typedef Solve<ExeSpace, Parms...> SolveT;
  typedef typename SolveT::Vector Vector;
  typedef typename SolveT::ConstVector ConstVector;
  typedef typename SolveT::Scalar Scalar;
  const auto n = a.nblkrows_matrix() * a.blocksize();
  Vector x("x", n, nrhs), b("b", n, nrhs);
  fill_vector(b);
  bool pass;
  if (impl::is_cuda<ExeSpace>::value) {
    solve(a_fac, b, x, TridiagImpl::team);
    pass = check_solve(a, b, x, 1e4*std::numeric_limits<Scalar>::epsilon());
    solve(a_fac, b, x, TridiagImpl::serial);
    pass = pass && check_solve(a, b, x, 1e4*std::numeric_limits<Scalar>::epsilon());
  } else {
    solve(a_fac, b, x);
    pass = check_solve(a, b, x, 1e4*std::numeric_limits<Scalar>::epsilon());
  }
  return pass;
}

// Run a bunch of unit tests.
void test (const Int ni, const Int nj, const Int nk, const Int bs, const Int nrhs) {
  bool success = true;
  StructuredBlock sb(ni, nj, nk);

  // Check StructuredBlock.
  for (Int c = 0; c < sb.size(); ++c) {
    Int i, j, k;
    sb.id2ijk(c, i, j, k);
    TEST_ASSERT(i >= 0 && i < sb.ni, success);
    TEST_ASSERT(j >= 0 && j < sb.nj, success);
    TEST_ASSERT(k >= 0 && k < sb.nk, success);
    TEST_ASSERT(sb.ijk2id(i, j, k) == c, success);
  }

  auto g = make_graph_for_structured_block(sb, StencilShape::cross);
  BcrsMatrix<ko::DefaultHostExecutionSpace> mh(g, bs);
  fill_matrix(mh);
  auto md = copy_to_device(mh);

  // Test MVP.
  {
    typedef BcrsMvp<>::Vector Vector;
    typedef BcrsMvp<>::ConstVector ConstVector;
    const auto n = md->g().nrow() * md->bs();
    Vector x("x", n, nrhs), y1("y1", n, nrhs), y2("y2", n, nrhs);
    fill_vector(x);
    mvp(*md, x, y1, MvpImpl::by_row);
    mvp(*md, x, y2, MvpImpl::by_block_row);
    TEST_ASSERT(reldif(y1, y2) <= 1e2*std::numeric_limits<Scalar>::epsilon(), success);
  }

  // Check A = L U.
  ko::View<Scalar**, ko::DefaultHostExecutionSpace> ac("ac", bs, bs);
  auto a = ko::subview(mh.v(), 0, ko::ALL(), ko::ALL());
  ko::View<Scalar*[1], ko::DefaultHostExecutionSpace> x("x", bs, 1), b("b", bs, 1);
  for (Int i = 0; i < bs; ++i) b(i,0) = 2*(urand() - 0.5);
  for (int trial = 0; trial < 3; ++trial) {
    for (Int i = 0; i < bs; ++i)
      for (Int j = 0; j < bs; ++j)
        ac(i,j) = a(i,j);
    switch (trial) {
    case 0:
      cholfac_serial(ac);
      cholslv_serial(ac, b, x);
      break;
    case 1:
      lunpfac_serial(ac);
      lunpslv_serial(ac, b, x);
      break;
    case 2:
#ifdef MIMIC_SPARC
      sparc::impl::factor_nopivot(ac.data(), bs);
      sparc::impl::fwdsub(ac.data(), bs, b.data(), x.data());
      sparc::impl::bcksub(ac.data(), bs, x.data(), x.data());
#endif
      break;
    default:
      assert(0);
    }
    Scalar num = 0, den = 0;
    for (Int i = 0; i < bs; ++i) {
      Scalar ax = 0;
      for (Int j = 0; j < bs; ++j)
        ax += a(i,j)*x(j,0);
      num += square(ax - b(i,0));
      den += square(b(i,0));
    }
    TEST_ASSERT(std::sqrt(num/den) <= 1e2*std::numeric_limits<Scalar>::epsilon(), success);
  }

  // Check block tridiag extraction.
  auto btm_fac = extract_block_tridiag(sb, *md);
  TEST_ASSERT(check_block_tridiag(sb, *md, *btm_fac), success);
  // Check factor and solve.
  auto btm_op = btm_fac->clone();
  factor(*btm_fac);
  TEST_ASSERT(check_factor_and_solve(*btm_op, *btm_fac, nrhs), success);
  if ( ! success)
    pr(puf(ni) pu(nj) pu(nk) pu(bs) pu(nrhs));
}

inline bool eq (const std::string& a, const char* const b1, const char* const b2 = 0) {
  return (a == std::string(b1) || (b2 && a == std::string(b2)) ||
          a == std::string("-") + std::string(b1));
}

// Command-line argument parser and holder.
struct Input {
  bool quiet, check;
  Int ni, nj, nk;
  Int bs; // block size
  Int nrhs; // #vectors in multivector
  StencilShape::Enum stencil_shape;

  Input (int argc, char** argv) {
    quiet = false;
    check = false;
    ni = nj = nk = 10;
    bs = 5;
    nrhs = 1;
    stencil_shape = StencilShape::cross;
    for (int i = 1; i < argc; ++i) {
      const std::string& token = argv[i];
      if (eq(token, "-nijk")) ni = nj = nk = std::atoi(argv[++i]);
      else if (eq(token, "-ni")) ni = std::atoi(argv[++i]);
      else if (eq(token, "-nj")) nj = std::atoi(argv[++i]);
      else if (eq(token, "-nk")) nk = std::atoi(argv[++i]);
      else if (eq(token, "-bs")) bs = std::atoi(argv[++i]);
      else if (eq(token, "-nrhs")) nrhs = std::atoi(argv[++i]);
      else if (eq(token, "-c", "-check")) check = true;
    }
    if (nk <= 1)
      throw std::runtime_error("k dimension is <= 1; must be >= 2.");
    if ( ! quiet) print(std::cout);
  }

  void print (std::ostream& os) const {
    os << "<I> ni " << ni << " nj " << nj << " nk " << nk
       << " bs " << bs
       << " nrhs " << nrhs
       << " sc " << stencil_shape << "\n";
  }
};

int run (const Input& in) {
  StructuredBlock sb(in.ni, in.nj, in.nk);
  std::shared_ptr<BcrsMatrix<> > A;
  {
    auto g = make_graph_for_structured_block(sb, in.stencil_shape);
    BcrsMatrix<ko::DefaultHostExecutionSpace> Ah(g, in.bs);
    fill_matrix(Ah);
    A = copy_to_device(Ah);
  }
  int dontopt = 0; // Make sure repetition loops don't get opt'ed away.
  double et_mvp;
  static const int nmvp = 50;
  {
    typedef BcrsMvp<>::Vector Vector;
    typedef BcrsMvp<>::ConstVector ConstVector;
    const auto n = A->g().nrow() * A->bs();
    Vector x("x", n, in.nrhs), y("y", n, in.nrhs);
    fill_vector(x);
    {
      SimpleTimer t("50 MVP");
      for (int i = 0; i < nmvp; ++i) {
        mvp(*A, x, y);
        dontopt += i;
      }
      ko::fence();
      et_mvp = t.elapsed_time();
    }
  }
  bool success = true;
  std::shared_ptr<BlockTridiagMatrices<> > T_fac;
  double et_extract;
  {
    SimpleTimer t("extract_block_tridiag");
    T_fac = extract_block_tridiag(sb, *A);
    ko::fence();
    et_extract = t.elapsed_time();
  }
  if (in.check) TEST_ASSERT(check_block_tridiag(sb, *A, *T_fac), success);
  double et_factor;
  {
    SimpleTimer t("factor");
    factor(*T_fac);
    ko::fence();
    et_factor = t.elapsed_time();
  }
  {
    typedef BcrsMatrix<>::Scalar Scalar;
    typedef BcrsMvp<>::Vector Vector;
    typedef BcrsMvp<>::ConstVector ConstVector;
    const auto n = T_fac->nblkrows_matrix() * T_fac->blocksize();
    Vector x("x", n, in.nrhs), b("b", n, in.nrhs);
    fill_vector(b);
    ConstVector bc(b);
    double et_solve;
    static const int nsolves = 50;
    {
      SimpleTimer t("50 solves");
      for (int i = 0; i < nsolves; ++i) {
        solve(*T_fac, bc, x);
        dontopt += i;
      }
      ko::fence();
      et_solve = t.elapsed_time();
    }
    std::cout << "   extract/mvp = " << et_extract / (et_mvp/nmvp) << "\n";
    std::cout << "    factor/mvp = " << et_factor / (et_mvp/nmvp) << "\n";
    std::cout << "     solve/mvp = " << (et_solve/nsolves) / (et_mvp/nmvp) << "\n";
    if (in.check) {
      const auto T_op = extract_block_tridiag(sb, *A);
      TEST_ASSERT(check_factor_and_solve(*T_op, *T_fac, in.nrhs), success);
    }
  }
  if (in.check && ! success) std::cerr << "FAILED: run\n";
  return dontopt;
}

int main (int argc, char** argv) {
  int ret;
  ko::initialize(argc, argv); {
    int nt = 0;
#pragma omp parallel 
    nt = omp_get_num_threads();
    
    std::cout << " nthreads = " << nt << std::endl;
    test(3, 4, 2, 5, 1);
    test(3, 4, 2, 25, 2);
    test(44, 63, 15, 4, 1);
    test(2, 2, 15, 3, 3);
    for (Int nrhs = 1; nrhs <= 33; ++nrhs)
      test(2, 2, 15, 3, nrhs);
    test(1, 1, 2, 63, 8);
    Input in(argc, argv);
    ret = run(in);
    ret = 0;
  } ko::finalize();
  return ret;
}
