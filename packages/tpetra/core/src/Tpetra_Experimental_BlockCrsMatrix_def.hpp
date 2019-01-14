// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_EXPERIMENTAL_BLOCKCRSMATRIX_DEF_HPP
#define TPETRA_EXPERIMENTAL_BLOCKCRSMATRIX_DEF_HPP

/// \file Tpetra_Experimental_BlockCrsMatrix_def.hpp
/// \brief Definition of Tpetra::Experimental::BlockCrsMatrix

#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_PackTraits.hpp"
#include "Tpetra_Details_Profiling.hpp"

#include "Teuchos_TimeMonitor.hpp"
#ifdef HAVE_TPETRA_DEBUG
#  include <set>
#endif // HAVE_TPETRA_DEBUG

//
// mfh 25 May 2016: Temporary fix for #393.
//
// Don't use lambdas in the BCRS mat-vec for GCC < 4.8, due to a GCC
// 4.7.2 compiler bug ("internal compiler error") when compiling them.
// Also, lambdas for Kokkos::parallel_* don't work with CUDA, so don't
// use them in that case, either.
//
// mfh 31 May 2016: GCC 4.9.[23] appears to be broken ("internal
// compiler error") too.  Ditto for GCC 5.1.  I'll just disable the
// thing for any GCC version.
//
#if defined(__CUDACC__)
   // Lambdas for Kokkos::parallel_* don't work with CUDA 7.5 either.
#  if defined(TPETRA_BLOCKCRSMATRIX_APPLY_USE_LAMBDA)
#    undef TPETRA_BLOCKCRSMATRIX_APPLY_USE_LAMBDA
#  endif // defined(TPETRA_BLOCKCRSMATRIX_APPLY_USE_LAMBDA)

#elif defined(__GNUC__)

#  if defined(TPETRA_BLOCKCRSMATRIX_APPLY_USE_LAMBDA)
#    undef TPETRA_BLOCKCRSMATRIX_APPLY_USE_LAMBDA
#  endif // defined(TPETRA_BLOCKCRSMATRIX_APPLY_USE_LAMBDA)

#else // some other compiler

   // Optimistically assume that other compilers aren't broken.
#  if ! defined(TPETRA_BLOCKCRSMATRIX_APPLY_USE_LAMBDA)
#    define TPETRA_BLOCKCRSMATRIX_APPLY_USE_LAMBDA 1
#  endif // ! defined(TPETRA_BLOCKCRSMATRIX_APPLY_USE_LAMBDA)
#endif // defined(__CUDACC__), defined(__GNUC__)


namespace Tpetra {
namespace Experimental {

namespace Impl {

  template<typename T>
  struct BlockCrsRowStruct {
    T totalNumEntries, totalNumBytes, maxRowLength;
    KOKKOS_INLINE_FUNCTION BlockCrsRowStruct() = default;
    KOKKOS_INLINE_FUNCTION BlockCrsRowStruct(const BlockCrsRowStruct &b) = default;
    KOKKOS_INLINE_FUNCTION BlockCrsRowStruct(const T& numEnt, const T& numBytes, const T& rowLength)
      : totalNumEntries(numEnt), totalNumBytes(numBytes), maxRowLength(rowLength) {}
  };

  template<typename T>
  static
  KOKKOS_INLINE_FUNCTION
  void operator+=(volatile BlockCrsRowStruct<T> &a,
                  volatile const BlockCrsRowStruct<T> &b) {
    a.totalNumEntries += b.totalNumEntries;
    a.totalNumBytes   += b.totalNumBytes;
    a.maxRowLength     = a.maxRowLength > b.maxRowLength ? a.maxRowLength : b.maxRowLength;
  }

  template<typename T>
  static
  KOKKOS_INLINE_FUNCTION
  void operator+=(BlockCrsRowStruct<T> &a, const BlockCrsRowStruct<T> &b) {
    a.totalNumEntries += b.totalNumEntries;
    a.totalNumBytes   += b.totalNumBytes;
    a.maxRowLength     = a.maxRowLength > b.maxRowLength ? a.maxRowLength : b.maxRowLength;
  }

  template<typename T, typename ExecSpace>
  struct BlockCrsReducer {
    typedef BlockCrsReducer reducer;
    typedef T value_type;
    typedef Kokkos::View<value_type,ExecSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;
    value_type *value;

    KOKKOS_INLINE_FUNCTION
    BlockCrsReducer(value_type &val) : value(&val) {}

    KOKKOS_INLINE_FUNCTION void join(value_type &dst, value_type &src) const { dst += src; }
    KOKKOS_INLINE_FUNCTION void join(volatile value_type &dst, const volatile value_type &src) const { dst += src; }
    KOKKOS_INLINE_FUNCTION void init(value_type &val) const { val = value_type(); }
    KOKKOS_INLINE_FUNCTION value_type& reference() { return *value; }
    KOKKOS_INLINE_FUNCTION result_view_type view() const { return result_view_type(value); }
  };

#if 0
template<class AlphaCoeffType,
         class GraphType,
         class MatrixValuesType,
         class InVecType,
         class BetaCoeffType,
         class OutVecType>
class BcrsApplyNoTrans1VecTeamFunctor {
private:
  static_assert (Kokkos::Impl::is_view<MatrixValuesType>::value,
                 "MatrixValuesType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<OutVecType>::value,
                 "OutVecType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<InVecType>::value,
                 "InVecType must be a Kokkos::View.");
  static_assert (std::is_same<MatrixValuesType,
                   typename MatrixValuesType::const_type>::value,
                 "MatrixValuesType must be a const Kokkos::View.");
  static_assert (std::is_same<OutVecType,
                   typename OutVecType::non_const_type>::value,
                 "OutVecType must be a nonconst Kokkos::View.");
  static_assert (std::is_same<InVecType, typename InVecType::const_type>::value,
                 "InVecType must be a const Kokkos::View.");
  static_assert (static_cast<int> (MatrixValuesType::rank) == 1,
                 "MatrixValuesType must be a rank-1 Kokkos::View.");
  static_assert (static_cast<int> (InVecType::rank) == 1,
                 "InVecType must be a rank-1 Kokkos::View.");
  static_assert (static_cast<int> (OutVecType::rank) == 1,
                 "OutVecType must be a rank-1 Kokkos::View.");
  typedef typename MatrixValuesType::non_const_value_type scalar_type;
  typedef typename GraphType::device_type device_type;
  typedef typename device_type::execution_space execution_space;
  typedef typename execution_space::scratch_memory_space shmem_space;

public:
  //! Type of the (mesh) column indices in the sparse graph / matrix.
  typedef typename std::remove_const<typename GraphType::data_type>::type
    local_ordinal_type;
  //! Use this for the Kokkos::parallel_for policy argument.
  typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Dynamic>,
                             execution_space,
                             local_ordinal_type> policy_type;
  /// \brief Set the current vector / current column of the input
  ///   (multi)vector X to use.
  ///
  /// This lets us handle multiple columns by iterating over them one
  /// column at a time, without needing to recreate the functor each
  /// time.
  void setX (const InVecType& X) { X_ = X; }

  /// \brief Set the current vector / current column of the output
  ///   (multi)vector Y to use.
  ///
  /// This lets us handle multiple columns by iterating over them one
  /// column at a time, without needing to recreate the functor each
  /// time.
  void setY (const OutVecType& Y) { Y_ = Y; }

  /// \brief Get the per-team scratch size (for setting up the TeamPolicy).
  ///
  /// This is used only for the team parallel_reduce.
  KOKKOS_INLINE_FUNCTION local_ordinal_type
  getScratchSizePerTeam () const
  {
    // WARNING (mfh 01 Jun 2016) This may not work for Scalar types
    // that have run-time sizes, like some of those in Stokhos.
    typedef typename std::decay<decltype (Y_(0))>::type y_value_type;
    return blockSize_ * sizeof (y_value_type);
  }

  /// \brief Get the per-thread scratch size (for setting up the TeamPolicy).
  ///
  /// This is used only for the team parallel_for.
  KOKKOS_INLINE_FUNCTION local_ordinal_type
  getScratchSizePerThread () const
  {
    // WARNING (mfh 01 Jun 2016) This may not work for Scalar types
    // that have run-time sizes, like some of those in Stokhos.
    typedef typename std::decay<decltype (Y_(0))>::type y_value_type;
    return blockSize_ * sizeof (y_value_type);
  }

private:
  KOKKOS_INLINE_FUNCTION local_ordinal_type
  getNumLclMeshRows () const
  {
    return ptr_.extent (0) == 0 ?
      static_cast<local_ordinal_type> (0) :
      static_cast<local_ordinal_type> (ptr_.extent (0) - 1);
  }

  static constexpr local_ordinal_type defaultRowsPerTeam = 20;

public:
  //! Get the number of teams (first argument of the TeamPolicy).
  local_ordinal_type getNumTeams () const {
    return 1;
    // const local_ordinal_type numLclMeshRows = getNumLclMeshRows ();
    // return (numLclMeshRows + rowsPerTeam_ - 1) / rowsPerTeam_;
  }

  //! Constructor.
  BcrsApplyNoTrans1VecTeamFunctor (const typename std::decay<AlphaCoeffType>::type& alpha,
                                   const GraphType& graph,
                                   const MatrixValuesType& val,
                                   const local_ordinal_type blockSize,
                                   const InVecType& X,
                                   const typename std::decay<BetaCoeffType>::type& beta,
                                   const OutVecType& Y,
                                   const local_ordinal_type rowsPerTeam = defaultRowsPerTeam) :
    alpha_ (alpha),
    ptr_ (graph.row_map),
    ind_ (graph.entries),
    val_ (val),
    blockSize_ (blockSize),
    X_ (X),
    beta_ (beta),
    Y_ (Y),
    rowsPerTeam_ (rowsPerTeam)
  {
    // std::ostringstream os;
    // os << "Functor ctor: numLclMeshRows = " << getNumLclMeshRows () << std::endl;
    // std::cerr << os.str ();
  }

  KOKKOS_INLINE_FUNCTION void
  operator () (const typename policy_type::member_type& member) const
  {
    using ::Tpetra::Experimental::COPY;
    using ::Tpetra::Experimental::FILL;
    using ::Tpetra::Experimental::SCAL;
    using ::Tpetra::Experimental::GEMV;
    using Kokkos::Details::ArithTraits;
    // I'm not writing 'using Kokkos::make_pair;' here, because that
    // may break builds for users who make the mistake of putting
    // 'using namespace std;' in the global namespace.  Please don't
    // ever do that!  But just in case you do, I'll take this
    // precaution.
    using Kokkos::parallel_for;
    using Kokkos::subview;
    typedef local_ordinal_type LO;
    typedef typename decltype (ptr_)::non_const_value_type offset_type;
    typedef Kokkos::View<typename OutVecType::non_const_value_type*,
      shmem_space, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      shared_array_type;
    typedef Kokkos::View<typename OutVecType::non_const_value_type*,
      Kokkos::LayoutRight,
      device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      out_little_vec_type;
    typedef Kokkos::View<typename MatrixValuesType::const_value_type**,
      Kokkos::LayoutRight,
      device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      little_block_type;

    const LO leagueRank = member.league_rank();

    // This looks wrong!  All the threads are writing into the same local storage.
    // shared_array_type threadLocalMem =
    //   shared_array_type (member.thread_scratch (1), blockSize_);
    shared_array_type threadLocalMem =
      shared_array_type (member.thread_scratch (1), blockSize_ * rowsPerTeam_);

    // This looks wrong!  All the threads are writing into the same local storage.
    //out_little_vec_type Y_tlm (threadLocalMem.data (), blockSize_, 1);

    const LO numLclMeshRows = getNumLclMeshRows ();
    const LO rowBeg = leagueRank * rowsPerTeam_;
    const LO rowTmp = rowBeg + rowsPerTeam_;
    const LO rowEnd = rowTmp < numLclMeshRows ? rowTmp : numLclMeshRows;

    // {
    //   std::ostringstream os;
    //   os << leagueRank << "," << member.team_rank () << ": "
    //      << rowBeg << "," << rowEnd << std::endl;
    //   std::cerr << os.str ();
    // }

    // Each team takes rowsPerTeam_ (block) rows.
    // Each thread in the team takes a single (block) row.
    parallel_for (Kokkos::TeamThreadRange (member, rowBeg, rowEnd),
                  [&] (const LO& lclRow) {
                    // Each thread in the team gets its own temp storage.
                    out_little_vec_type Y_tlm (threadLocalMem.data () + blockSize_ * member.team_rank (), blockSize_);

                    const offset_type Y_ptBeg = lclRow * blockSize_;
                    const offset_type Y_ptEnd = Y_ptBeg + blockSize_;
                    auto Y_cur =
                      subview (Y_, ::Kokkos::make_pair (Y_ptBeg, Y_ptEnd));
                    if (beta_ == ArithTraits<BetaCoeffType>::zero ()) {
                      FILL (Y_tlm, ArithTraits<BetaCoeffType>::zero ());
                    }
                    else if (beta_ == ArithTraits<BetaCoeffType>::one ()) {
                      COPY (Y_cur, Y_tlm);
                    }
                    else { // beta != 0 && beta != 1
                      COPY (Y_cur, Y_tlm);
                      SCAL (beta_, Y_tlm);
                    }

                    if (alpha_ != ArithTraits<AlphaCoeffType>::zero ()) {
                      const offset_type blkBeg = ptr_[lclRow];
                      const offset_type blkEnd = ptr_[lclRow+1];
                      // Precompute to save integer math in the inner loop.
                      const offset_type bs2 = blockSize_ * blockSize_;
                      for (offset_type absBlkOff = blkBeg; absBlkOff < blkEnd;
                           ++absBlkOff) {
                        little_block_type A_cur (val_.data () + absBlkOff * bs2,
                                                 blockSize_, blockSize_);
                        const offset_type X_blkCol = ind_[absBlkOff];
                        const offset_type X_ptBeg = X_blkCol * blockSize_;
                        const offset_type X_ptEnd = X_ptBeg + blockSize_;
                        auto X_cur =
                          subview (X_, ::Kokkos::make_pair (X_ptBeg, X_ptEnd));
                        // Y_tlm += alpha*A_cur*X_cur
                        GEMV (alpha_, A_cur, X_cur, Y_tlm);
                      } // for each entry in current local block row of matrix
                      COPY (Y_tlm, Y_cur);
                    }
                  });
  }

private:
  typename std::decay<AlphaCoeffType>::type alpha_;
  typename GraphType::row_map_type::const_type ptr_;
  typename GraphType::entries_type::const_type ind_;
  MatrixValuesType val_;
  local_ordinal_type blockSize_;
  InVecType X_;
  typename std::decay<BetaCoeffType>::type beta_;
  OutVecType Y_;
  local_ordinal_type rowsPerTeam_;
};
#endif // 0

template<class AlphaCoeffType,
         class GraphType,
         class MatrixValuesType,
         class InVecType,
         class BetaCoeffType,
         class OutVecType>
class BcrsApplyNoTrans1VecFunctor {
private:
  static_assert (Kokkos::Impl::is_view<MatrixValuesType>::value,
                 "MatrixValuesType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<OutVecType>::value,
                 "OutVecType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<InVecType>::value,
                 "InVecType must be a Kokkos::View.");
  static_assert (std::is_same<MatrixValuesType,
                   typename MatrixValuesType::const_type>::value,
                 "MatrixValuesType must be a const Kokkos::View.");
  static_assert (std::is_same<OutVecType,
                   typename OutVecType::non_const_type>::value,
                 "OutVecType must be a nonconst Kokkos::View.");
  static_assert (std::is_same<InVecType, typename InVecType::const_type>::value,
                 "InVecType must be a const Kokkos::View.");
  static_assert (static_cast<int> (MatrixValuesType::rank) == 1,
                 "MatrixValuesType must be a rank-1 Kokkos::View.");
  static_assert (static_cast<int> (InVecType::rank) == 1,
                 "InVecType must be a rank-1 Kokkos::View.");
  static_assert (static_cast<int> (OutVecType::rank) == 1,
                 "OutVecType must be a rank-1 Kokkos::View.");
  typedef typename MatrixValuesType::non_const_value_type scalar_type;

public:
  typedef typename GraphType::device_type device_type;

  //! Type of the (mesh) column indices in the sparse graph / matrix.
  typedef typename std::remove_const<typename GraphType::data_type>::type
    local_ordinal_type;
  //! Use this for the Kokkos::parallel_for policy argument.
  typedef Kokkos::RangePolicy<Kokkos::Schedule<Kokkos::Dynamic>,
                              typename device_type::execution_space,
                              local_ordinal_type> policy_type;
  /// \brief Set the current vector / current column of the input
  ///   (multi)vector X to use.
  ///
  /// This lets us handle multiple columns by iterating over them one
  /// column at a time, without needing to recreate the functor each
  /// time.
  void setX (const InVecType& X) { X_ = X; }

  /// \brief Set the current vector / current column of the output
  ///   (multi)vector Y to use.
  ///
  /// This lets us handle multiple columns by iterating over them one
  /// column at a time, without needing to recreate the functor each
  /// time.
  void setY (const OutVecType& Y) { Y_ = Y; }

  typedef typename Kokkos::ArithTraits<typename std::decay<AlphaCoeffType>::type>::val_type alpha_coeff_type;
  typedef typename Kokkos::ArithTraits<typename std::decay<BetaCoeffType>::type>::val_type beta_coeff_type;

  //! Constructor.
  BcrsApplyNoTrans1VecFunctor (const alpha_coeff_type& alpha,
                               const GraphType& graph,
                               const MatrixValuesType& val,
                               const local_ordinal_type blockSize,
                               const InVecType& X,
                               const beta_coeff_type& beta,
                               const OutVecType& Y) :
    alpha_ (alpha),
    ptr_ (graph.row_map),
    ind_ (graph.entries),
    val_ (val),
    blockSize_ (blockSize),
    X_ (X),
    beta_ (beta),
    Y_ (Y)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator () (const local_ordinal_type& lclRow) const
  {
    using ::Tpetra::Experimental::COPY;
    using ::Tpetra::Experimental::FILL;
    using ::Tpetra::Experimental::SCAL;
    using ::Tpetra::Experimental::GEMV;
    using Kokkos::Details::ArithTraits;
    // I'm not writing 'using Kokkos::make_pair;' here, because that
    // may break builds for users who make the mistake of putting
    // 'using namespace std;' in the global namespace.  Please don't
    // ever do that!  But just in case you do, I'll take this
    // precaution.
    using Kokkos::parallel_for;
    using Kokkos::subview;
    typedef typename decltype (ptr_)::non_const_value_type offset_type;
    typedef Kokkos::View<typename MatrixValuesType::const_value_type**,
      Kokkos::LayoutRight,
      device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      little_block_type;

    const offset_type Y_ptBeg = lclRow * blockSize_;
    const offset_type Y_ptEnd = Y_ptBeg + blockSize_;
    auto Y_cur = subview (Y_, ::Kokkos::make_pair (Y_ptBeg, Y_ptEnd));

    // This version of the code does not use temporary storage.
    // Each thread writes to its own block of the target vector.
    if (beta_ == ArithTraits<beta_coeff_type>::zero ()) {
      FILL (Y_cur, ArithTraits<beta_coeff_type>::zero ());
    }
    else if (beta_ != ArithTraits<beta_coeff_type>::one ()) { // beta != 0 && beta != 1
      SCAL (beta_, Y_cur);
    }

    if (alpha_ != ArithTraits<alpha_coeff_type>::zero ()) {
      const offset_type blkBeg = ptr_[lclRow];
      const offset_type blkEnd = ptr_[lclRow+1];
      // Precompute to save integer math in the inner loop.
      const offset_type bs2 = blockSize_ * blockSize_;
      for (offset_type absBlkOff = blkBeg; absBlkOff < blkEnd;
           ++absBlkOff) {
        little_block_type A_cur (val_.data () + absBlkOff * bs2,
                                 blockSize_, blockSize_);
        const offset_type X_blkCol = ind_[absBlkOff];
        const offset_type X_ptBeg = X_blkCol * blockSize_;
        const offset_type X_ptEnd = X_ptBeg + blockSize_;
        auto X_cur = subview (X_, ::Kokkos::make_pair (X_ptBeg, X_ptEnd));

        GEMV (alpha_, A_cur, X_cur, Y_cur); // Y_cur += alpha*A_cur*X_cur
      } // for each entry in current local block row of matrix
    }
  }

private:
  alpha_coeff_type alpha_;
  typename GraphType::row_map_type::const_type ptr_;
  typename GraphType::entries_type::const_type ind_;
  MatrixValuesType val_;
  local_ordinal_type blockSize_;
  InVecType X_;
  beta_coeff_type beta_;
  OutVecType Y_;
};

template<class AlphaCoeffType,
         class GraphType,
         class MatrixValuesType,
         class InMultiVecType,
         class BetaCoeffType,
         class OutMultiVecType>
void
bcrsLocalApplyNoTrans (const AlphaCoeffType& alpha,
                       const GraphType& graph,
                       const MatrixValuesType& val,
                       const typename std::remove_const<typename GraphType::data_type>::type blockSize,
                       const InMultiVecType& X,
                       const BetaCoeffType& beta,
                       const OutMultiVecType& Y
#if 0
                       , const typename std::remove_const<typename GraphType::data_type>::type rowsPerTeam = 20
#endif // 0
                       )
{
  static_assert (Kokkos::Impl::is_view<MatrixValuesType>::value,
                 "MatrixValuesType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<OutMultiVecType>::value,
                 "OutMultiVecType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<InMultiVecType>::value,
                 "InMultiVecType must be a Kokkos::View.");
  static_assert (static_cast<int> (MatrixValuesType::rank) == 1,
                 "MatrixValuesType must be a rank-1 Kokkos::View.");
  static_assert (static_cast<int> (OutMultiVecType::rank) == 2,
                 "OutMultiVecType must be a rank-2 Kokkos::View.");
  static_assert (static_cast<int> (InMultiVecType::rank) == 2,
                 "InMultiVecType must be a rank-2 Kokkos::View.");

  typedef typename MatrixValuesType::const_type matrix_values_type;
  typedef typename OutMultiVecType::non_const_type out_multivec_type;
  typedef typename InMultiVecType::const_type in_multivec_type;
  typedef typename Kokkos::ArithTraits<typename std::decay<AlphaCoeffType>::type>::val_type alpha_type;
  typedef typename Kokkos::ArithTraits<typename std::decay<BetaCoeffType>::type>::val_type beta_type;
  typedef typename std::remove_const<typename GraphType::data_type>::type LO;

  const LO numLocalMeshRows = graph.row_map.extent (0) == 0 ?
    static_cast<LO> (0) :
    static_cast<LO> (graph.row_map.extent (0) - 1);
  const LO numVecs = Y.extent (1);
  if (numLocalMeshRows == 0 || numVecs == 0) {
    return; // code below doesn't handle numVecs==0 correctly
  }

  // These assignments avoid instantiating the functor extra times
  // unnecessarily, e.g., for X const vs. nonconst.  We only need the
  // X const case, so only instantiate for that case.
  in_multivec_type X_in = X;
  out_multivec_type Y_out = Y;

  // The functor only knows how to handle one vector at a time, and it
  // expects 1-D Views.  Thus, we need to know the type of each column
  // of X and Y.
  typedef decltype (Kokkos::subview (X_in, Kokkos::ALL (), 0)) in_vec_type;
  typedef decltype (Kokkos::subview (Y_out, Kokkos::ALL (), 0)) out_vec_type;
#if 0
  typedef BcrsApplyNoTrans1VecTeamFunctor<alpha_type, GraphType,
    matrix_values_type, in_vec_type, beta_type, out_vec_type> functor_type;
#else
  typedef BcrsApplyNoTrans1VecFunctor<alpha_type, GraphType,
    matrix_values_type, in_vec_type, beta_type, out_vec_type> functor_type;
#endif // 0
  typedef typename functor_type::policy_type policy_type;

  auto X_0 = Kokkos::subview (X_in, Kokkos::ALL (), 0);
  auto Y_0 = Kokkos::subview (Y_out, Kokkos::ALL (), 0);
#if 0
  functor_type functor (alpha, graph, val, blockSize, X_0, beta, Y_0, rowsPerTeam);
  const LO numTeams = functor.getNumTeams ();
  policy_type policy (numTeams, Kokkos::AUTO ());
  {
    // KJ : hierarchy level of memory allocated e.g., cache (1),
    // HBM (2), DDR (3), not used for now
    const LO level = 1;
    // KJ : for now provide two options for parallelizing (for vs. reduce)
    const LO scratchSizePerTeam   = functor.getScratchSizePerTeam (); // used for team parallel_red
    const LO scratchSizePerThread = functor.getScratchSizePerThread (); // used for team parallel_for
    policy =
      policy.set_scratch_size (level,
                               Kokkos::PerTeam (scratchSizePerTeam),
                               Kokkos::PerThread (scratchSizePerThread));
  }
#else
  functor_type functor (alpha, graph, val, blockSize, X_0, beta, Y_0);
  policy_type policy (0, numLocalMeshRows);
#endif // 0

  // Compute the first column of Y.
  Kokkos::parallel_for (policy, functor);

  // Compute the remaining columns of Y.
  for (LO j = 1; j < numVecs; ++j) {
    auto X_j = Kokkos::subview (X_in, Kokkos::ALL (), j);
    auto Y_j = Kokkos::subview (Y_out, Kokkos::ALL (), j);
    functor.setX (X_j);
    functor.setY (Y_j);
    Kokkos::parallel_for (policy, functor);
  }
}

} // namespace Impl

namespace { // (anonymous)

// Implementation of BlockCrsMatrix::getLocalDiagCopy (non-deprecated
// version that takes two Kokkos::View arguments).
template<class Scalar, class LO, class GO, class Node>
class GetLocalDiagCopy {
public:
  typedef typename Node::device_type device_type;
  typedef size_t diag_offset_type;
  typedef Kokkos::View<const size_t*, device_type,
                       Kokkos::MemoryUnmanaged> diag_offsets_type;
  typedef typename ::Tpetra::CrsGraph<LO, GO, Node> global_graph_type;
  typedef typename global_graph_type::local_graph_type local_graph_type;
  typedef typename local_graph_type::row_map_type row_offsets_type;
  typedef typename ::Tpetra::Experimental::BlockMultiVector<Scalar, LO, GO, Node>::impl_scalar_type IST;
  typedef Kokkos::View<IST***, device_type, Kokkos::MemoryUnmanaged> diag_type;
  typedef Kokkos::View<const IST*, device_type, Kokkos::MemoryUnmanaged> values_type;

  // Constructor
  GetLocalDiagCopy (const diag_type& diag,
                    const values_type& val,
                    const diag_offsets_type& diagOffsets,
                    const row_offsets_type& ptr,
                    const LO blockSize) :
    diag_ (diag),
    diagOffsets_ (diagOffsets),
    ptr_ (ptr),
    blockSize_ (blockSize),
    offsetPerBlock_ (blockSize_*blockSize_),
    val_(val)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator() (const LO& lclRowInd) const
  {
    using Kokkos::ALL;

    // Get row offset
    const size_t absOffset = ptr_[lclRowInd];

    // Get offset relative to start of row
    const size_t relOffset = diagOffsets_[lclRowInd];

    // Get the total offset
    const size_t pointOffset = (absOffset+relOffset)*offsetPerBlock_;

    // Get a view of the block.  BCRS currently uses LayoutRight
    // regardless of the device.
    typedef Kokkos::View<const IST**, Kokkos::LayoutRight,
      device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      const_little_block_type;
    const_little_block_type D_in (val_.data () + pointOffset,
                                  blockSize_, blockSize_);
    auto D_out = Kokkos::subview (diag_, lclRowInd, ALL (), ALL ());
    COPY (D_in, D_out);
  }

  private:
    diag_type diag_;
    diag_offsets_type diagOffsets_;
    row_offsets_type ptr_;
    LO blockSize_;
    LO offsetPerBlock_;
    values_type val_;
  };
} // namespace (anonymous)

  template<class Scalar, class LO, class GO, class Node>
  std::ostream&
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  markLocalErrorAndGetStream ()
  {
    * (this->localError_) = true;
    if ((*errs_).is_null ()) {
      *errs_ = Teuchos::rcp (new std::ostringstream ());
    }
    return **errs_;
  }

  template<class Scalar, class LO, class GO, class Node>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  BlockCrsMatrix () :
    dist_object_type (Teuchos::rcp (new map_type ())), // nonnull, so DistObject doesn't throw
    graph_ (Teuchos::rcp (new map_type ()), 0), // FIXME (mfh 16 May 2014) no empty ctor yet
    blockSize_ (static_cast<LO> (0)),
    X_colMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    Y_rowMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    pointImporter_ (new Teuchos::RCP<typename crs_graph_type::import_type> ()),
    offsetPerBlock_ (0),
    localError_ (new bool (false)),
    errs_ (new Teuchos::RCP<std::ostringstream> ()) // ptr to a null ptr
  {
  }

  template<class Scalar, class LO, class GO, class Node>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  BlockCrsMatrix (const crs_graph_type& graph,
                  const LO blockSize) :
    dist_object_type (graph.getMap ()),
    graph_ (graph),
    rowMeshMap_ (* (graph.getRowMap ())),
    blockSize_ (blockSize),
    X_colMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    Y_rowMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    pointImporter_ (new Teuchos::RCP<typename crs_graph_type::import_type> ()),
    offsetPerBlock_ (blockSize * blockSize),
    localError_ (new bool (false)),
    errs_ (new Teuchos::RCP<std::ostringstream> ()) // ptr to a null ptr
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph_.isSorted (), std::invalid_argument, "Tpetra::Experimental::"
      "BlockCrsMatrix constructor: The input CrsGraph does not have sorted "
      "rows (isSorted() is false).  This class assumes sorted rows.");

    graphRCP_ = Teuchos::rcpFromRef(graph_);

    // Trick to test whether LO is nonpositive, without a compiler
    // warning in case LO is unsigned (which is generally a bad idea
    // anyway).  I don't promise that the trick works, but it
    // generally does with gcc at least, in my experience.
    const bool blockSizeIsNonpositive = (blockSize + 1 <= 1);
    TEUCHOS_TEST_FOR_EXCEPTION(
      blockSizeIsNonpositive, std::invalid_argument, "Tpetra::Experimental::"
      "BlockCrsMatrix constructor: The input blockSize = " << blockSize <<
      " <= 0.  The block size must be positive.");

    domainPointMap_ = BMV::makePointMap (* (graph.getDomainMap ()), blockSize);
    rangePointMap_ = BMV::makePointMap (* (graph.getRangeMap ()), blockSize);

    {
      typedef typename crs_graph_type::local_graph_type::row_map_type row_map_type;
      typedef typename row_map_type::HostMirror::non_const_type nc_host_row_map_type;

      row_map_type ptr_d = graph.getLocalGraph ().row_map;
      nc_host_row_map_type ptr_h_nc = Kokkos::create_mirror_view (ptr_d);
      Kokkos::deep_copy (ptr_h_nc, ptr_d);
      ptrHost_ = ptr_h_nc;
    }
    {
      typedef typename crs_graph_type::local_graph_type::entries_type entries_type;
      typedef typename entries_type::HostMirror::non_const_type nc_host_entries_type;

      entries_type ind_d = graph.getLocalGraph ().entries;
      nc_host_entries_type ind_h_nc = Kokkos::create_mirror_view (ind_d);
      Kokkos::deep_copy (ind_h_nc, ind_d);
      indHost_ = ind_h_nc;
    }

    const auto numValEnt = graph.getNodeNumEntries () * offsetPerBlock ();
    val_ = decltype (val_) ("val", numValEnt);
  }

  template<class Scalar, class LO, class GO, class Node>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  BlockCrsMatrix (const crs_graph_type& graph,
                  const map_type& domainPointMap,
                  const map_type& rangePointMap,
                  const LO blockSize) :
    dist_object_type (graph.getMap ()),
    graph_ (graph),
    rowMeshMap_ (* (graph.getRowMap ())),
    domainPointMap_ (domainPointMap),
    rangePointMap_ (rangePointMap),
    blockSize_ (blockSize),
    X_colMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    Y_rowMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    pointImporter_ (new Teuchos::RCP<typename crs_graph_type::import_type> ()),
    offsetPerBlock_ (blockSize * blockSize),
    localError_ (new bool (false)),
    errs_ (new Teuchos::RCP<std::ostringstream> ()) // ptr to a null ptr
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph_.isSorted (), std::invalid_argument, "Tpetra::Experimental::"
      "BlockCrsMatrix constructor: The input CrsGraph does not have sorted "
      "rows (isSorted() is false).  This class assumes sorted rows.");

    graphRCP_ = Teuchos::rcpFromRef(graph_);

    // Trick to test whether LO is nonpositive, without a compiler
    // warning in case LO is unsigned (which is generally a bad idea
    // anyway).  I don't promise that the trick works, but it
    // generally does with gcc at least, in my experience.
    const bool blockSizeIsNonpositive = (blockSize + 1 <= 1);
    TEUCHOS_TEST_FOR_EXCEPTION(
      blockSizeIsNonpositive, std::invalid_argument, "Tpetra::Experimental::"
      "BlockCrsMatrix constructor: The input blockSize = " << blockSize <<
      " <= 0.  The block size must be positive.");

    {
      typedef typename crs_graph_type::local_graph_type::row_map_type row_map_type;
      typedef typename row_map_type::HostMirror::non_const_type nc_host_row_map_type;

      row_map_type ptr_d = graph.getLocalGraph ().row_map;
      nc_host_row_map_type ptr_h_nc = Kokkos::create_mirror_view (ptr_d);
      Kokkos::deep_copy (ptr_h_nc, ptr_d);
      ptrHost_ = ptr_h_nc;
    }
    {
      typedef typename crs_graph_type::local_graph_type::entries_type entries_type;
      typedef typename entries_type::HostMirror::non_const_type nc_host_entries_type;

      entries_type ind_d = graph.getLocalGraph ().entries;
      nc_host_entries_type ind_h_nc = Kokkos::create_mirror_view (ind_d);
      Kokkos::deep_copy (ind_h_nc, ind_d);
      indHost_ = ind_h_nc;
    }

    const auto numValEnt = graph.getNodeNumEntries () * offsetPerBlock ();
    val_ = decltype (val_) ("val", numValEnt);
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const typename BlockCrsMatrix<Scalar, LO, GO, Node>::map_type>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getDomainMap () const
  { // Copy constructor of map_type does a shallow copy.
    // We're only returning by RCP for backwards compatibility.
    return Teuchos::rcp (new map_type (domainPointMap_));
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const typename BlockCrsMatrix<Scalar, LO, GO, Node>::map_type>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getRangeMap () const
  { // Copy constructor of map_type does a shallow copy.
    // We're only returning by RCP for backwards compatibility.
    return Teuchos::rcp (new map_type (rangePointMap_));
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const typename BlockCrsMatrix<Scalar, LO, GO, Node>::map_type>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getRowMap () const
  {
    return graph_.getRowMap();
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const typename BlockCrsMatrix<Scalar, LO, GO, Node>::map_type>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getColMap () const
  {
    return graph_.getColMap();
  }

  template<class Scalar, class LO, class GO, class Node>
  global_size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalNumRows() const
  {
    return graph_.getGlobalNumRows();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNodeNumRows() const
  {
    return graph_.getNodeNumRows();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNodeMaxNumRowEntries() const
  {
    return graph_.getNodeMaxNumRowEntries();
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  apply (const mv_type& X,
         mv_type& Y,
         Teuchos::ETransp mode,
         Scalar alpha,
         Scalar beta) const
  {
    typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_type;
    TEUCHOS_TEST_FOR_EXCEPTION(
      mode != Teuchos::NO_TRANS && mode != Teuchos::TRANS && mode != Teuchos::CONJ_TRANS,
      std::invalid_argument, "Tpetra::Experimental::BlockCrsMatrix::apply: "
      "Invalid 'mode' argument.  Valid values are Teuchos::NO_TRANS, "
      "Teuchos::TRANS, and Teuchos::CONJ_TRANS.");

    BMV X_view;
    BMV Y_view;
    const LO blockSize = getBlockSize ();
    try {
      X_view = BMV (X, * (graph_.getDomainMap ()), blockSize);
      Y_view = BMV (Y, * (graph_.getRangeMap ()), blockSize);
    }
    catch (std::invalid_argument& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::Experimental::BlockCrsMatrix::"
        "apply: Either the input MultiVector X or the output MultiVector Y "
        "cannot be viewed as a BlockMultiVector, given this BlockCrsMatrix's "
        "graph.  BlockMultiVector's constructor threw the following exception: "
        << e.what ());
    }

    try {
      // mfh 16 May 2014: Casting 'this' to nonconst is icky, but it's
      // either that or mark fields of this class as 'mutable'.  The
      // problem is that applyBlock wants to do lazy initialization of
      // temporary block multivectors.
      const_cast<this_type*> (this)->applyBlock (X_view, Y_view, mode, alpha, beta);
    } catch (std::invalid_argument& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::Experimental::BlockCrsMatrix::"
        "apply: The implementation method applyBlock complained about having "
        "an invalid argument.  It reported the following: " << e.what ());
    } catch (std::logic_error& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::Experimental::BlockCrsMatrix::"
        "apply: The implementation method applyBlock complained about a "
        "possible bug in its implementation.  It reported the following: "
        << e.what ());
    } catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::Experimental::BlockCrsMatrix::"
        "apply: The implementation method applyBlock threw an exception which "
        "is neither std::invalid_argument nor std::logic_error, but is a "
        "subclass of std::exception.  It reported the following: "
        << e.what ());
    } catch (...) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
        "apply: The implementation method applyBlock threw an exception which "
        "is not an instance of a subclass of std::exception.  This probably "
        "indicates a bug.  Please report this to the Tpetra developers.");
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  applyBlock (const BlockMultiVector<Scalar, LO, GO, Node>& X,
              BlockMultiVector<Scalar, LO, GO, Node>& Y,
              Teuchos::ETransp mode,
              const Scalar alpha,
              const Scalar beta)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getBlockSize () != Y.getBlockSize (), std::invalid_argument,
      "Tpetra::Experimental::BlockCrsMatrix::applyBlock: "
      "X and Y have different block sizes.  X.getBlockSize() = "
      << X.getBlockSize () << " != Y.getBlockSize() = "
      << Y.getBlockSize () << ".");

    if (mode == Teuchos::NO_TRANS) {
      applyBlockNoTrans (X, Y, alpha, beta);
    } else if (mode == Teuchos::TRANS || mode == Teuchos::CONJ_TRANS) {
      applyBlockTrans (X, Y, mode, alpha, beta);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::Experimental::BlockCrsMatrix::"
        "applyBlock: Invalid 'mode' argument.  Valid values are "
        "Teuchos::NO_TRANS, Teuchos::TRANS, and Teuchos::CONJ_TRANS.");
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  setAllToScalar (const Scalar& alpha)
  {
#ifdef HAVE_TPETRA_DEBUG
    const char prefix[] = "Tpetra::Experimental::BlockCrsMatrix::setAllToScalar: ";
#endif // HAVE_TPETRA_DEBUG

    if (this->template need_sync<device_type> ()) {
      // If we need to sync to device, then the data were last
      // modified on host.  In that case, we should again modify them
      // on host.
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION
        (this->need_sync_host (), std::runtime_error,
         prefix << "The matrix's values need sync on both device and host.");
#endif // HAVE_TPETRA_DEBUG
      this->modify_host ();
      Kokkos::deep_copy (getValuesHost (), alpha);
    }
    else if (this->need_sync_host ()) {
      // If we need to sync to host, then the data were last modified
      // on device.  In that case, we should again modify them on
      // device.
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION
        (this->template need_sync<device_type> (), std::runtime_error,
         prefix << "The matrix's values need sync on both host and device.");
#endif // HAVE_TPETRA_DEBUG
      this->template modify<device_type> ();
      Kokkos::deep_copy (this->template getValues<device_type> (), alpha);
    }
    else { // neither host nor device marked as modified, so modify on device
      this->template modify<device_type> ();
      Kokkos::deep_copy (this->template getValues<device_type> (), alpha);
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  replaceLocalValues (const LO localRowInd,
                      const LO colInds[],
                      const Scalar vals[],
                      const LO numColInds) const
  {
#ifdef HAVE_TPETRA_DEBUG
    const char prefix[] =
      "Tpetra::Experimental::BlockCrsMatrix::replaceLocalValues: ";
#endif // HAVE_TPETRA_DEBUG

    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We modified no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }
    const impl_scalar_type* const vIn =
      reinterpret_cast<const impl_scalar_type*> (vals);
    const size_t absRowBlockOffset = ptrHost_[localRowInd];
    const LO LINV = Teuchos::OrdinalTraits<LO>::invalid ();
    const LO perBlockSize = this->offsetPerBlock ();
    LO hint = 0; // Guess for the relative offset into the current row
    LO pointOffset = 0; // Current offset into input values
    LO validCount = 0; // number of valid column indices in colInds

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION
      (this->need_sync_host (), std::runtime_error,
       prefix << "The matrix's data were last modified on device, but have "
       "not been sync'd to host.  Please sync to host (by calling "
       "sync<Kokkos::HostSpace>() on this matrix) before calling this "
       "method.");
#endif // HAVE_TPETRA_DEBUG

    auto vals_host_out = getValuesHost ();
    impl_scalar_type* vals_host_out_raw = vals_host_out.data ();

    for (LO k = 0; k < numColInds; ++k, pointOffset += perBlockSize) {
      const LO relBlockOffset =
        this->findRelOffsetOfColumnIndex (localRowInd, colInds[k], hint);
      if (relBlockOffset != LINV) {
        // mfh 21 Dec 2015: Here we encode the assumption that blocks
        // are stored contiguously, with no padding.  "Contiguously"
        // means that all memory between the first and last entries
        // belongs to the block (no striding).  "No padding" means
        // that getBlockSize() * getBlockSize() is exactly the number
        // of entries that the block uses.  For another place where
        // this assumption is encoded, see sumIntoLocalValues.

        const size_t absBlockOffset = absRowBlockOffset + relBlockOffset;
        // little_block_type A_old =
        //   getNonConstLocalBlockFromAbsOffset (absBlockOffset);
        impl_scalar_type* const A_old =
          vals_host_out_raw + absBlockOffset * perBlockSize;
        // const_little_block_type A_new =
        //   getConstLocalBlockFromInput (vIn, pointOffset);
        const impl_scalar_type* const A_new = vIn + pointOffset;
        // COPY (A_new, A_old);
        for (LO i = 0; i < perBlockSize; ++i) {
          A_old[i] = A_new[i];
        }
        hint = relBlockOffset + 1;
        ++validCount;
      }
    }
    return validCount;
  }

  template <class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  getLocalDiagOffsets (const Kokkos::View<size_t*, device_type,
                         Kokkos::MemoryUnmanaged>& offsets) const
  {
    graph_.getLocalDiagOffsets (offsets);
  }

  template <class Scalar, class LO, class GO, class Node>
  void TPETRA_DEPRECATED
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  getLocalDiagOffsets (Teuchos::ArrayRCP<size_t>& offsets) const
  {
    // mfh 19 Mar 2016: We plan to deprecate the ArrayRCP version of
    // this method in CrsGraph too, so don't call it (otherwise build
    // warnings will show up and annoy users).  Instead, copy results
    // in and out, if the memory space requires it.

    const size_t lclNumRows = graph_.getNodeNumRows ();
    if (static_cast<size_t> (offsets.size ()) < lclNumRows) {
      offsets.resize (lclNumRows);
    }

    // The input ArrayRCP must always be a host pointer.  Thus, if
    // device_type::memory_space is Kokkos::HostSpace, it's OK for us
    // to write to that allocation directly as a Kokkos::View.
    typedef typename device_type::memory_space memory_space;
    if (std::is_same<memory_space, Kokkos::HostSpace>::value) {
      // It is always syntactically correct to assign a raw host
      // pointer to a device View, so this code will compile correctly
      // even if this branch never runs.
      typedef Kokkos::View<size_t*, device_type,
                           Kokkos::MemoryUnmanaged> output_type;
      output_type offsetsOut (offsets.getRawPtr (), lclNumRows);
      graph_.getLocalDiagOffsets (offsetsOut);
    }
    else {
      Kokkos::View<size_t*, device_type> offsetsTmp ("diagOffsets", lclNumRows);
      graph_.getLocalDiagOffsets (offsetsTmp);
      typedef Kokkos::View<size_t*, Kokkos::HostSpace,
                           Kokkos::MemoryUnmanaged> output_type;
      output_type offsetsOut (offsets.getRawPtr (), lclNumRows);
      Kokkos::deep_copy (offsetsOut, offsetsTmp);
    }
  }

  template <class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  localGaussSeidel (const BlockMultiVector<Scalar, LO, GO, Node>& B,
                    BlockMultiVector<Scalar, LO, GO, Node>& X,
                    const Kokkos::View<impl_scalar_type***, device_type,
                      Kokkos::MemoryUnmanaged>& D_inv,
                    const Scalar& omega,
                    const ESweepDirection direction) const
  {
    using Kokkos::ALL;
    const impl_scalar_type zero =
      Kokkos::Details::ArithTraits<impl_scalar_type>::zero ();
    const impl_scalar_type one =
      Kokkos::Details::ArithTraits<impl_scalar_type>::one ();
    const LO numLocalMeshRows =
      static_cast<LO> (rowMeshMap_.getNodeNumElements ());
    const LO numVecs = static_cast<LO> (X.getNumVectors ());

    // If using (new) Kokkos, replace localMem with thread-local
    // memory.  Note that for larger block sizes, this will affect the
    // two-level parallelization.  Look to Stokhos for best practice
    // on making this fast for GPUs.
    const LO blockSize = getBlockSize ();
    Teuchos::Array<impl_scalar_type> localMem (blockSize);
    Teuchos::Array<impl_scalar_type> localMat (blockSize*blockSize);
    little_vec_type X_lcl (localMem.getRawPtr (), blockSize);

    // FIXME (mfh 12 Aug 2014) This probably won't work if LO is unsigned.
    LO rowBegin = 0, rowEnd = 0, rowStride = 0;
    if (direction == Forward) {
      rowBegin = 1;
      rowEnd = numLocalMeshRows+1;
      rowStride = 1;
    }
    else if (direction == Backward) {
      rowBegin = numLocalMeshRows;
      rowEnd = 0;
      rowStride = -1;
    }
    else if (direction == Symmetric) {
      this->localGaussSeidel (B, X, D_inv, omega, Forward);
      this->localGaussSeidel (B, X, D_inv, omega, Backward);
      return;
    }

    const Scalar one_minus_omega = Teuchos::ScalarTraits<Scalar>::one()-omega;
    const Scalar     minus_omega = -omega;

    if (numVecs == 1) {
      for (LO lclRow = rowBegin; lclRow != rowEnd; lclRow += rowStride) {
        const LO actlRow = lclRow - 1;

        little_vec_type B_cur = B.getLocalBlock (actlRow, 0);
        COPY (B_cur, X_lcl);
        SCAL (static_cast<impl_scalar_type> (omega), X_lcl);

        const size_t meshBeg = ptrHost_[actlRow];
        const size_t meshEnd = ptrHost_[actlRow+1];
        for (size_t absBlkOff = meshBeg; absBlkOff < meshEnd; ++absBlkOff) {
          const LO meshCol = indHost_[absBlkOff];
          const_little_block_type A_cur =
            getConstLocalBlockFromAbsOffset (absBlkOff);
          little_vec_type X_cur = X.getLocalBlock (meshCol, 0);

          // X_lcl += alpha*A_cur*X_cur
          const Scalar alpha = meshCol == actlRow ? one_minus_omega : minus_omega;
          //X_lcl.matvecUpdate (alpha, A_cur, X_cur);
          GEMV (static_cast<impl_scalar_type> (alpha), A_cur, X_cur, X_lcl);
        } // for each entry in the current local row of the matrix

        // NOTE (mfh 20 Jan 2016) The two input Views here are
        // unmanaged already, so we don't have to take unmanaged
        // subviews first.
        auto D_lcl = Kokkos::subview (D_inv, actlRow, ALL (), ALL ());
        little_vec_type X_update = X.getLocalBlock (actlRow, 0);
        FILL (X_update, zero);
        GEMV (one, D_lcl, X_lcl, X_update); // overwrite X_update
      } // for each local row of the matrix
    }
    else {
      for (LO lclRow = rowBegin; lclRow != rowEnd; lclRow += rowStride) {
        for (LO j = 0; j < numVecs; ++j) {
          LO actlRow = lclRow-1;

          little_vec_type B_cur = B.getLocalBlock (actlRow, j);
          COPY (B_cur, X_lcl);
          SCAL (static_cast<impl_scalar_type> (omega), X_lcl);

          const size_t meshBeg = ptrHost_[actlRow];
          const size_t meshEnd = ptrHost_[actlRow+1];
          for (size_t absBlkOff = meshBeg; absBlkOff < meshEnd; ++absBlkOff) {
            const LO meshCol = indHost_[absBlkOff];
            const_little_block_type A_cur =
              getConstLocalBlockFromAbsOffset (absBlkOff);
            little_vec_type X_cur = X.getLocalBlock (meshCol, j);

            // X_lcl += alpha*A_cur*X_cur
            const Scalar alpha = meshCol == actlRow ? one_minus_omega : minus_omega;
            GEMV (static_cast<impl_scalar_type> (alpha), A_cur, X_cur, X_lcl);
          } // for each entry in the current local row of the matrx

          auto D_lcl = Kokkos::subview (D_inv, actlRow, ALL (), ALL ());
          auto X_update = X.getLocalBlock (actlRow, j);
          FILL (X_update, zero);
          GEMV (one, D_lcl, X_lcl, X_update); // overwrite X_update
        } // for each entry in the current local row of the matrix
      } // for each local row of the matrix
    }
  }

  template <class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  gaussSeidelCopy (MultiVector<Scalar,LO,GO,Node> &X,
                   const MultiVector<Scalar,LO,GO,Node> &B,
                   const MultiVector<Scalar,LO,GO,Node> &D,
                   const Scalar& dampingFactor,
                   const ESweepDirection direction,
                   const int numSweeps,
                   const bool zeroInitialGuess) const
  {
    // FIXME (mfh 12 Aug 2014) This method has entirely the wrong
    // interface for block Gauss-Seidel.
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
      "gaussSeidelCopy: Not implemented.");
  }

  template <class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  reorderedGaussSeidelCopy (MultiVector<Scalar,LO,GO,Node>& X,
                            const MultiVector<Scalar,LO,GO,Node>& B,
                            const MultiVector<Scalar,LO,GO,Node>& D,
                            const Teuchos::ArrayView<LO>& rowIndices,
                            const Scalar& dampingFactor,
                            const ESweepDirection direction,
                            const int numSweeps,
                            const bool zeroInitialGuess) const
  {
    // FIXME (mfh 12 Aug 2014) This method has entirely the wrong
    // interface for block Gauss-Seidel.
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
      "reorderedGaussSeidelCopy: Not implemented.");
  }

  template <class Scalar, class LO, class GO, class Node>
  void TPETRA_DEPRECATED
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  getLocalDiagCopy (BlockCrsMatrix<Scalar,LO,GO,Node>& diag,
                    const Teuchos::ArrayView<const size_t>& offsets) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero ();

    const size_t myNumRows = rowMeshMap_.getNodeNumElements();
    const LO* columnIndices;
    Scalar* vals;
    LO numColumns;
    Teuchos::Array<LO> cols(1);

    // FIXME (mfh 12 Aug 2014) Should use a "little block" for this instead.
    Teuchos::Array<Scalar> zeroMat (blockSize_*blockSize_, ZERO);
    for (size_t i = 0; i < myNumRows; ++i) {
      cols[0] = i;
      if (offsets[i] == Teuchos::OrdinalTraits<size_t>::invalid ()) {
        diag.replaceLocalValues (i, cols.getRawPtr (), zeroMat.getRawPtr (), 1);
      }
      else {
        getLocalRowView (i, columnIndices, vals, numColumns);
        diag.replaceLocalValues (i, cols.getRawPtr(), &vals[offsets[i]*blockSize_*blockSize_], 1);
      }
    }
  }


  template <class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  getLocalDiagCopy (const Kokkos::View<impl_scalar_type***, device_type,
                                       Kokkos::MemoryUnmanaged>& diag,
                    const Kokkos::View<const size_t*, device_type,
                                       Kokkos::MemoryUnmanaged>& offsets) const
  {
    using Kokkos::parallel_for;
    typedef typename device_type::execution_space execution_space;
    const char prefix[] = "Tpetra::BlockCrsMatrix::getLocalDiagCopy (2-arg): ";

    const LO lclNumMeshRows = static_cast<LO> (rowMeshMap_.getNodeNumElements ());
    const LO blockSize = this->getBlockSize ();
    TEUCHOS_TEST_FOR_EXCEPTION
      (static_cast<LO> (diag.extent (0)) < lclNumMeshRows ||
       static_cast<LO> (diag.extent (1)) < blockSize ||
       static_cast<LO> (diag.extent (2)) < blockSize,
       std::invalid_argument, prefix <<
       "The input Kokkos::View is not big enough to hold all the data.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (static_cast<LO> (offsets.size ()) < lclNumMeshRows, std::invalid_argument,
       prefix << "offsets.size() = " << offsets.size () << " < local number of "
       "diagonal blocks " << lclNumMeshRows << ".");

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION
      (this->template need_sync<device_type> (), std::runtime_error,
       prefix << "The matrix's data were last modified on host, but have "
       "not been sync'd to device.  Please sync to device (by calling "
       "sync<device_type>() on this matrix) before calling this method.");
#endif // HAVE_TPETRA_DEBUG

    typedef Kokkos::RangePolicy<execution_space, LO> policy_type;
    typedef GetLocalDiagCopy<Scalar, LO, GO, Node> functor_type;

    // FIXME (mfh 26 May 2016) Not really OK to const_cast here, since
    // we reserve the right to do lazy allocation of device data.  (We
    // don't plan to do lazy allocation for host data; the host
    // version of the data always exists.)
    typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_type;
    auto vals_dev =
      const_cast<this_type*> (this)->template getValues<device_type> ();

    parallel_for (policy_type (0, lclNumMeshRows),
                  functor_type (diag, vals_dev, offsets,
                                graph_.getLocalGraph ().row_map, blockSize_));
  }


  template <class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  getLocalDiagCopy (const Kokkos::View<impl_scalar_type***, device_type,
                                       Kokkos::MemoryUnmanaged>& diag,
                    const Teuchos::ArrayView<const size_t>& offsets) const
  {
    using Kokkos::ALL;
    using Kokkos::parallel_for;
    typedef typename Kokkos::View<impl_scalar_type***, device_type,
      Kokkos::MemoryUnmanaged>::HostMirror::execution_space host_exec_space;

    const LO lclNumMeshRows = static_cast<LO> (rowMeshMap_.getNodeNumElements ());
    const LO blockSize = this->getBlockSize ();
    TEUCHOS_TEST_FOR_EXCEPTION
      (static_cast<LO> (diag.extent (0)) < lclNumMeshRows ||
       static_cast<LO> (diag.extent (1)) < blockSize ||
       static_cast<LO> (diag.extent (2)) < blockSize,
       std::invalid_argument, "Tpetra::BlockCrsMatrix::getLocalDiagCopy: "
       "The input Kokkos::View is not big enough to hold all the data.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (static_cast<LO> (offsets.size ()) < lclNumMeshRows,
       std::invalid_argument, "Tpetra::BlockCrsMatrix::getLocalDiagCopy: "
       "offsets.size() = " << offsets.size () << " < local number of diagonal "
       "blocks " << lclNumMeshRows << ".");

    // mfh 12 Dec 2015: Use the host execution space, since we haven't
    // quite made everything work with CUDA yet.
    typedef Kokkos::RangePolicy<host_exec_space, LO> policy_type;
    parallel_for (policy_type (0, lclNumMeshRows), [=] (const LO& lclMeshRow) {
        auto D_in = this->getConstLocalBlockFromRelOffset (lclMeshRow, offsets[lclMeshRow]);
        auto D_out = Kokkos::subview (diag, lclMeshRow, ALL (), ALL ());
        COPY (D_in, D_out);
      });
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  absMaxLocalValues (const LO localRowInd,
                     const LO colInds[],
                     const Scalar vals[],
                     const LO numColInds) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We modified no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }
    const impl_scalar_type* const vIn =
      reinterpret_cast<const impl_scalar_type*> (vals);
    const size_t absRowBlockOffset = ptrHost_[localRowInd];
    const LO LINV = Teuchos::OrdinalTraits<LO>::invalid ();
    const LO perBlockSize = this->offsetPerBlock ();
    LO hint = 0; // Guess for the relative offset into the current row
    LO pointOffset = 0; // Current offset into input values
    LO validCount = 0; // number of valid column indices in colInds

    for (LO k = 0; k < numColInds; ++k, pointOffset += perBlockSize) {
      const LO relBlockOffset =
        this->findRelOffsetOfColumnIndex (localRowInd, colInds[k], hint);
      if (relBlockOffset != LINV) {
        const size_t absBlockOffset = absRowBlockOffset + relBlockOffset;
        little_block_type A_old =
          getNonConstLocalBlockFromAbsOffset (absBlockOffset);
        const_little_block_type A_new =
          getConstLocalBlockFromInput (vIn, pointOffset);

        ::Tpetra::Experimental::Impl::absMax (A_old, A_new);
        hint = relBlockOffset + 1;
        ++validCount;
      }
    }
    return validCount;
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  sumIntoLocalValues (const LO localRowInd,
                      const LO colInds[],
                      const Scalar vals[],
                      const LO numColInds) const
  {
#ifdef HAVE_TPETRA_DEBUG
    const char prefix[] =
      "Tpetra::Experimental::BlockCrsMatrix::sumIntoLocalValues: ";
#endif // HAVE_TPETRA_DEBUG

    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We modified no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }
    //const impl_scalar_type ONE = static_cast<impl_scalar_type> (1.0);
    const impl_scalar_type* const vIn =
      reinterpret_cast<const impl_scalar_type*> (vals);
    const size_t absRowBlockOffset = ptrHost_[localRowInd];
    const LO LINV = Teuchos::OrdinalTraits<LO>::invalid ();
    const LO perBlockSize = this->offsetPerBlock ();
    LO hint = 0; // Guess for the relative offset into the current row
    LO pointOffset = 0; // Current offset into input values
    LO validCount = 0; // number of valid column indices in colInds

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION
      (this->need_sync_host (), std::runtime_error,
       prefix << "The matrix's data were last modified on device, but have not "
       "been sync'd to host.  Please sync to host (by calling "
       "sync<Kokkos::HostSpace>() on this matrix) before calling this method.");
#endif // HAVE_TPETRA_DEBUG

    auto vals_host_out =
      getValuesHost ();
    impl_scalar_type* vals_host_out_raw = vals_host_out.data ();

    for (LO k = 0; k < numColInds; ++k, pointOffset += perBlockSize) {
      const LO relBlockOffset =
        this->findRelOffsetOfColumnIndex (localRowInd, colInds[k], hint);
      if (relBlockOffset != LINV) {
        // mfh 21 Dec 2015: Here we encode the assumption that blocks
        // are stored contiguously, with no padding.  "Contiguously"
        // means that all memory between the first and last entries
        // belongs to the block (no striding).  "No padding" means
        // that getBlockSize() * getBlockSize() is exactly the number
        // of entries that the block uses.  For another place where
        // this assumption is encoded, see replaceLocalValues.

        const size_t absBlockOffset = absRowBlockOffset + relBlockOffset;
        // little_block_type A_old =
        //   getNonConstLocalBlockFromAbsOffset (absBlockOffset);
        impl_scalar_type* const A_old =
          vals_host_out_raw + absBlockOffset * perBlockSize;
        // const_little_block_type A_new =
        //   getConstLocalBlockFromInput (vIn, pointOffset);
        const impl_scalar_type* const A_new = vIn + pointOffset;
        // AXPY (ONE, A_new, A_old);
        for (LO i = 0; i < perBlockSize; ++i) {
          A_old[i] += A_new[i];
        }
        hint = relBlockOffset + 1;
        ++validCount;
      }
    }
    return validCount;
  }

  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalRowView (const LO localRowInd,
                   const LO*& colInds,
                   Scalar*& vals,
                   LO& numInds) const
  {
#ifdef HAVE_TPETRA_DEBUG
    const char prefix[] =
      "Tpetra::Experimental::BlockCrsMatrix::getLocalRowView: ";
#endif // HAVE_TPETRA_DEBUG

    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      colInds = NULL;
      vals = NULL;
      numInds = 0;
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    else {
      const size_t absBlockOffsetStart = ptrHost_[localRowInd];
      colInds = indHost_.data () + absBlockOffsetStart;

#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION
        (this->need_sync_host (), std::runtime_error,
         prefix << "The matrix's data were last modified on device, but have "
         "not been sync'd to host.  Please sync to host (by calling "
         "sync<Kokkos::HostSpace>() on this matrix) before calling this "
         "method.");
#endif // HAVE_TPETRA_DEBUG

      auto vals_host_out = getValuesHost ();
      impl_scalar_type* vals_host_out_raw = vals_host_out.data ();
      impl_scalar_type* const vOut = vals_host_out_raw +
        absBlockOffsetStart * offsetPerBlock ();
      vals = reinterpret_cast<Scalar*> (vOut);

      numInds = ptrHost_[localRowInd + 1] - absBlockOffsetStart;
      return 0; // indicates no error
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalRowCopy (LO LocalRow,
                   const Teuchos::ArrayView<LO>& Indices,
                   const Teuchos::ArrayView<Scalar>& Values,
                   size_t &NumEntries) const
  {
    const LO *colInds;
    Scalar *vals;
    LO numInds;
    getLocalRowView(LocalRow,colInds,vals,numInds);
    if (numInds > Indices.size() || numInds*blockSize_*blockSize_ > Values.size()) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                  "Tpetra::BlockCrsMatrix::getLocalRowCopy : Column and/or values array is not large enough to hold "
                  << numInds << " row entries");
    }
    for (LO i=0; i<numInds; ++i) {
      Indices[i] = colInds[i];
    }
    for (LO i=0; i<numInds*blockSize_*blockSize_; ++i) {
      Values[i] = vals[i];
    }
    NumEntries = numInds;
  }

  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalRowOffsets (const LO localRowInd,
                      ptrdiff_t offsets[],
                      const LO colInds[],
                      const LO numColInds) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We got no offsets, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }

    const LO LINV = Teuchos::OrdinalTraits<LO>::invalid ();
    LO hint = 0; // Guess for the relative offset into the current row
    LO validCount = 0; // number of valid column indices in colInds

    for (LO k = 0; k < numColInds; ++k) {
      const LO relBlockOffset =
        this->findRelOffsetOfColumnIndex (localRowInd, colInds[k], hint);
      if (relBlockOffset != LINV) {
        offsets[k] = static_cast<ptrdiff_t> (relBlockOffset);
        hint = relBlockOffset + 1;
        ++validCount;
      }
    }
    return validCount;
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  replaceLocalValuesByOffsets (const LO localRowInd,
                               const ptrdiff_t offsets[],
                               const Scalar vals[],
                               const LO numOffsets) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We modified no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }
    const impl_scalar_type* const vIn = reinterpret_cast<const impl_scalar_type*> (vals);

    const size_t absRowBlockOffset = ptrHost_[localRowInd];
    const size_t perBlockSize = static_cast<LO> (offsetPerBlock ());
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    size_t pointOffset = 0; // Current offset into input values
    LO validCount = 0; // number of valid offsets

    for (LO k = 0; k < numOffsets; ++k, pointOffset += perBlockSize) {
      const size_t relBlockOffset = offsets[k];
      if (relBlockOffset != STINV) {
        const size_t absBlockOffset = absRowBlockOffset + relBlockOffset;
        little_block_type A_old =
          getNonConstLocalBlockFromAbsOffset (absBlockOffset);
        const_little_block_type A_new =
          getConstLocalBlockFromInput (vIn, pointOffset);
        COPY (A_new, A_old);
        ++validCount;
      }
    }
    return validCount;
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  absMaxLocalValuesByOffsets (const LO localRowInd,
                              const ptrdiff_t offsets[],
                              const Scalar vals[],
                              const LO numOffsets) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We modified no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }
    const impl_scalar_type* const vIn = reinterpret_cast<const impl_scalar_type*> (vals);

    const size_t absRowBlockOffset = ptrHost_[localRowInd];
    const size_t perBlockSize = static_cast<LO> (offsetPerBlock ());
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    size_t pointOffset = 0; // Current offset into input values
    LO validCount = 0; // number of valid offsets

    for (LO k = 0; k < numOffsets; ++k, pointOffset += perBlockSize) {
      const size_t relBlockOffset = offsets[k];
      if (relBlockOffset != STINV) {
        const size_t absBlockOffset = absRowBlockOffset + relBlockOffset;
        little_block_type A_old =
          getNonConstLocalBlockFromAbsOffset (absBlockOffset);
        const_little_block_type A_new =
          getConstLocalBlockFromInput (vIn, pointOffset);
        ::Tpetra::Experimental::Impl::absMax (A_old, A_new);
        ++validCount;
      }
    }
    return validCount;
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  sumIntoLocalValuesByOffsets (const LO localRowInd,
                               const ptrdiff_t offsets[],
                               const Scalar vals[],
                               const LO numOffsets) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We modified no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }
    const impl_scalar_type ONE = static_cast<impl_scalar_type> (1.0);
    const impl_scalar_type* const vIn = reinterpret_cast<const impl_scalar_type*> (vals);

    const size_t absRowBlockOffset = ptrHost_[localRowInd];
    const size_t perBlockSize = static_cast<LO> (offsetPerBlock ());
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    size_t pointOffset = 0; // Current offset into input values
    LO validCount = 0; // number of valid offsets

    for (LO k = 0; k < numOffsets; ++k, pointOffset += perBlockSize) {
      const size_t relBlockOffset = offsets[k];
      if (relBlockOffset != STINV) {
        const size_t absBlockOffset = absRowBlockOffset + relBlockOffset;
        little_block_type A_old =
          getNonConstLocalBlockFromAbsOffset (absBlockOffset);
        const_little_block_type A_new =
          getConstLocalBlockFromInput (vIn, pointOffset);
        //A_old.update (ONE, A_new);
        AXPY (ONE, A_new, A_old);
        ++validCount;
      }
    }
    return validCount;
  }


  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNumEntriesInLocalRow (const LO localRowInd) const
  {
    const size_t numEntInGraph = graph_.getNumEntriesInLocalRow (localRowInd);
    if (numEntInGraph == Teuchos::OrdinalTraits<size_t>::invalid ()) {
      return static_cast<LO> (0); // the calling process doesn't have that row
    } else {
      return static_cast<LO> (numEntInGraph);
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  applyBlockTrans (const BlockMultiVector<Scalar, LO, GO, Node>& X,
                   BlockMultiVector<Scalar, LO, GO, Node>& Y,
                   const Teuchos::ETransp mode,
                   const Scalar alpha,
                   const Scalar beta)
  {
    (void) X;
    (void) Y;
    (void) mode;
    (void) alpha;
    (void) beta;

    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::apply: "
      "transpose and conjugate transpose modes are not yet implemented.");
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  applyBlockNoTrans (const BlockMultiVector<Scalar, LO, GO, Node>& X,
                     BlockMultiVector<Scalar, LO, GO, Node>& Y,
                     const Scalar alpha,
                     const Scalar beta)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef ::Tpetra::Import<LO, GO, Node> import_type;
    typedef ::Tpetra::Export<LO, GO, Node> export_type;
    const Scalar zero = STS::zero ();
    const Scalar one = STS::one ();
    RCP<const import_type> import = graph_.getImporter ();
    // "export" is a reserved C++ keyword, so we can't use it.
    RCP<const export_type> theExport = graph_.getExporter ();
    const char prefix[] = "Tpetra::BlockCrsMatrix::applyBlockNoTrans: ";

    if (alpha == zero) {
      if (beta == zero) {
        Y.putScalar (zero); // replace Inf or NaN (BLAS rules)
      }
      else if (beta != one) {
        Y.scale (beta);
      }
    }
    else { // alpha != 0
      const BMV* X_colMap = NULL;
      if (import.is_null ()) {
        try {
          X_colMap = &X;
        }
        catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::logic_error, prefix << "Tpetra::MultiVector::"
             "operator= threw an exception: " << std::endl << e.what ());
        }
      }
      else {
        // X_colMap_ is a pointer to a pointer to BMV.  Ditto for
        // Y_rowMap_ below.  This lets us do lazy initialization
        // correctly with view semantics of BlockCrsMatrix.  All views
        // of this BlockCrsMatrix have the same outer pointer.  That
        // way, we can set the inner pointer in one view, and all
        // other views will see it.
        if ((*X_colMap_).is_null () ||
            (**X_colMap_).getNumVectors () != X.getNumVectors () ||
            (**X_colMap_).getBlockSize () != X.getBlockSize ()) {
          *X_colMap_ = rcp (new BMV (* (graph_.getColMap ()), getBlockSize (),
                                     static_cast<LO> (X.getNumVectors ())));
        }
#ifdef HAVE_TPETRA_BCRS_DO_POINT_IMPORT
        if (pointImporter_->is_null ()) {
          // The Import ctor needs RCPs. Map's copy ctor does a shallow copy, so
          // these are small operations.
          const auto domainPointMap = rcp (new typename BMV::map_type (domainPointMap_));
          const auto colPointMap = rcp (new typename BMV::map_type (
                                          BMV::makePointMap (*graph_.getColMap(),
                                                             blockSize_)));
          *pointImporter_ = rcp (new typename crs_graph_type::import_type (
                                   domainPointMap, colPointMap));
        }
        (*X_colMap_)->getMultiVectorView().doImport (X.getMultiVectorView (),
                                                     **pointImporter_,
                                                     ::Tpetra::REPLACE);
#else
        (**X_colMap_).doImport (X, *import, ::Tpetra::REPLACE);
#endif
        try {
          X_colMap = &(**X_colMap_);
        }
        catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::logic_error, prefix << "Tpetra::MultiVector::"
            "operator= threw an exception: " << std::endl << e.what ());
        }
      }

      BMV* Y_rowMap = NULL;
      if (theExport.is_null ()) {
        Y_rowMap = &Y;
      }
      else if ((*Y_rowMap_).is_null () ||
                 (**Y_rowMap_).getNumVectors () != Y.getNumVectors () ||
                 (**Y_rowMap_).getBlockSize () != Y.getBlockSize ()) {
        *Y_rowMap_ = rcp (new BMV (* (graph_.getRowMap ()), getBlockSize (),
                                   static_cast<LO> (X.getNumVectors ())));
        try {
          Y_rowMap = &(**Y_rowMap_);
        }
        catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::logic_error, prefix << "Tpetra::MultiVector::"
            "operator= threw an exception: " << std::endl << e.what ());
        }
      }

      try {
        localApplyBlockNoTrans (*X_colMap, *Y_rowMap, alpha, beta);
      }
      catch (std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::runtime_error, prefix << "localApplyBlockNoTrans threw "
           "an exception: " << e.what ());
      }
      catch (...) {
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::runtime_error, prefix << "localApplyBlockNoTrans threw "
           "an exception not a subclass of std::exception.");
      }

      if (! theExport.is_null ()) {
        Y.doExport (*Y_rowMap, *theExport, ::Tpetra::REPLACE);
      }
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  localApplyBlockNoTrans (const BlockMultiVector<Scalar, LO, GO, Node>& X,
                          BlockMultiVector<Scalar, LO, GO, Node>& Y,
                          const Scalar alpha,
                          const Scalar beta)
  {
    using ::Tpetra::Experimental::Impl::bcrsLocalApplyNoTrans;

    const impl_scalar_type alpha_impl = alpha;
    const auto graph = this->graph_.getLocalGraph ();
    const impl_scalar_type beta_impl = beta;
    const LO blockSize = this->getBlockSize ();

    auto X_mv = X.getMultiVectorView ();
    auto Y_mv = Y.getMultiVectorView ();
    Y_mv.template modify<device_type> ();

    auto X_lcl = X_mv.template getLocalView<device_type> ();
    auto Y_lcl = Y_mv.template getLocalView<device_type> ();
    auto val = this->val_.template view<device_type> ();

    bcrsLocalApplyNoTrans (alpha_impl, graph, val, blockSize, X_lcl,
                           beta_impl, Y_lcl);
  }

  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  findRelOffsetOfColumnIndex (const LO localRowIndex,
                              const LO colIndexToFind,
                              const LO hint) const
  {
    const size_t absStartOffset = ptrHost_[localRowIndex];
    const size_t absEndOffset = ptrHost_[localRowIndex+1];
    const LO numEntriesInRow = static_cast<LO> (absEndOffset - absStartOffset);
    // Amortize pointer arithmetic over the search loop.
    const LO* const curInd = indHost_.data () + absStartOffset;

    // If the hint was correct, then the hint is the offset to return.
    if (hint < numEntriesInRow && curInd[hint] == colIndexToFind) {
      // Always return the offset relative to the current row.
      return hint;
    }

    // The hint was wrong, so we must search for the given column
    // index in the column indices for the given row.
    LO relOffset = Teuchos::OrdinalTraits<LO>::invalid ();

    // We require that the graph have sorted rows.  However, binary
    // search only pays if the current row is longer than a certain
    // amount.  We set this to 32, but you might want to tune this.
    const LO maxNumEntriesForLinearSearch = 32;
    if (numEntriesInRow > maxNumEntriesForLinearSearch) {
      // Use binary search.  It would probably be better for us to
      // roll this loop by hand.  If we wrote it right, a smart
      // compiler could perhaps use conditional loads and avoid
      // branches (according to Jed Brown on May 2014).
      const LO* beg = curInd;
      const LO* end = curInd + numEntriesInRow;
      std::pair<const LO*, const LO*> p =
        std::equal_range (beg, end, colIndexToFind);
      if (p.first != p.second) {
        // offset is relative to the current row
        relOffset = static_cast<LO> (p.first - beg);
      }
    }
    else { // use linear search
      for (LO k = 0; k < numEntriesInRow; ++k) {
        if (colIndexToFind == curInd[k]) {
          relOffset = k; // offset is relative to the current row
          break;
        }
      }
    }

    return relOffset;
  }

  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  offsetPerBlock () const
  {
    return offsetPerBlock_;
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::const_little_block_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getConstLocalBlockFromInput (const impl_scalar_type* val,
                               const size_t pointOffset) const
  {
    // Row major blocks
    const LO rowStride = blockSize_;
    return const_little_block_type (val + pointOffset, blockSize_, rowStride);
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::little_block_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNonConstLocalBlockFromInput (impl_scalar_type* val,
                                  const size_t pointOffset) const
  {
    // Row major blocks
    const LO rowStride = blockSize_;
    return little_block_type (val + pointOffset, blockSize_, rowStride);
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::const_little_block_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getConstLocalBlockFromAbsOffset (const size_t absBlockOffset) const
  {
#ifdef HAVE_TPETRA_DEBUG
    const char prefix[] =
      "Tpetra::Experimental::BlockCrsMatrix::getConstLocalBlockFromAbsOffset: ";
#endif // HAVE_TPETRA_DEBUG

    if (absBlockOffset >= ptrHost_[rowMeshMap_.getNodeNumElements ()]) {
      // An empty block signifies an error.  We don't expect to see
      // this error in correct code, but it's helpful for avoiding
      // memory corruption in case there is a bug.
      return const_little_block_type ();
    }
    else {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION
        (this->need_sync_host (), std::runtime_error,
         prefix << "The matrix's data were last modified on device, but have "
         "not been sync'd to host.  Please sync to host (by calling "
         "sync<Kokkos::HostSpace>() on this matrix) before calling this "
         "method.");
#endif // HAVE_TPETRA_DEBUG
      const size_t absPointOffset = absBlockOffset * offsetPerBlock ();

      auto vals_host = getValuesHost ();
      const impl_scalar_type* vals_host_raw = vals_host.data ();

      return getConstLocalBlockFromInput (vals_host_raw, absPointOffset);
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::const_little_block_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getConstLocalBlockFromRelOffset (const LO lclMeshRow,
                                   const size_t relMeshOffset) const
  {
    typedef impl_scalar_type IST;

    const LO* lclColInds = NULL;
    Scalar* lclVals = NULL;
    LO numEnt = 0;

    LO err = this->getLocalRowView (lclMeshRow, lclColInds, lclVals, numEnt);
    if (err != 0) {
      // An empty block signifies an error.  We don't expect to see
      // this error in correct code, but it's helpful for avoiding
      // memory corruption in case there is a bug.
      return const_little_block_type ();
    }
    else {
      const size_t relPointOffset = relMeshOffset * this->offsetPerBlock ();
      IST* lclValsImpl = reinterpret_cast<IST*> (lclVals);
      return this->getConstLocalBlockFromInput (const_cast<const IST*> (lclValsImpl),
                                                relPointOffset);
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::little_block_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNonConstLocalBlockFromAbsOffset (const size_t absBlockOffset) const
  {
#ifdef HAVE_TPETRA_DEBUG
    const char prefix[] =
      "Tpetra::Experimental::BlockCrsMatrix::getNonConstLocalBlockFromAbsOffset: ";
#endif // HAVE_TPETRA_DEBUG

    if (absBlockOffset >= ptrHost_[rowMeshMap_.getNodeNumElements ()]) {
      // An empty block signifies an error.  We don't expect to see
      // this error in correct code, but it's helpful for avoiding
      // memory corruption in case there is a bug.
      return little_block_type ();
    }
    else {
      const size_t absPointOffset = absBlockOffset * offsetPerBlock ();
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION
        (this->need_sync_host (), std::runtime_error,
         prefix << "The matrix's data were last modified on device, but have "
         "not been sync'd to host.  Please sync to host (by calling "
         "sync<Kokkos::HostSpace>() on this matrix) before calling this "
         "method.");
#endif // HAVE_TPETRA_DEBUG
      auto vals_host = getValuesHost();
      impl_scalar_type* vals_host_raw = vals_host.data ();
      return getNonConstLocalBlockFromInput (vals_host_raw, absPointOffset);
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::little_block_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalBlock (const LO localRowInd, const LO localColInd) const
  {
    const size_t absRowBlockOffset = ptrHost_[localRowInd];
    const LO relBlockOffset =
      this->findRelOffsetOfColumnIndex (localRowInd, localColInd);

    if (relBlockOffset != Teuchos::OrdinalTraits<LO>::invalid ()) {
      const size_t absBlockOffset = absRowBlockOffset + relBlockOffset;
      return getNonConstLocalBlockFromAbsOffset (absBlockOffset);
    }
    else {
      return little_block_type ();
    }
  }

  // template<class Scalar, class LO, class GO, class Node>
  // void
  // BlockCrsMatrix<Scalar, LO, GO, Node>::
  // clearLocalErrorStateAndStream ()
  // {
  //   typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_type;
  //   * (const_cast<this_type*> (this)->localError_) = false;
  //   *errs_ = Teuchos::null;
  // }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  checkSizes (const ::Tpetra::SrcDistObject& source)
  {
    using std::endl;
    typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_type;
    const this_type* src = dynamic_cast<const this_type* > (&source);

    if (src == NULL) {
      std::ostream& err = markLocalErrorAndGetStream ();
      err << "checkSizes: The source object of the Import or Export "
        "must be a BlockCrsMatrix with the same template parameters as the "
        "target object." << endl;
    }
    else {
      // Use a string of ifs, not if-elseifs, because we want to know
      // all the errors.
      if (src->getBlockSize () != this->getBlockSize ()) {
        std::ostream& err = markLocalErrorAndGetStream ();
        err << "checkSizes: The source and target objects of the Import or "
            << "Export must have the same block sizes.  The source's block "
            << "size = " << src->getBlockSize () << " != the target's block "
            << "size = " << this->getBlockSize () << "." << endl;
      }
      if (! src->graph_.isFillComplete ()) {
        std::ostream& err = markLocalErrorAndGetStream ();
        err << "checkSizes: The source object of the Import or Export is "
          "not fill complete.  Both source and target objects must be fill "
          "complete." << endl;
      }
      if (! this->graph_.isFillComplete ()) {
        std::ostream& err = markLocalErrorAndGetStream ();
        err << "checkSizes: The target object of the Import or Export is "
          "not fill complete.  Both source and target objects must be fill "
          "complete." << endl;
      }
      if (src->graph_.getColMap ().is_null ()) {
        std::ostream& err = markLocalErrorAndGetStream ();
        err << "checkSizes: The source object of the Import or Export does "
          "not have a column Map.  Both source and target objects must have "
          "column Maps." << endl;
      }
      if (this->graph_.getColMap ().is_null ()) {
        std::ostream& err = markLocalErrorAndGetStream ();
        err << "checkSizes: The target object of the Import or Export does "
          "not have a column Map.  Both source and target objects must have "
          "column Maps." << endl;
      }
    }

    return ! (* (this->localError_));
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  copyAndPermuteNew (const ::Tpetra::SrcDistObject& source,
                     const size_t numSameIDs,
                     const Kokkos::DualView<const local_ordinal_type*, device_type>& permuteToLIDs,
                     const Kokkos::DualView<const local_ordinal_type*, device_type>& permuteFromLIDs)
  {
    using ::Tpetra::Details::Behavior;
    using ::Tpetra::Details::dualViewStatusToString;
    using ::Tpetra::Details::ProfilingRegion;

    typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_type;

    ProfilingRegion profile_region("Tpetra::BlockCrsMatrix::copyAndPermuteNew");

    const bool debug = Behavior::debug();
    const bool verbose = Behavior::verbose();

    // Define this function prefix
    std::string prefix;
    {
      std::ostringstream os;
      const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
      os << "Proc " << myRank << ": BlockCrsMatrix::copyAndPermuteNew : " << std::endl;
      prefix = os.str();
    }

    // check if this already includes a local error
    if (* (this->localError_)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix
          << "The target object of the Import or Export is already in an error state."
          << std::endl;
      return;
    }

    //
    // Verbose input dual view status
    //
    if (verbose) {
      std::ostringstream os;
      os << prefix << std::endl
         << prefix << "  " << dualViewStatusToString (permuteToLIDs, "permuteToLIDs") << std::endl
         << prefix << "  " << dualViewStatusToString (permuteFromLIDs, "permuteFromLIDs") << std::endl;
      std::cerr << os.str ();
    }

    ///
    /// Check input valid
    ///
    if (permuteToLIDs.extent (0) != permuteFromLIDs.extent (0)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix
          << "permuteToLIDs.extent(0) = " << permuteToLIDs.extent (0)
          << " != permuteFromLIDs.extent(0) = " << permuteFromLIDs.extent(0)
          << "." << std::endl;
      return;
    }

    const this_type* src = dynamic_cast<const this_type* > (&source);
    if (src == NULL) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix
          << "The source (input) object of the Import or "
        "Export is either not a BlockCrsMatrix, or does not have the right "
        "template parameters.  checkSizes() should have caught this.  "
        "Please report this bug to the Tpetra developers." << std::endl;
      return;
    }

    bool lclErr = false;
#ifdef HAVE_TPETRA_DEBUG
    std::set<LO> invalidSrcCopyRows;
    std::set<LO> invalidDstCopyRows;
    std::set<LO> invalidDstCopyCols;
    std::set<LO> invalidDstPermuteCols;
    std::set<LO> invalidPermuteFromRows;
#endif // HAVE_TPETRA_DEBUG

    // Copy the initial sequence of rows that are the same.
    //
    // The two graphs might have different column Maps, so we need to
    // do this using global column indices.  This is purely local, so
    // we only need to check for local sameness of the two column
    // Maps.

#ifdef HAVE_TPETRA_DEBUG
    const map_type& srcRowMap = * (src->graph_.getRowMap ());
#endif // HAVE_TPETRA_DEBUG
    const map_type& dstRowMap = * (this->graph_.getRowMap ());
    const map_type& srcColMap = * (src->graph_.getColMap ());
    const map_type& dstColMap = * (this->graph_.getColMap ());
    const bool canUseLocalColumnIndices = srcColMap.locallySameAs (dstColMap);

    const size_t numPermute = static_cast<size_t> (permuteFromLIDs.extent(0));
    if (verbose) {
      std::ostringstream os;
      os << prefix
         << "canUseLocalColumnIndices: "
         << (canUseLocalColumnIndices ? "true" : "false")
         << ", numPermute: " << numPermute
         << std::endl;
      std::cerr << os.str ();
    }

    // work around for const object sync
    if (permuteToLIDs.need_sync_host()) {
      Kokkos::DualView<const local_ordinal_type*, device_type> permuteToLIDsTemp = permuteToLIDs;
      permuteToLIDsTemp.sync_host();
    }
    const auto permuteToLIDsHost = permuteToLIDs.view_host();
    if (permuteFromLIDs.need_sync_host()) {
      Kokkos::DualView<const local_ordinal_type*, device_type> permuteFromLIDsTemp = permuteFromLIDs;
      permuteFromLIDsTemp.sync_host();
    }
    const auto permuteFromLIDsHost = permuteFromLIDs.view_host();

    if (canUseLocalColumnIndices) {
      // Copy local rows that are the "same" in both source and target.
      for (LO localRow = 0; localRow < static_cast<LO> (numSameIDs); ++localRow) {
#ifdef HAVE_TPETRA_DEBUG
        if (! srcRowMap.isNodeLocalElement (localRow)) {
          lclErr = true;
          invalidSrcCopyRows.insert (localRow);
          continue; // skip invalid rows
        }
#endif // HAVE_TPETRA_DEBUG

        const LO* lclSrcCols;
        Scalar* vals;
        LO numEntries;
        // If this call fails, that means the mesh row local index is
        // invalid.  That means the Import or Export is invalid somehow.
        LO err = src->getLocalRowView (localRow, lclSrcCols, vals, numEntries);
        if (err != 0) {
          lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
          (void) invalidSrcCopyRows.insert (localRow);
#endif // HAVE_TPETRA_DEBUG
        }
        else {
          err = this->replaceLocalValues (localRow, lclSrcCols, vals, numEntries);
          if (err != numEntries) {
            lclErr = true;
            if (! dstRowMap.isNodeLocalElement (localRow)) {
#ifdef HAVE_TPETRA_DEBUG
              invalidDstCopyRows.insert (localRow);
#endif // HAVE_TPETRA_DEBUG
            }
            else {
              // Once there's an error, there's no sense in saving
              // time, so we check whether the column indices were
              // invalid.  However, only remember which ones were
              // invalid in a debug build, because that might take a
              // lot of space.
              for (LO k = 0; k < numEntries; ++k) {
                if (! dstColMap.isNodeLocalElement (lclSrcCols[k])) {
                  lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
                  (void) invalidDstCopyCols.insert (lclSrcCols[k]);
#endif // HAVE_TPETRA_DEBUG
                }
              }
            }
          }
        }
      } // for each "same" local row

      // Copy the "permute" local rows.
      for (size_t k = 0; k < numPermute; ++k) {
        const LO srcLclRow = static_cast<LO> (permuteFromLIDsHost(k));
        const LO dstLclRow = static_cast<LO> (permuteToLIDsHost(k));

        const LO* lclSrcCols;
        Scalar* vals;
        LO numEntries;
        LO err = src->getLocalRowView (srcLclRow, lclSrcCols, vals, numEntries);
        if (err != 0) {
          lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
          invalidPermuteFromRows.insert (srcLclRow);
#endif // HAVE_TPETRA_DEBUG
        }
        else {
          err = this->replaceLocalValues (dstLclRow, lclSrcCols, vals, numEntries);
          if (err != numEntries) {
            lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
            for (LO c = 0; c < numEntries; ++c) {
              if (! dstColMap.isNodeLocalElement (lclSrcCols[c])) {
                invalidDstPermuteCols.insert (lclSrcCols[c]);
              }
            }
#endif // HAVE_TPETRA_DEBUG
          }
        }
      }
    }
    else { // must convert column indices to global
      // Reserve space to store the destination matrix's local column indices.
      const size_t maxNumEnt = src->graph_.getNodeMaxNumRowEntries ();
      Teuchos::Array<LO> lclDstCols (maxNumEnt);

      // Copy local rows that are the "same" in both source and target.
      for (LO localRow = 0; localRow < static_cast<LO> (numSameIDs); ++localRow) {
        const LO* lclSrcCols;
        Scalar* vals;
        LO numEntries;
        // If this call fails, that means the mesh row local index is
        // invalid.  That means the Import or Export is invalid somehow.
        LO err = 0;
        try {
          err = src->getLocalRowView (localRow, lclSrcCols, vals, numEntries);
        } catch (std::exception& e) {
          if (debug) {
            std::ostringstream os;
            const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
            os << "Proc " << myRank << ": copyAndPermute: At \"same\" localRow "
               << localRow << ", src->getLocalRowView() threw an exception: "
               << e.what ();
            std::cerr << os.str ();
          }
          throw e;
        }

        if (err != 0) {
          lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
          invalidSrcCopyRows.insert (localRow);
#endif // HAVE_TPETRA_DEBUG
        }
        else {
          if (static_cast<size_t> (numEntries) > static_cast<size_t> (lclDstCols.size ())) {
            lclErr = true;
            if (debug) {
              std::ostringstream os;
              const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
              os << "Proc " << myRank << ": copyAndPermute: At \"same\" localRow "
                 << localRow << ", numEntries = " << numEntries << " > maxNumEnt = "
                 << maxNumEnt << std::endl;
              std::cerr << os.str ();
            }
          }
          else {
            // Convert the source matrix's local column indices to the
            // destination matrix's local column indices.
            Teuchos::ArrayView<LO> lclDstColsView = lclDstCols.view (0, numEntries);
            for (LO j = 0; j < numEntries; ++j) {
              lclDstColsView[j] = dstColMap.getLocalElement (srcColMap.getGlobalElement (lclSrcCols[j]));
              if (lclDstColsView[j] == Teuchos::OrdinalTraits<LO>::invalid ()) {
                lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
                invalidDstCopyCols.insert (lclDstColsView[j]);
#endif // HAVE_TPETRA_DEBUG
              }
            }
            try {
              err = this->replaceLocalValues (localRow, lclDstColsView.getRawPtr (), vals, numEntries);
            } catch (std::exception& e) {
              if (debug) {
                std::ostringstream os;
                const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
                os << "Proc " << myRank << ": copyAndPermute: At \"same\" localRow "
                   << localRow << ", this->replaceLocalValues() threw an exception: "
                   << e.what ();
                std::cerr << os.str ();
              }
              throw e;
            }
            if (err != numEntries) {
              lclErr = true;
              if (debug) {
                std::ostringstream os;
                const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
                os << "Proc " << myRank << ": copyAndPermute: At \"same\" "
                  "localRow " << localRow << ", this->replaceLocalValues "
                  "returned " << err << " instead of numEntries = "
                   << numEntries << std::endl;
                std::cerr << os.str ();
              }
            }
          }
        }
      }

      // Copy the "permute" local rows.
      for (size_t k = 0; k < numPermute; ++k) {
        const LO srcLclRow = static_cast<LO> (permuteFromLIDsHost(k));
        const LO dstLclRow = static_cast<LO> (permuteToLIDsHost(k));

        const LO* lclSrcCols;
        Scalar* vals;
        LO numEntries;
        LO err = 0;
        try {
          err = src->getLocalRowView (srcLclRow, lclSrcCols, vals, numEntries);
        } catch (std::exception& e) {
          if (debug) {
            std::ostringstream os;
            const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
            os << "Proc " << myRank << ": copyAndPermute: At \"permute\" "
              "srcLclRow " << srcLclRow << " and dstLclRow " << dstLclRow
               << ", src->getLocalRowView() threw an exception: " << e.what ();
            std::cerr << os.str ();
          }
          throw e;
        }

        if (err != 0) {
          lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
          invalidPermuteFromRows.insert (srcLclRow);
#endif // HAVE_TPETRA_DEBUG
        }
        else {
          if (static_cast<size_t> (numEntries) > static_cast<size_t> (lclDstCols.size ())) {
            lclErr = true;
          }
          else {
            // Convert the source matrix's local column indices to the
            // destination matrix's local column indices.
            Teuchos::ArrayView<LO> lclDstColsView = lclDstCols.view (0, numEntries);
            for (LO j = 0; j < numEntries; ++j) {
              lclDstColsView[j] = dstColMap.getLocalElement (srcColMap.getGlobalElement (lclSrcCols[j]));
              if (lclDstColsView[j] == Teuchos::OrdinalTraits<LO>::invalid ()) {
                lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
                invalidDstPermuteCols.insert (lclDstColsView[j]);
#endif // HAVE_TPETRA_DEBUG
              }
            }
            err = this->replaceLocalValues (dstLclRow, lclDstColsView.getRawPtr (), vals, numEntries);
            if (err != numEntries) {
              lclErr = true;
            }
          }
        }
      }
    }

    if (lclErr) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
#ifdef HAVE_TPETRA_DEBUG
      err << "copyAndPermute: The graph structure of the source of the "
        "Import or Export must be a subset of the graph structure of the "
        "target.  ";
      err << "invalidSrcCopyRows = [";
      for (typename std::set<LO>::const_iterator it = invalidSrcCopyRows.begin ();
           it != invalidSrcCopyRows.end (); ++it) {
        err << *it;
        typename std::set<LO>::const_iterator itp1 = it;
        itp1++;
        if (itp1 != invalidSrcCopyRows.end ()) {
          err << ",";
        }
      }
      err << "], invalidDstCopyRows = [";
      for (typename std::set<LO>::const_iterator it = invalidDstCopyRows.begin ();
           it != invalidDstCopyRows.end (); ++it) {
        err << *it;
        typename std::set<LO>::const_iterator itp1 = it;
        itp1++;
        if (itp1 != invalidDstCopyRows.end ()) {
          err << ",";
        }
      }
      err << "], invalidDstCopyCols = [";
      for (typename std::set<LO>::const_iterator it = invalidDstCopyCols.begin ();
           it != invalidDstCopyCols.end (); ++it) {
        err << *it;
        typename std::set<LO>::const_iterator itp1 = it;
        itp1++;
        if (itp1 != invalidDstCopyCols.end ()) {
          err << ",";
        }
      }
      err << "], invalidDstPermuteCols = [";
      for (typename std::set<LO>::const_iterator it = invalidDstPermuteCols.begin ();
           it != invalidDstPermuteCols.end (); ++it) {
        err << *it;
        typename std::set<LO>::const_iterator itp1 = it;
        itp1++;
        if (itp1 != invalidDstPermuteCols.end ()) {
          err << ",";
        }
      }
      err << "], invalidPermuteFromRows = [";
      for (typename std::set<LO>::const_iterator it = invalidPermuteFromRows.begin ();
           it != invalidPermuteFromRows.end (); ++it) {
        err << *it;
        typename std::set<LO>::const_iterator itp1 = it;
        itp1++;
        if (itp1 != invalidPermuteFromRows.end ()) {
          err << ",";
        }
      }
      err << "]" << std::endl;
#else
      err << "copyAndPermute: The graph structure of the source of the "
        "Import or Export must be a subset of the graph structure of the "
        "target." << std::endl;
#endif // HAVE_TPETRA_DEBUG
    }

    if (debug) {
      std::ostringstream os;
      const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
      const bool lclSuccess = ! (* (this->localError_));
      os << "*** Proc " << myRank << ": copyAndPermute "
         << (lclSuccess ? "succeeded" : "FAILED");
      if (lclSuccess) {
        os << std::endl;
      } else {
        os << ": error messages: " << this->errorMessages (); // comes w/ endl
      }
      std::cerr << os.str ();
    }
  }

  namespace { // (anonymous)

    /// \brief Return the (maximum) number of bytes required to pack a
    ///   block row's entries.
    ///
    /// \param numEnt [in] Number of block entries in the row.
    ///
    /// \param numBytesPerValue [in] Maximum number of bytes per
    ///   scalar (not block!) entry (value) of the row.
    ///
    /// \param blkSize [in] Block size of the block sparse matrix.
    ///
    /// If \c Scalar (the type of entries in the matrix) is a plain
    /// old data (POD) type like \c float or \c double, or a struct of
    /// POD (like <tt>std::complex<double></tt>), then the second
    /// argument is just <tt>sizeof(Scalar)</tt>.  If a \c Scalar
    /// instance has a size determined at run time (e.g., when calling
    /// its constructor), then the second argument is the result of
    /// <tt>PackTraits<Scalar>::packValueCount</tt>, called on a
    /// <tt>Scalar</tt> value with the correct run-time size.
    template<class LO, class GO, class D>
    size_t
    packRowCount (const size_t numEnt,
                  const size_t numBytesPerValue,
                  const size_t blkSize)
    {
      using ::Tpetra::Details::PackTraits;

      if (numEnt == 0) {
        // Empty rows always take zero bytes, to ensure sparsity.
        return 0;
      }
      else {
        // We store the number of entries as a local index (LO).
        LO numEntLO = 0; // packValueCount wants this.
        GO gid;
        const size_t numEntLen = PackTraits<LO, D>::packValueCount (numEntLO);
        const size_t gidsLen = numEnt * PackTraits<GO, D>::packValueCount (gid);
        const size_t valsLen = numEnt * numBytesPerValue * blkSize * blkSize;
        return numEntLen + gidsLen + valsLen;
      }
    }

    /// \brief Unpack and return the number of (block) entries in the
    ///   packed row.
    ///
    /// \param imports [in] All the packed data.
    /// \param offset [in] Index of \c imports at which the row starts.
    /// \param numBytes [in] Number of bytes in the packed row.
    /// \param numBytesPerValue [in] Maximum number of bytes per
    ///   scalar (not block!) entry (value) of the row.
    ///
    /// \return Number of (block) entries in the packed row.
    template<class ST, class LO, class GO, class D>
    size_t
    unpackRowCount (const typename ::Tpetra::Details::PackTraits<LO, D>::input_buffer_type& imports,
                    const size_t offset,
                    const size_t numBytes,
                    const size_t numBytesPerValue)
    {
      using Kokkos::subview;
      using ::Tpetra::Details::PackTraits;

      if (numBytes == 0) {
        // Empty rows always take zero bytes, to ensure sparsity.
        return static_cast<size_t> (0);
      }
      else {
        LO numEntLO = 0;
#ifdef HAVE_TPETRA_DEBUG
        const size_t theNumBytes = PackTraits<LO, D>::packValueCount (numEntLO);
        TEUCHOS_TEST_FOR_EXCEPTION(
          theNumBytes > numBytes, std::logic_error, "unpackRowCount: "
          "theNumBytes = " << theNumBytes << " < numBytes = " << numBytes
          << ".");
#endif // HAVE_TPETRA_DEBUG
        const char* const inBuf = imports.data () + offset;
        const size_t actualNumBytes = PackTraits<LO, D>::unpackValue (numEntLO, inBuf);
#ifdef HAVE_TPETRA_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION(
          actualNumBytes > numBytes, std::logic_error, "unpackRowCount: "
          "actualNumBytes = " << actualNumBytes << " < numBytes = " << numBytes
          << ".");
#else
        (void) actualNumBytes;
#endif // HAVE_TPETRA_DEBUG
        return static_cast<size_t> (numEntLO);
      }
    }

    /// \brief Pack the block row (stored in the input arrays).
    ///
    /// \return The number of bytes packed.
    ///
    /// \note This function is not called packRow, because Intel 16
    /// has a bug that makes it confuse this packRow with
    /// Tpetra::RowMatrix::packRow.
    template<class ST, class LO, class GO, class D>
    size_t
    packRowForBlockCrs (const typename ::Tpetra::Details::PackTraits<LO, D>::output_buffer_type exports,
                        const size_t offset,
                        const size_t numEnt,
                        const typename ::Tpetra::Details::PackTraits<GO, D>::input_array_type& gidsIn,
                        const typename ::Tpetra::Details::PackTraits<ST, D>::input_array_type& valsIn,
                        const size_t numBytesPerValue,
                        const size_t blockSize)
    {
      using Kokkos::subview;
      using ::Tpetra::Details::PackTraits;

      if (numEnt == 0) {
        // Empty rows always take zero bytes, to ensure sparsity.
        return 0;
      }
      const size_t numScalarEnt = numEnt * blockSize * blockSize;

      const GO gid = 0; // packValueCount wants this
      const LO numEntLO = static_cast<size_t> (numEnt);

      const size_t numEntBeg = offset;
      const size_t numEntLen = PackTraits<LO, D>::packValueCount (numEntLO);
      const size_t gidsBeg = numEntBeg + numEntLen;
      const size_t gidsLen = numEnt * PackTraits<GO, D>::packValueCount (gid);
      const size_t valsBeg = gidsBeg + gidsLen;
      const size_t valsLen = numScalarEnt * numBytesPerValue;

      char* const numEntOut = exports.data () + numEntBeg;
      char* const gidsOut = exports.data () + gidsBeg;
      char* const valsOut = exports.data () + valsBeg;

      size_t numBytesOut = 0;
      int errorCode = 0;
      numBytesOut += PackTraits<LO, D>::packValue (numEntOut, numEntLO);

      {
        Kokkos::pair<int, size_t> p;
        p = PackTraits<GO, D>::packArray (gidsOut, gidsIn.data (), numEnt);
        errorCode += p.first;
        numBytesOut += p.second;

        p = PackTraits<ST, D>::packArray (valsOut, valsIn.data (), numScalarEnt);
        errorCode += p.first;
        numBytesOut += p.second;
      }

      const size_t expectedNumBytes = numEntLen + gidsLen + valsLen;
      TEUCHOS_TEST_FOR_EXCEPTION(
        numBytesOut != expectedNumBytes, std::logic_error, "packRow: "
        "numBytesOut = " << numBytesOut << " != expectedNumBytes = "
        << expectedNumBytes << ".");

      TEUCHOS_TEST_FOR_EXCEPTION(
        errorCode != 0, std::runtime_error, "packRow: "
        "PackTraits::packArray returned a nonzero error code");

      return numBytesOut;
    }

    // Return the number of bytes actually read / used.
    template<class ST, class LO, class GO, class D>
    size_t
    unpackRowForBlockCrs (const typename ::Tpetra::Details::PackTraits<GO, D>::output_array_type& gidsOut,
                          const typename ::Tpetra::Details::PackTraits<ST, D>::output_array_type& valsOut,
                          const typename ::Tpetra::Details::PackTraits<int, D>::input_buffer_type& imports,
                          const size_t offset,
                          const size_t numBytes,
                          const size_t numEnt,
                          const size_t numBytesPerValue,
                          const size_t blockSize)
    {
      using ::Tpetra::Details::PackTraits;

      if (numBytes == 0) {
        // Rows with zero bytes always have zero entries.
        return 0;
      }
      const size_t numScalarEnt = numEnt * blockSize * blockSize;
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (imports.extent (0)) <= offset,
        std::logic_error, "unpackRow: imports.extent(0) = "
        << imports.extent (0) << " <= offset = " << offset << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (imports.extent (0)) < offset + numBytes,
        std::logic_error, "unpackRow: imports.extent(0) = "
        << imports.extent (0) << " < offset + numBytes = "
        << (offset + numBytes) << ".");

      const GO gid = 0; // packValueCount wants this
      const LO lid = 0; // packValueCount wants this

      const size_t numEntBeg = offset;
      const size_t numEntLen = PackTraits<LO, D>::packValueCount (lid);
      const size_t gidsBeg = numEntBeg + numEntLen;
      const size_t gidsLen = numEnt * PackTraits<GO, D>::packValueCount (gid);
      const size_t valsBeg = gidsBeg + gidsLen;
      const size_t valsLen = numScalarEnt * numBytesPerValue;

      const char* const numEntIn = imports.data () + numEntBeg;
      const char* const gidsIn = imports.data () + gidsBeg;
      const char* const valsIn = imports.data () + valsBeg;

      size_t numBytesOut = 0;
      int errorCode = 0;
      LO numEntOut;
      numBytesOut += PackTraits<LO, D>::unpackValue (numEntOut, numEntIn);
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (numEntOut) != numEnt, std::logic_error,
        "unpackRow: Expected number of entries " << numEnt
        << " != actual number of entries " << numEntOut << ".");

      {
        Kokkos::pair<int, size_t> p;
        p = PackTraits<GO, D>::unpackArray (gidsOut.data (), gidsIn, numEnt);
        errorCode += p.first;
        numBytesOut += p.second;

        p = PackTraits<ST, D>::unpackArray (valsOut.data (), valsIn, numScalarEnt);
        errorCode += p.first;
        numBytesOut += p.second;
      }

      TEUCHOS_TEST_FOR_EXCEPTION(
        numBytesOut != numBytes, std::logic_error, "unpackRow: numBytesOut = "
        << numBytesOut << " != numBytes = " << numBytes << ".");

      const size_t expectedNumBytes = numEntLen + gidsLen + valsLen;
      TEUCHOS_TEST_FOR_EXCEPTION(
        numBytesOut != expectedNumBytes, std::logic_error, "unpackRow: "
        "numBytesOut = " << numBytesOut << " != expectedNumBytes = "
        << expectedNumBytes << ".");

      TEUCHOS_TEST_FOR_EXCEPTION(
        errorCode != 0, std::runtime_error, "unpackRow: "
        "PackTraits::unpackArray returned a nonzero error code");

      return numBytesOut;
    }
  } // namespace (anonymous)

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  packAndPrepareNew (const ::Tpetra::SrcDistObject& source,
                     const Kokkos::DualView<const local_ordinal_type*, device_type>& exportLIDs,
                     Kokkos::DualView<impl_scalar_type*, buffer_device_type>& exports,
                     const Kokkos::DualView<size_t*, buffer_device_type>& numPacketsPerLID,
                     size_t& constantNumPackets,
                     Distributor& /* distor */)
  {
    using ::Tpetra::Details::Behavior;
    using ::Tpetra::Details::dualViewStatusToString;
    using ::Tpetra::Details::ProfilingRegion;
    using ::Tpetra::Details::PackTraits;

    typedef typename Kokkos::View<int*, device_type>::HostMirror::execution_space host_exec;

    typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_type;

    ProfilingRegion profile_region("Tpetra::BlockCrsMatrix::packAndPrepareNew");

    const bool debug = Behavior::debug();
    const bool verbose = Behavior::verbose();

    // Define this function prefix
    std::string prefix;
    {
      std::ostringstream os;
      const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
      os << "Proc " << myRank << ": BlockCrsMatrix::packAndPrepareNew : " << std::endl;
      prefix = os.str();
    }

    // check if this already includes a local error
    if (* (this->localError_)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix
          << "The target object of the Import or Export is already in an error state."
          << std::endl;
      return;
    }

    //
    // Verbose input dual view status
    //
    if (verbose) {
      std::ostringstream os;
      os << prefix << std::endl
         << prefix << "  " << dualViewStatusToString (exportLIDs, "exportLIDs") << std::endl
         << prefix << "  " << dualViewStatusToString (exports, "exports") << std::endl
         << prefix << "  " << dualViewStatusToString (numPacketsPerLID, "numPacketsPerLID") << std::endl;
      std::cerr << os.str ();
    }

    ///
    /// Check input valid
    ///
    if (exportLIDs.extent (0) != numPacketsPerLID.extent (0)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix
          << "exportLIDs.extent(0) = " << exportLIDs.extent (0)
          << " != numPacketsPerLID.extent(0) = " << numPacketsPerLID.extent(0)
          << "." << std::endl;
      return;
    }

    const this_type* src = dynamic_cast<const this_type* > (&source);
    if (src == NULL) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix
          << "The source (input) object of the Import or "
        "Export is either not a BlockCrsMatrix, or does not have the right "
        "template parameters.  checkSizes() should have caught this.  "
        "Please report this bug to the Tpetra developers." << std::endl;
      return;
    }

    // Graphs and matrices are allowed to have a variable number of
    // entries per row.  We could test whether all rows have the same
    // number of entries, but DistObject can only use this
    // optimization if all rows on _all_ processes have the same
    // number of entries.  Rather than do the all-reduce necessary to
    // test for this unlikely case, we tell DistObject (by setting
    // constantNumPackets to zero) to assume that different rows may
    // have different numbers of entries.
    constantNumPackets = 0;

    // const values
    const crs_graph_type& srcGraph = src->graph_;
    const size_t blockSize = static_cast<size_t> (src->getBlockSize ());
    const size_t numExportLIDs = exportLIDs.extent (0);
    const size_t numBytesPerValue =
      PackTraits<impl_scalar_type, host_exec>
      ::packValueCount(this->val_.extent(0) ? this->val_.view_host()(0) : impl_scalar_type());

    // Compute the number of bytes ("packets") per row to pack.  While
    // we're at it, compute the total # of block entries to send, and
    // the max # of block entries in any of the rows we're sending.

    Impl::BlockCrsRowStruct<size_t> rowReducerStruct;

    // Graph information is on host; let's do this on host parallel reduce
    // Sync necessary data to host

    if (exportLIDs.need_sync_host()) {
      Kokkos::DualView<const local_ordinal_type*, device_type> exportLIDsTemp = exportLIDs;
      exportLIDsTemp.sync_host();
    }
    auto exportLIDsHost = exportLIDs.view_host();
    auto numPacketsPerLIDHost = numPacketsPerLID.view_host(); // we will modify this
    {
      typedef Impl::BlockCrsReducer<Impl::BlockCrsRowStruct<size_t>,host_exec> reducer_type;
      const auto policy = Kokkos::RangePolicy<host_exec>(size_t(0), numExportLIDs);
      Kokkos::parallel_reduce
        (policy,
         [=](const size_t &i, typename reducer_type::value_type &update) {
          const LO lclRow = exportLIDsHost(i);
          size_t numEnt = srcGraph.getNumEntriesInLocalRow (lclRow);
          numEnt = (numEnt == Teuchos::OrdinalTraits<size_t>::invalid () ? 0 : numEnt);

          const size_t numBytes = packRowCount<LO, GO, host_exec> (numEnt, numBytesPerValue, blockSize);
          numPacketsPerLIDHost(i) = numBytes;
          update += typename reducer_type::value_type(numEnt, numBytes, numEnt);
        }, rowReducerStruct);
      {
        Kokkos::DualView<size_t*, buffer_device_type> numPacketsPerLIDTemp = numPacketsPerLID;
        numPacketsPerLIDTemp.modify_host();
      }
    }

    // Compute the number of bytes ("packets") per row to pack.  While
    // we're at it, compute the total # of block entries to send, and
    // the max # of block entries in any of the rows we're sending.
    const size_t totalNumBytes   = rowReducerStruct.totalNumBytes;
    const size_t totalNumEntries = rowReducerStruct.totalNumEntries;
    const size_t maxRowLength    = rowReducerStruct.maxRowLength;

    if (verbose) {
      std::ostringstream os;
      os << prefix
         << "totalNumBytes = " << totalNumBytes << ", totalNumEntries = " << totalNumEntries
         << std::endl;
      std::cerr << os.str ();
    }

    // We use a "struct of arrays" approach to packing each row's
    // entries.  All the column indices (as global indices) go first,
    // then all their owning process ranks, and then the values.
    exports.resize (totalNumBytes/numBytesPerValue);
    if (totalNumEntries > 0) {
      // exports is resized in the above and we assume the data on device is invalidated.
      Kokkos::View<char*,host_exec> exportsByteHost ((char*)exports.view_host().data(), totalNumBytes);

      // Current position (in bytes) in the 'exports' output array.
      Kokkos::View<size_t*, host_exec> offset("offset", numExportLIDs+1);
      {
        const auto policy = Kokkos::RangePolicy<host_exec>(size_t(0), numExportLIDs+1);
        Kokkos::parallel_scan
          (policy,
           [=](const size_t &i, size_t &update, const bool &final) {
            if (final) offset(i) = update;
            update += (i == numExportLIDs ? 0 : numPacketsPerLIDHost(i));
          });
      }
      if (offset(numExportLIDs) != totalNumBytes) {
        std::ostream& err = this->markLocalErrorAndGetStream ();
        err << prefix
            << "At end of method, the final offset (in bytes) "
            << offset(numExportLIDs) << " does not equal the total number of bytes packed "
            << totalNumBytes << ".  "
            << "Please report this bug to the Tpetra developers." << std::endl;
        return;
      }

      // For each block row of the matrix owned by the calling
      // process, pack that block row's column indices and values into
      // the exports array.

      // Source matrix's column Map.  We verified in checkSizes() that
      // the column Map exists (is not null).
      const map_type& srcColMap = * (srcGraph.getColMap ());

      // Pack the data for each row to send, into the 'exports' buffer.
      {
        typedef Kokkos::TeamPolicy<host_exec> policy_type;
        const auto policy =
          policy_type(numExportLIDs, 1, 1)
          .set_scratch_size(0, Kokkos::PerTeam(sizeof(GO)*maxRowLength));
        Kokkos::parallel_for
          (policy,
           [=](const typename policy_type::member_type &member) {
            const size_t i = member.league_rank();
            Kokkos::View<GO*, typename host_exec::scratch_memory_space>
              gblColInds(member.team_scratch(0), maxRowLength);

            const LO  lclRowInd = exportLIDsHost(i);
            const LO* lclColIndsRaw;
            Scalar* valsRaw;
            LO numEntLO;
            // It's OK to ignore the return value, since if the calling
            // process doesn't own that local row, then the number of
            // entries in that row on the calling process is zero.
            (void) src->getLocalRowView (lclRowInd, lclColIndsRaw, valsRaw, numEntLO);

            const size_t numEnt = static_cast<size_t> (numEntLO);
            Kokkos::View<const LO*,host_exec> lclColInds (lclColIndsRaw, numEnt);

            // Convert column indices from local to global.
            for (size_t j = 0; j < numEnt; ++j)
              gblColInds(j) = srcColMap.getGlobalElement (lclColInds(j));

            // Kyungjoo: additional wrapping scratch view is necessary
            //   the following function interface need the same execution space
            //   host scratch space somehow is not considered same as the host_exec
            // Copy the row's data into the current spot in the exports array.
            const size_t numBytes = packRowForBlockCrs<impl_scalar_type,LO,GO,host_exec>
              (exportsByteHost,
               offset(i),
               numEnt,
               Kokkos::View<const GO*,host_exec>(gblColInds.data(), maxRowLength),
               Kokkos::View<const impl_scalar_type*,host_exec>(reinterpret_cast<const impl_scalar_type*>(valsRaw), numEnt*blockSize*blockSize),
               numBytesPerValue,
               blockSize);

            // numBytes should be same as the difference between offsets
            if (debug) {
              const size_t offsetDiff = offset(i+1) - offset(i);
              if (numBytes != offsetDiff) {
                std::ostringstream os;
                os << prefix
                   << "numBytes computed from packRowForBlockCrs is different from "
                   << "precomputed offset values, LID = " << i << std::endl;
                std::cerr << os.str ();
              }
            }
          }); // for each LID (of a row) to send
      }
    } // if totalNumEntries > 0

    if (debug) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      const bool lclSuccess = ! (* (this->localError_));
      err << prefix
          << (lclSuccess ? "succeeded" : "FAILED")
          << " (totalNumEntries = " << totalNumEntries << ") ***" << std::endl;
    }
  }


  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  unpackAndCombineNew (const Kokkos::DualView<const local_ordinal_type*, device_type>& importLIDs,
                       const Kokkos::DualView<const impl_scalar_type*, buffer_device_type>& imports,
                       const Kokkos::DualView<const size_t*, buffer_device_type>& numPacketsPerLID,
                       const size_t /* constantNumPackets */,
                       Distributor& /* distor */,
                       const CombineMode combineMode)
  {
    using ::Tpetra::Details::Behavior;
    using ::Tpetra::Details::dualViewStatusToString;
    using ::Tpetra::Details::ProfilingRegion;
    using ::Tpetra::Details::PackTraits;

    typedef typename Kokkos::View<int*, device_type>::HostMirror::execution_space host_exec;

    ProfilingRegion profile_region("Tpetra::BlockCrsMatrix::unpackAndCombineNew");

    const bool debug = Behavior::debug();
    const bool verbose = Behavior::verbose();

    // Define this function prefix
    std::string prefix;
    {
      std::ostringstream os;
      const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
      os << "Proc " << myRank << ": BlockCrsMatrix::unpackAndCombineNew : " << std::endl;
      prefix = os.str();
    }

    // check if this already includes a local error
    if (* (this->localError_)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix
          << "The target object of the Import or Export is already in an error state."
          << std::endl;
      return;
    }

    //
    // Verbose input dual view status
    //
    if (verbose) {
      std::ostringstream os;
      os << prefix << std::endl
         << prefix << "  " << dualViewStatusToString (importLIDs, "importLIDs") << std::endl
         << prefix << "  " << dualViewStatusToString (imports, "imports") << std::endl
         << prefix << "  " << dualViewStatusToString (numPacketsPerLID, "numPacketsPerLID") << std::endl;
      std::cerr << os.str ();
    }

    ///
    /// Check input valid
    ///
    if (importLIDs.extent (0) != numPacketsPerLID.extent (0)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix
          << "importLIDs.extent(0) = " << importLIDs.extent (0)
          << " != numPacketsPerLID.extent(0) = " << numPacketsPerLID.extent(0)
          << "." << std::endl;
      return;
    }

    if (combineMode != ADD     && combineMode != INSERT &&
        combineMode != REPLACE && combineMode != ABSMAX &&
        combineMode != ZERO) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix
          << "Invalid CombineMode value " << combineMode << ".  Valid "
          << "values include ADD, INSERT, REPLACE, ABSMAX, and ZERO."
          << std::endl;
      return;
    }

    // Target matrix's column Map.  Use to convert the global column
    // indices in the receive buffer to local indices.  We verified in
    // checkSizes() that the column Map exists (is not null).
    const map_type& tgtColMap = * (this->graph_.getColMap ());

    // Const values
    const size_t blockSize = this->getBlockSize ();
    const size_t numImportLIDs = importLIDs.extent(0);
    const size_t numBytesPerValue
      = PackTraits<impl_scalar_type, host_exec>
      ::packValueCount(this->val_.extent(0) ? this->val_.view_host()(0) : impl_scalar_type());
    const size_t maxRowNumEnt = graph_.getNodeMaxNumRowEntries ();
    const size_t maxRowNumScalarEnt = maxRowNumEnt * blockSize * blockSize;

    // Early returns
    if (combineMode == ZERO || numImportLIDs == 0) {
      if (debug) {
        std::ostringstream os;
        os << prefix << "Nothing to do" << std::endl;
        std::cerr << os.str ();
      }
      return; // nothing to do; no need to combine entries
    }

    if (verbose) {
      std::ostringstream os;
      os << prefix << "Getting ready" << std::endl;
      std::cerr << os.str ();
    }

    Kokkos::View<char*,host_exec> importsByteHost((char*)imports.view_host().data(), imports.extent(0)*numBytesPerValue);

    if (importLIDs.need_sync_host()) {
      Kokkos::DualView<const local_ordinal_type*, device_type> importLIDsTemp = importLIDs;
      importLIDsTemp.sync_host();
    }
    const auto importLIDsHost = importLIDs.view_host();

    if (numPacketsPerLID.need_sync_host()) {
      Kokkos::DualView<const size_t*, buffer_device_type> numPacketsPerLIDTemp = numPacketsPerLID;
      numPacketsPerLIDTemp.sync_host();
    }
    const auto numPacketsPerLIDHost = numPacketsPerLID.view_host();

    Kokkos::View<size_t*,host_exec> offset("offset", numImportLIDs+1);
    {
      const auto policy = Kokkos::RangePolicy<host_exec>(size_t(0), numImportLIDs+1);
      Kokkos::parallel_scan
        (policy,
         [=](const size_t &i, size_t &update, const bool &final) {
          if (final) offset(i) = update;
          update += (i == numImportLIDs ? 0 : numPacketsPerLIDHost(i));
        });
    }

    // this variable does not matter with a race condition (just error flag)
    Kokkos::View<bool,host_exec,Kokkos::MemoryTraits<Kokkos::Atomic> > errorDuringUnpack("errorDuringUnpack");
    errorDuringUnpack() = false;
    {
      typedef Kokkos::TeamPolicy<host_exec> policy_type;
      const auto policy =
        policy_type(numImportLIDs, 1, 1)
        .set_scratch_size(0, Kokkos::PerTeam(sizeof(GO(0))*maxRowNumEnt+
                                             sizeof(LO(0))*maxRowNumEnt+
                                             sizeof(impl_scalar_type(0))*maxRowNumScalarEnt));

      typedef typename host_exec::scratch_memory_space host_scratch_space;
      Kokkos::parallel_for
        (policy,
         [=](const typename policy_type::member_type &member) {
          const size_t i = member.league_rank();
          Kokkos::View<GO*,host_scratch_space> gblColInds(member.team_scratch(0), maxRowNumEnt);
          Kokkos::View<LO*,host_scratch_space> lclColInds(member.team_scratch(0), maxRowNumEnt);
          Kokkos::View<impl_scalar_type*,host_scratch_space> vals(member.team_scratch(0), maxRowNumScalarEnt);

          const size_t offval = offset(i);
          const LO lclRow = importLIDsHost(i);
          const size_t numBytes = numPacketsPerLIDHost(i);
          const size_t numEnt =
            unpackRowCount<impl_scalar_type, LO, GO, host_exec>
            (importsByteHost, offval, numBytes, numBytesPerValue);

          if (numBytes > 0) {
            if (numEnt > maxRowNumEnt) {
              errorDuringUnpack() = true;
              if (debug) {
                std::ostringstream os;
                os << prefix
                   << "At i = " << i << ", numEnt = " << numEnt
                   << " > maxRowNumEnt = " << maxRowNumEnt
                   << std::endl;
                std::cerr << os.str();
              }
            }
          }
          const size_t numScalarEnt = numEnt*blockSize*blockSize;
          auto gidsOut = Kokkos::subview(gblColInds, Kokkos::pair<size_t,size_t>(0, numEnt));
          auto lidsOut = Kokkos::subview(lclColInds, Kokkos::pair<size_t,size_t>(0, numEnt));
          auto valsOut = Kokkos::subview(vals,       Kokkos::pair<size_t,size_t>(0, numScalarEnt));

          // Kyungjoo: additional wrapping scratch view is necessary
          //   the following function interface need the same execution space
          //   host scratch space somehow is not considered same as the host_exec
          const size_t numBytesOut =
            unpackRowForBlockCrs<impl_scalar_type, LO, GO, host_exec>
            (Kokkos::View<GO*,host_exec>(gidsOut.data(), numEnt),
             Kokkos::View<impl_scalar_type*,host_exec>(valsOut.data(), numScalarEnt),
             importsByteHost,
             offval, numBytes, numEnt,
             numBytesPerValue, blockSize);

          if (numBytes != numBytesOut) {
            errorDuringUnpack() = true;
            if (debug) {
              std::ostringstream os;
              os << prefix
                 << "At i = " << i << ", numBytes = " << numBytes
                 << " != numBytesOut = " << numBytesOut << "."
                 << std::endl;
              std::cerr << os.str();
            }
          }

          // Convert incoming global indices to local indices.
          for (size_t k = 0; k < numEnt; ++k) {
            lidsOut(k) = tgtColMap.getLocalElement (gidsOut(k));
            if (lidsOut(k) == Teuchos::OrdinalTraits<LO>::invalid ()) {
              errorDuringUnpack() = true;
              if (debug) {
                std::ostringstream os;
                os << prefix
                    << "At i = " << i << ", GID " << gidsOut(k)
                    << " is not owned by the calling process."
                    << std::endl;
                std::cerr << os.str();
              }
            }
          }

          // Combine the incoming data with the matrix's current data.
          LO numCombd = 0;
          const LO* const lidsRaw = const_cast<const LO*> (lidsOut.data ());
          const Scalar* const valsRaw =
            reinterpret_cast<const Scalar*> (const_cast<const impl_scalar_type*> (valsOut.data ()));
          if (combineMode == ADD) {
            numCombd = this->sumIntoLocalValues (lclRow, lidsRaw, valsRaw, numEnt);
          } else if (combineMode == INSERT || combineMode == REPLACE) {
            numCombd = this->replaceLocalValues (lclRow, lidsRaw, valsRaw, numEnt);
          } else if (combineMode == ABSMAX) {
            numCombd = this->absMaxLocalValues (lclRow, lidsRaw, valsRaw, numEnt);
          }

          if (static_cast<LO> (numEnt) != numCombd) {
            errorDuringUnpack() = true;
            if (debug) {
            std::ostream& err = this->markLocalErrorAndGetStream ();
            err << prefix
                << "At i = " << i << ", numEnt = " << numEnt
                << " != numCombd = " << numCombd << "."
                << std::endl;
            }
          }
        }); // for each import LID i
    }

    if (errorDuringUnpack()) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << "Unpacking failed.";
      if (!debug)
        err << "  Please run again with a debug build to get more verbose diagnostic output.";
      err << std::endl;
    }

    if (debug) {
      std::ostringstream os;
      const bool lclSuccess = ! (* (this->localError_));
      os << prefix
         << (lclSuccess ? "succeeded" : "FAILED")
         << std::endl;
      std::cerr << os.str ();
    }
  }


  template<class Scalar, class LO, class GO, class Node>
  std::string
  BlockCrsMatrix<Scalar, LO, GO, Node>::description () const
  {
    using Teuchos::TypeNameTraits;
    std::ostringstream os;
    os << "\"Tpetra::BlockCrsMatrix\": { "
       << "Template parameters: { "
       << "Scalar: " << TypeNameTraits<Scalar>::name ()
       << "LO: " << TypeNameTraits<LO>::name ()
       << "GO: " << TypeNameTraits<GO>::name ()
       << "Node: " << TypeNameTraits<Node>::name ()
       << " }"
       << ", Label: \"" << this->getObjectLabel () << "\""
       << ", Global dimensions: ["
       << graph_.getDomainMap ()->getGlobalNumElements () << ", "
       << graph_.getRangeMap ()->getGlobalNumElements () << "]"
       << ", Block size: " << getBlockSize ()
       << " }";
    return os.str ();
  }


  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::CommRequest;
    using Teuchos::FancyOStream;
    using Teuchos::getFancyOStream;
    using Teuchos::ireceive;
    using Teuchos::isend;
    using Teuchos::outArg;
    using Teuchos::TypeNameTraits;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    // using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    using Teuchos::RCP;
    using Teuchos::wait;
    using std::endl;
#ifdef HAVE_TPETRA_DEBUG
    const char prefix[] = "Tpetra::Experimental::BlockCrsMatrix::describe: ";
#endif // HAVE_TPETRA_DEBUG

    // Set default verbosity if applicable.
    const Teuchos::EVerbosityLevel vl =
      (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

    if (vl == VERB_NONE) {
      return; // print nothing
    }

    // describe() always starts with a tab before it prints anything.
    Teuchos::OSTab tab0 (out);

    out << "\"Tpetra::BlockCrsMatrix\":" << endl;
    Teuchos::OSTab tab1 (out);

    out << "Template parameters:" << endl;
    {
      Teuchos::OSTab tab2 (out);
      out << "Scalar: " << TypeNameTraits<Scalar>::name () << endl
          << "LO: " << TypeNameTraits<LO>::name () << endl
          << "GO: " << TypeNameTraits<GO>::name () << endl
          << "Node: " << TypeNameTraits<Node>::name () << endl;
    }
    out << "Label: \"" << this->getObjectLabel () << "\"" << endl
        << "Global dimensions: ["
        << graph_.getDomainMap ()->getGlobalNumElements () << ", "
        << graph_.getRangeMap ()->getGlobalNumElements () << "]" << endl;

    const LO blockSize = getBlockSize ();
    out << "Block size: " << blockSize << endl;

    // constituent objects
    if (vl >= VERB_MEDIUM) {
      const Teuchos::Comm<int>& comm = * (graph_.getMap ()->getComm ());
      const int myRank = comm.getRank ();
      if (myRank == 0) {
        out << "Row Map:" << endl;
      }
      getRowMap()->describe (out, vl);

      if (! getColMap ().is_null ()) {
        if (getColMap() == getRowMap()) {
          if (myRank == 0) {
            out << "Column Map: same as row Map" << endl;
          }
        }
        else {
          if (myRank == 0) {
            out << "Column Map:" << endl;
          }
          getColMap ()->describe (out, vl);
        }
      }
      if (! getDomainMap ().is_null ()) {
        if (getDomainMap () == getRowMap ()) {
          if (myRank == 0) {
            out << "Domain Map: same as row Map" << endl;
          }
        }
        else if (getDomainMap () == getColMap ()) {
          if (myRank == 0) {
            out << "Domain Map: same as column Map" << endl;
          }
        }
        else {
          if (myRank == 0) {
            out << "Domain Map:" << endl;
          }
          getDomainMap ()->describe (out, vl);
        }
      }
      if (! getRangeMap ().is_null ()) {
        if (getRangeMap () == getDomainMap ()) {
          if (myRank == 0) {
            out << "Range Map: same as domain Map" << endl;
          }
        }
        else if (getRangeMap () == getRowMap ()) {
          if (myRank == 0) {
            out << "Range Map: same as row Map" << endl;
          }
        }
        else {
          if (myRank == 0) {
            out << "Range Map: " << endl;
          }
          getRangeMap ()->describe (out, vl);
        }
      }
    }

    if (vl >= VERB_EXTREME) {
      // FIXME (mfh 26 May 2016) It's not nice for this method to sync
      // to host, since it's supposed to be const.  However, that's
      // the easiest and least memory-intensive way to implement this
      // method.
      typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_type;
      const_cast<this_type*> (this)->sync_host ();

#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION
        (this->need_sync_host (), std::logic_error,
         prefix << "Right after sync to host, the matrix claims that it needs "
         "sync to host.  Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

      const Teuchos::Comm<int>& comm = * (graph_.getMap ()->getComm ());
      const int myRank = comm.getRank ();
      const int numProcs = comm.getSize ();

      // Print the calling process' data to the given output stream.
      RCP<std::ostringstream> lclOutStrPtr (new std::ostringstream ());
      RCP<FancyOStream> osPtr = getFancyOStream (lclOutStrPtr);
      FancyOStream& os = *osPtr;
      os << "Process " << myRank << ":" << endl;
      Teuchos::OSTab tab2 (os);

      const map_type& meshRowMap = * (graph_.getRowMap ());
      const map_type& meshColMap = * (graph_.getColMap ());
      for (LO meshLclRow = meshRowMap.getMinLocalIndex ();
           meshLclRow <= meshRowMap.getMaxLocalIndex ();
           ++meshLclRow) {
        const GO meshGblRow = meshRowMap.getGlobalElement (meshLclRow);
        os << "Row " << meshGblRow << ": {";

        const LO* lclColInds = NULL;
        Scalar* vals = NULL;
        LO numInds = 0;
        this->getLocalRowView (meshLclRow, lclColInds, vals, numInds);

        for (LO k = 0; k < numInds; ++k) {
          const GO gblCol = meshColMap.getGlobalElement (lclColInds[k]);

          os << "Col " << gblCol << ": [";
          for (LO i = 0; i < blockSize; ++i) {
            for (LO j = 0; j < blockSize; ++j) {
              os << vals[blockSize*blockSize*k + i*blockSize + j];
              if (j + 1 < blockSize) {
                os << ", ";
              }
            }
            if (i + 1 < blockSize) {
              os << "; ";
            }
          }
          os << "]";
          if (k + 1 < numInds) {
            os << ", ";
          }
        }
        os << "}" << endl;
      }

      // Print data on Process 0.  This will automatically respect the
      // current indentation level.
      if (myRank == 0) {
        out << lclOutStrPtr->str ();
        lclOutStrPtr = Teuchos::null; // clear it to save space
      }

      const int sizeTag = 1337;
      const int dataTag = 1338;

      ArrayRCP<char> recvDataBuf; // only used on Process 0

      // Send string sizes and data from each process in turn to
      // Process 0, and print on that process.
      for (int p = 1; p < numProcs; ++p) {
        if (myRank == 0) {
          // Receive the incoming string's length.
          ArrayRCP<size_t> recvSize (1);
          recvSize[0] = 0;
          RCP<CommRequest<int> > recvSizeReq =
            ireceive<int, size_t> (recvSize, p, sizeTag, comm);
          wait<int> (comm, outArg (recvSizeReq));
          const size_t numCharsToRecv = recvSize[0];

          // Allocate space for the string to receive.  Reuse receive
          // buffer space if possible.  We can do this because in the
          // current implementation, we only have one receive in
          // flight at a time.  Leave space for the '\0' at the end,
          // in case the sender doesn't send it.
          if (static_cast<size_t>(recvDataBuf.size()) < numCharsToRecv + 1) {
            recvDataBuf.resize (numCharsToRecv + 1);
          }
          ArrayRCP<char> recvData = recvDataBuf.persistingView (0, numCharsToRecv);
          // Post the receive of the actual string data.
          RCP<CommRequest<int> > recvDataReq =
            ireceive<int, char> (recvData, p, dataTag, comm);
          wait<int> (comm, outArg (recvDataReq));

          // Print the received data.  This will respect the current
          // indentation level.  Make sure that the string is
          // null-terminated.
          recvDataBuf[numCharsToRecv] = '\0';
          out << recvDataBuf.getRawPtr ();
        }
        else if (myRank == p) { // if I am not Process 0, and my rank is p
          // This deep-copies the string at most twice, depending on
          // whether std::string reference counts internally (it
          // generally does, so this won't deep-copy at all).
          const std::string stringToSend = lclOutStrPtr->str ();
          lclOutStrPtr = Teuchos::null; // clear original to save space

          // Send the string's length to Process 0.
          const size_t numCharsToSend = stringToSend.size ();
          ArrayRCP<size_t> sendSize (1);
          sendSize[0] = numCharsToSend;
          RCP<CommRequest<int> > sendSizeReq =
            isend<int, size_t> (sendSize, 0, sizeTag, comm);
          wait<int> (comm, outArg (sendSizeReq));

          // Send the actual string to Process 0.  We know that the
          // string has length > 0, so it's save to take the address
          // of the first entry.  Make a nonowning ArrayRCP to hold
          // the string.  Process 0 will add a null termination
          // character at the end of the string, after it receives the
          // message.
          ArrayRCP<const char> sendData (&stringToSend[0], 0, numCharsToSend, false);
          RCP<CommRequest<int> > sendDataReq =
            isend<int, char> (sendData, 0, dataTag, comm);
          wait<int> (comm, outArg (sendDataReq));
        }
      } // for each process rank p other than 0
    } // extreme verbosity level (print the whole matrix)
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const Teuchos::Comm<int> >
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getComm() const
  {
    return graph_.getComm();
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<Node>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNode() const
  {
    return graph_.getNode();

  }

  template<class Scalar, class LO, class GO, class Node>
  global_size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalNumCols() const
  {
    return graph_.getGlobalNumCols();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNodeNumCols() const
  {
    return graph_.getNodeNumCols();
  }

  template<class Scalar, class LO, class GO, class Node>
  GO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getIndexBase() const
  {
    return graph_.getIndexBase();
  }

  template<class Scalar, class LO, class GO, class Node>
  global_size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalNumEntries() const
  {
    return graph_.getGlobalNumEntries();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNodeNumEntries() const
  {
    return graph_.getNodeNumEntries();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNumEntriesInGlobalRow (GO globalRow) const
  {
    return graph_.getNumEntriesInGlobalRow(globalRow);
  }

  template<class Scalar, class LO, class GO, class Node>
  global_size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalNumDiags() const
  {
    using HDM = Details::HasDeprecatedMethods2630_WarningThisClassIsNotForUsers;
    return dynamic_cast<const HDM&> (this->graph_).getGlobalNumDiagsImpl ();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNodeNumDiags() const
  {
    using HDM = Details::HasDeprecatedMethods2630_WarningThisClassIsNotForUsers;
    return dynamic_cast<const HDM&> (this->graph_).getNodeNumDiagsImpl ();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalMaxNumRowEntries() const
  {
    return graph_.getGlobalMaxNumRowEntries();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  hasColMap() const
  {
    return graph_.hasColMap();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  isLowerTriangular () const
  {
    using HDM = ::Tpetra::Details::HasDeprecatedMethods2630_WarningThisClassIsNotForUsers;
    return dynamic_cast<const HDM&> (this->graph_).isLowerTriangularImpl ();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  isUpperTriangular () const
  {
    using HDM = ::Tpetra::Details::HasDeprecatedMethods2630_WarningThisClassIsNotForUsers;
    return dynamic_cast<const HDM&> (this->graph_).isUpperTriangularImpl ();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  isLocallyIndexed() const
  {
    return graph_.isLocallyIndexed();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  isGloballyIndexed() const
  {
    return graph_.isGloballyIndexed();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  isFillComplete() const
  {
    return graph_.isFillComplete ();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  supportsRowViews() const
  {
    return false;
  }


  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalRowCopy (GO GlobalRow,
                    const Teuchos::ArrayView<GO> &Indices,
                    const Teuchos::ArrayView<Scalar> &Values,
                    size_t &NumEntries) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::getGlobalRowCopy: "
      "This class doesn't support global matrix indexing.");

  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalRowView (GO GlobalRow,
                    Teuchos::ArrayView<const GO> &indices,
                    Teuchos::ArrayView<const Scalar> &values) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::getGlobalRowView: "
      "This class doesn't support global matrix indexing.");

  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalRowView (LO LocalRow,
                   Teuchos::ArrayView<const LO>& indices,
                   Teuchos::ArrayView<const Scalar>& values) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::getLocalRowView: "
      "This class doesn't support local matrix indexing.");
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalDiagCopy (::Tpetra::Vector<Scalar,LO,GO,Node>& diag) const
  {
#ifdef HAVE_TPETRA_DEBUG
    const char prefix[] =
      "Tpetra::Experimental::BlockCrsMatrix::getLocalDiagCopy: ";
#endif // HAVE_TPETRA_DEBUG

    const size_t lclNumMeshRows = graph_.getNodeNumRows ();

    Kokkos::View<size_t*, device_type> diagOffsets ("diagOffsets", lclNumMeshRows);
    graph_.getLocalDiagOffsets (diagOffsets);

    // The code below works on host, so use a host View.
    auto diagOffsetsHost = Kokkos::create_mirror_view (diagOffsets);
    Kokkos::deep_copy (diagOffsetsHost, diagOffsets);
    // We're filling diag on host for now.
    diag.template modify<typename decltype (diagOffsetsHost)::memory_space> ();

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION
      (this->need_sync_host (), std::runtime_error,
       prefix << "The matrix's data were last modified on device, but have "
       "not been sync'd to host.  Please sync to host (by calling "
       "sync<Kokkos::HostSpace>() on this matrix) before calling this "
       "method.");
#endif // HAVE_TPETRA_DEBUG

    auto vals_host_out = getValuesHost ();
    Scalar* vals_host_out_raw =
      reinterpret_cast<Scalar*> (vals_host_out.data ());

    // TODO amk: This is a temporary measure to make the code run with Ifpack2
    size_t rowOffset = 0;
    size_t offset = 0;
    LO bs = getBlockSize();
    for(size_t r=0; r<getNodeNumRows(); r++)
    {
      // move pointer to start of diagonal block
      offset = rowOffset + diagOffsetsHost(r)*bs*bs;
      for(int b=0; b<bs; b++)
      {
        diag.replaceLocalValue(r*bs+b, vals_host_out_raw[offset+b*(bs+1)]);
      }
      // move pointer to start of next block row
      rowOffset += getNumEntriesInLocalRow(r)*bs*bs;
    }

    diag.template sync<memory_space> (); // sync vec of diag entries back to dev
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  leftScale (const ::Tpetra::Vector<Scalar, LO, GO, Node>& x)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::leftScale: "
      "not implemented.");

  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  rightScale (const ::Tpetra::Vector<Scalar, LO, GO, Node>& x)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::rightScale: "
      "not implemented.");

  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const ::Tpetra::RowGraph<LO, GO, Node> >
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGraph() const
  {
    return graphRCP_;
  }

  template<class Scalar, class LO, class GO, class Node>
  typename ::Tpetra::RowMatrix<Scalar, LO, GO, Node>::mag_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getFrobeniusNorm () const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::getFrobeniusNorm: "
      "not implemented.");
  }

} // namespace Experimental
} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_EXPERIMENTAL_BLOCKCRSMATRIX_INSTANT(S,LO,GO,NODE) \
  namespace Experimental { \
    template class BlockCrsMatrix< S, LO, GO, NODE >; \
  }

#endif // TPETRA_EXPERIMENTAL_BLOCKCRSMATRIX_DEF_HPP
