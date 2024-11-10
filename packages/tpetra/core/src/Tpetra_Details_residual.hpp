// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_RESIDUAL_HPP
#define TPETRA_DETAILS_RESIDUAL_HPP

#include "TpetraCore_config.h"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_LocalCrsMatrixOperator.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_Profiling.hpp"
#include "KokkosSparse_spmv_impl.hpp"

/// \file Tpetra_Details_residual.hpp
/// \brief Functions that allow for fused residual calculation.
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.-

namespace Tpetra {
  namespace Details {


/// \brief Functor for computing the residual
///
/// The template parameters decide what is actually computed.
/// is_MV:          whether the multivectors have actually more than 1 vector.
/// restrictedMode: when true, read on-rank RHS data out of X_domainmap_lcl instead of X_colmap_lcl
/// skipOffRank:    when true, only compute the local part of the residual
template<class AMatrix, class MV, class ConstMV, class Offsets, bool is_MV, bool restrictedMode, bool skipOffRank>
struct LocalResidualFunctor {

  using execution_space = typename AMatrix::execution_space;
  using LO = typename AMatrix::non_const_ordinal_type;
  using value_type = typename AMatrix::non_const_value_type;
  using team_policy = typename Kokkos::TeamPolicy<execution_space>;
  using team_member = typename team_policy::member_type;
  using ATV = Kokkos::ArithTraits<value_type>;

  AMatrix A_lcl;
  ConstMV X_colmap_lcl;
  ConstMV B_lcl;
  MV R_lcl;
  int rows_per_team;
  Offsets offsets;
  ConstMV X_domainmap_lcl;


  LocalResidualFunctor (const AMatrix& A_lcl_,
                        const ConstMV& X_colmap_lcl_,
                        const ConstMV& B_lcl_,
                        const MV& R_lcl_,
                        const int rows_per_team_,
                        const Offsets& offsets_,
                        const ConstMV& X_domainmap_lcl_) :
    A_lcl(A_lcl_),
    X_colmap_lcl(X_colmap_lcl_),
    B_lcl(B_lcl_),
    R_lcl(R_lcl_),
    rows_per_team(rows_per_team_),
    offsets(offsets_),
    X_domainmap_lcl(X_domainmap_lcl_)
  { }

  KOKKOS_INLINE_FUNCTION
  void operator() (const team_member& dev) const
  {

    Kokkos::parallel_for(Kokkos::TeamThreadRange (dev, 0, rows_per_team),[&] (const LO& loop) {
      const LO lclRow = static_cast<LO> (dev.league_rank ()) * rows_per_team + loop;

      if (lclRow >= A_lcl.numRows ()) {
        return;
      }

      if (!is_MV) { // MultiVectors only have a single column

        value_type A_x = ATV::zero ();

        if (!restrictedMode) {
          const auto A_row = A_lcl.rowConst(lclRow);
          const LO row_length = static_cast<LO> (A_row.length);

          Kokkos::parallel_reduce(Kokkos::ThreadVectorRange (dev, row_length), [&] (const LO iEntry, value_type& lsum) {
            const auto A_val = A_row.value(iEntry);
            lsum += A_val * X_colmap_lcl(A_row.colidx(iEntry),0);
          }, A_x);

        }
        else {

          const LO offRankOffset = offsets(lclRow);
          const size_t start = A_lcl.graph.row_map(lclRow);
          const size_t end = A_lcl.graph.row_map(lclRow+1);

          Kokkos::parallel_reduce(Kokkos::ThreadVectorRange (dev, start, end), [&] (const LO iEntry, value_type& lsum) {
            const auto A_val = A_lcl.values(iEntry);
            const auto lclCol = A_lcl.graph.entries(iEntry);
            if (iEntry < offRankOffset)
              lsum += A_val * X_domainmap_lcl(lclCol,0);
            else if (!skipOffRank)
              lsum += A_val * X_colmap_lcl(lclCol,0);
          }, A_x);
        }

        Kokkos::single(Kokkos::PerThread(dev),[&] () {
          R_lcl(lclRow,0) = B_lcl(lclRow,0) - A_x;
        });
      }
      else { // MultiVectors have more than one column

        const LO numVectors = static_cast<LO>(X_colmap_lcl.extent(1));

        for(LO v=0; v<numVectors; v++) {

          value_type A_x = ATV::zero ();

          if (!restrictedMode) {

            const auto A_row = A_lcl.rowConst(lclRow);
            const LO row_length = static_cast<LO> (A_row.length);

            Kokkos::parallel_reduce(Kokkos::ThreadVectorRange (dev, row_length), [&] (const LO iEntry, value_type& lsum) {
              const auto A_val = A_row.value(iEntry);
              lsum += A_val * X_colmap_lcl(A_row.colidx(iEntry),v);
            }, A_x);
          }
          else {
            const LO offRankOffset = offsets(lclRow);
            const size_t start = A_lcl.graph.row_map(lclRow);
            const size_t end = A_lcl.graph.row_map(lclRow+1);

            Kokkos::parallel_reduce(Kokkos::ThreadVectorRange (dev, start, end), [&] (const LO iEntry, value_type& lsum) {
              const auto A_val = A_lcl.values(iEntry);
              const auto lclCol = A_lcl.graph.entries(iEntry);
              if (iEntry < offRankOffset)
                lsum += A_val * X_domainmap_lcl(lclCol,v);
              else if (!skipOffRank)
                lsum += A_val * X_colmap_lcl(lclCol,v);
            }, A_x);
          }

          Kokkos::single(Kokkos::PerThread(dev),[&] () {
            R_lcl(lclRow,v) = B_lcl(lclRow,v) - A_x;
          });

        }//end for numVectors
      }
    });//end parallel_for TeamThreadRange
  }
};


/// \brief Functor for computing R -= A_offRank*X_colmap
template<class AMatrix, class MV, class ConstMV, class Offsets, bool is_MV>
struct OffRankUpdateFunctor {

  using execution_space = typename AMatrix::execution_space;
  using LO = typename AMatrix::non_const_ordinal_type;
  using value_type = typename AMatrix::non_const_value_type;
  using team_policy = typename Kokkos::TeamPolicy<execution_space>;
  using team_member = typename team_policy::member_type;
  using ATV = Kokkos::ArithTraits<value_type>;

  AMatrix A_lcl;
  ConstMV X_colmap_lcl;
  MV R_lcl;
  int rows_per_team;
  Offsets offsets;


  OffRankUpdateFunctor (const AMatrix& A_lcl_,
                        const ConstMV& X_colmap_lcl_,
                        const MV& R_lcl_,
                        const int rows_per_team_,
                        const Offsets& offsets_) :
    A_lcl(A_lcl_),
    X_colmap_lcl(X_colmap_lcl_),
    R_lcl(R_lcl_),
    rows_per_team(rows_per_team_),
    offsets(offsets_)
  { }

  KOKKOS_INLINE_FUNCTION
  void operator() (const team_member& dev) const
  {

    Kokkos::parallel_for(Kokkos::TeamThreadRange (dev, 0, rows_per_team),[&] (const LO& loop) {
      const LO lclRow = static_cast<LO> (dev.league_rank ()) * rows_per_team + loop;

      if (lclRow >= A_lcl.numRows ()) {
        return;
      }

      const LO offRankOffset = offsets(lclRow);
      const size_t end = A_lcl.graph.row_map(lclRow+1);

      if (!is_MV) { // MultiVectors only have a single column

        value_type A_x = ATV::zero ();

        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange (dev, offRankOffset, end), [&] (const LO iEntry, value_type& lsum) {
          const auto A_val = A_lcl.values(iEntry);
          const auto lclCol = A_lcl.graph.entries(iEntry);
          lsum += A_val * X_colmap_lcl(lclCol,0);
        }, A_x);

        Kokkos::single(Kokkos::PerThread(dev),[&] () {
          R_lcl(lclRow,0) -=  A_x;
        });
      }
      else { // MultiVectors have more than one column

        const LO numVectors = static_cast<LO>(X_colmap_lcl.extent(1));

        for(LO v=0; v<numVectors; v++) {

          value_type A_x = ATV::zero ();

          Kokkos::parallel_reduce(Kokkos::ThreadVectorRange (dev, offRankOffset, end), [&] (const LO iEntry, value_type& lsum) {
            const auto A_val = A_lcl.values(iEntry);
            const auto lclCol = A_lcl.graph.entries(iEntry);
            lsum += A_val * X_colmap_lcl(lclCol,v);
          }, A_x);

          Kokkos::single(Kokkos::PerThread(dev),[&] () {
            R_lcl(lclRow,v) -= A_x;
          });

        }//end for numVectors
      }
    });
  }
};

template<class SC, class LO, class GO, class NO>
void localResidual(const CrsMatrix<SC,LO,GO,NO> &  A,
                   const MultiVector<SC,LO,GO,NO> & X_colmap,
                   const MultiVector<SC,LO,GO,NO> & B,
                   MultiVector<SC,LO,GO,NO> & R,
                   const Kokkos::View<const size_t*, typename NO::device_type>& offsets,
                   const MultiVector<SC,LO,GO,NO> * X_domainmap=nullptr) {
  using Tpetra::Details::ProfilingRegion;
  using Teuchos::NO_TRANS;
  ProfilingRegion regionLocalApply ("Tpetra::CrsMatrix::localResidual");

  using local_matrix_device_type = typename CrsMatrix<SC,LO,GO,NO>::local_matrix_device_type;
  using local_view_device_type = typename  MultiVector<SC,LO,GO,NO>::dual_view_type::t_dev::non_const_type;
  using const_local_view_device_type = typename  MultiVector<SC,LO,GO,NO>::dual_view_type::t_dev::const_type;
  using offset_type = Kokkos::View<const size_t*, typename NO::device_type>;

  local_matrix_device_type     A_lcl = A.getLocalMatrixDevice ();
  const_local_view_device_type X_colmap_lcl = X_colmap.getLocalViewDevice(Access::ReadOnly);
  const_local_view_device_type B_lcl = B.getLocalViewDevice(Access::ReadOnly);
  local_view_device_type       R_lcl = R.getLocalViewDevice(Access::OverwriteAll);
  const_local_view_device_type X_domainmap_lcl;

  const bool debug = ::Tpetra::Details::Behavior::debug ();
  if (debug) {
    TEUCHOS_TEST_FOR_EXCEPTION
      (X_colmap.getNumVectors () != R.getNumVectors (), std::runtime_error,
       "X.getNumVectors() = " << X_colmap.getNumVectors () << " != "
       "R.getNumVectors() = " << R.getNumVectors () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (X_colmap.getLocalLength () !=
       A.getColMap ()->getLocalNumElements (), std::runtime_error,
       "X has the wrong number of local rows.  "
       "X.getLocalLength() = " << X_colmap.getLocalLength () << " != "
       "A.getColMap()->getLocalNumElements() = " <<
       A.getColMap ()->getLocalNumElements () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (R.getLocalLength () !=
       A.getRowMap ()->getLocalNumElements (), std::runtime_error,
       "R has the wrong number of local rows.  "
       "R.getLocalLength() = " << R.getLocalLength () << " != "
       "A.getRowMap()->getLocalNumElements() = " <<
       A.getRowMap ()->getLocalNumElements () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (B.getLocalLength () !=
       A.getRowMap ()->getLocalNumElements (), std::runtime_error,
       "B has the wrong number of local rows.  "
       "B.getLocalLength() = " << B.getLocalLength () << " != "
       "A.getRowMap()->getLocalNumElements() = " <<
       A.getRowMap ()->getLocalNumElements () << ".");

    TEUCHOS_TEST_FOR_EXCEPTION
      (! A.isFillComplete (), std::runtime_error, "The matrix A is not "
       "fill complete.  You must call fillComplete() (possibly with "
       "domain and range Map arguments) without an intervening "
       "resumeFill() call before you may call this method.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (! X_colmap.isConstantStride () || ! R.isConstantStride () || ! B.isConstantStride (),
       std::runtime_error, "X, Y and B must be constant stride.");
    // If the two pointers are NULL, then they don't alias one
    // another, even though they are equal.
    TEUCHOS_TEST_FOR_EXCEPTION
      ((X_colmap_lcl.data () == R_lcl.data () && X_colmap_lcl.data () != nullptr) ||
       (X_colmap_lcl.data () == B_lcl.data () && X_colmap_lcl.data () != nullptr),
       std::runtime_error, "X, Y and R may not alias one another.");
  }

  const bool fusedResidual = ::Tpetra::Details::Behavior::fusedResidual ();
  if (!fusedResidual) {
    SC one = Teuchos::ScalarTraits<SC>::one();
    SC negone = -one;
    SC zero = Teuchos::ScalarTraits<SC>::zero();
    // This is currently a "reference implementation" waiting until Kokkos Kernels provides
    // a residual kernel.
    A.localApply(X_colmap,R,Teuchos::NO_TRANS, one, zero);
    R.update(one,B,negone);
    return;
  }

  if (A_lcl.numRows() == 0) {
    return;
  }

  int64_t numLocalRows = A_lcl.numRows ();
  int64_t myNnz = A_lcl.nnz();

  int team_size = -1;
  int vector_length = -1;
  int64_t rows_per_thread = -1;

  using execution_space = typename CrsMatrix<SC,LO,GO,NO>::execution_space;
  using policy_type = typename Kokkos::TeamPolicy<execution_space>;

  int64_t rows_per_team = KokkosSparse::Impl::spmv_launch_parameters<execution_space>(numLocalRows, myNnz, rows_per_thread, team_size, vector_length);
  int64_t worksets = (B_lcl.extent (0) + rows_per_team - 1) / rows_per_team;

  policy_type policy (1, 1);
  if (team_size < 0) {
    policy = policy_type (worksets, Kokkos::AUTO, vector_length);
  }
  else {
    policy = policy_type (worksets, team_size, vector_length);
  }

  bool is_vector = (X_colmap_lcl.extent(1) == 1);

  if(is_vector) {

    if (X_domainmap == nullptr) {

      using functor_type = LocalResidualFunctor<local_matrix_device_type,local_view_device_type,const_local_view_device_type,offset_type,false,false,false>;
      functor_type func (A_lcl, X_colmap_lcl, B_lcl, R_lcl, rows_per_team, offsets, X_domainmap_lcl);
      Kokkos::parallel_for("residual-vector",policy,func);

    }
    else {

      X_domainmap_lcl = X_domainmap->getLocalViewDevice(Access::ReadOnly);
      using functor_type = LocalResidualFunctor<local_matrix_device_type,local_view_device_type,const_local_view_device_type,offset_type,false,true,false>;
      functor_type func (A_lcl, X_colmap_lcl, B_lcl, R_lcl, rows_per_team, offsets, X_domainmap_lcl);
      Kokkos::parallel_for("residual-vector",policy,func);

    }
  }
  else {

    if (X_domainmap == nullptr) {

      using functor_type = LocalResidualFunctor<local_matrix_device_type,local_view_device_type,const_local_view_device_type,offset_type,true,false,false>;
      functor_type func (A_lcl, X_colmap_lcl, B_lcl, R_lcl, rows_per_team, offsets, X_domainmap_lcl);
      Kokkos::parallel_for("residual-multivector",policy,func);

    }
    else {

      X_domainmap_lcl = X_domainmap->getLocalViewDevice(Access::ReadOnly);
      using functor_type = LocalResidualFunctor<local_matrix_device_type,local_view_device_type,const_local_view_device_type,offset_type,true,true,false>;
      functor_type func (A_lcl, X_colmap_lcl, B_lcl, R_lcl, rows_per_team, offsets, X_domainmap_lcl);
      Kokkos::parallel_for("residual-multivector",policy,func);

    }
  }
}


template<class SC, class LO, class GO, class NO>
void localResidualWithCommCompOverlap(const CrsMatrix<SC,LO,GO,NO> &  A,
                                      MultiVector<SC,LO,GO,NO> & X_colmap,
                                      const MultiVector<SC,LO,GO,NO> & B,
                                      MultiVector<SC,LO,GO,NO> & R,
                                      const Kokkos::View<const size_t*, typename NO::device_type>& offsets,
                                      const MultiVector<SC,LO,GO,NO> & X_domainmap) {
  using Tpetra::Details::ProfilingRegion;
  using Teuchos::NO_TRANS;
  using Teuchos::RCP;
  using import_type = typename CrsMatrix<SC,LO,GO,NO>::import_type;

  ProfilingRegion regionLocalApply ("Tpetra::CrsMatrix::localResidualWithCommCompOverlap");

  using local_matrix_device_type = typename CrsMatrix<SC,LO,GO,NO>::local_matrix_device_type;
  using local_view_device_type = typename  MultiVector<SC,LO,GO,NO>::dual_view_type::t_dev::non_const_type;
  using const_local_view_device_type = typename  MultiVector<SC,LO,GO,NO>::dual_view_type::t_dev::const_type;
  using offset_type = Kokkos::View<const size_t*, typename NO::device_type>;

  local_matrix_device_type     A_lcl = A.getLocalMatrixDevice ();
  const_local_view_device_type X_colmap_lcl = X_colmap.getLocalViewDevice(Access::ReadOnly);
  const_local_view_device_type B_lcl = B.getLocalViewDevice(Access::ReadOnly);
  local_view_device_type       R_lcl = R.getLocalViewDevice(Access::OverwriteAll);
  const_local_view_device_type X_domainmap_lcl = X_domainmap.getLocalViewDevice(Access::ReadOnly);;

  const bool debug = ::Tpetra::Details::Behavior::debug ();
  if (debug) {
    TEUCHOS_TEST_FOR_EXCEPTION
      (X_colmap.getNumVectors () != R.getNumVectors (), std::runtime_error,
       "X.getNumVectors() = " << X_colmap.getNumVectors () << " != "
       "R.getNumVectors() = " << R.getNumVectors () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (X_colmap.getLocalLength () !=
       A.getColMap ()->getLocalNumElements (), std::runtime_error,
       "X has the wrong number of local rows.  "
       "X.getLocalLength() = " << X_colmap.getLocalLength () << " != "
       "A.getColMap()->getLocalNumElements() = " <<
       A.getColMap ()->getLocalNumElements () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (R.getLocalLength () !=
       A.getRowMap ()->getLocalNumElements (), std::runtime_error,
       "R has the wrong number of local rows.  "
       "R.getLocalLength() = " << R.getLocalLength () << " != "
       "A.getRowMap()->getLocalNumElements() = " <<
       A.getRowMap ()->getLocalNumElements () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (B.getLocalLength () !=
       A.getRowMap ()->getLocalNumElements (), std::runtime_error,
       "B has the wrong number of local rows.  "
       "B.getLocalLength() = " << B.getLocalLength () << " != "
       "A.getRowMap()->getLocalNumElements() = " <<
       A.getRowMap ()->getLocalNumElements () << ".");

    TEUCHOS_TEST_FOR_EXCEPTION
      (! A.isFillComplete (), std::runtime_error, "The matrix A is not "
       "fill complete.  You must call fillComplete() (possibly with "
       "domain and range Map arguments) without an intervening "
       "resumeFill() call before you may call this method.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (! X_colmap.isConstantStride () || ! R.isConstantStride () || ! B.isConstantStride (),
       std::runtime_error, "X, Y and B must be constant stride.");
    // If the two pointers are NULL, then they don't alias one
    // another, even though they are equal.
    TEUCHOS_TEST_FOR_EXCEPTION
      ((X_colmap_lcl.data () == R_lcl.data () && X_colmap_lcl.data () != nullptr) ||
       (X_colmap_lcl.data () == B_lcl.data () && X_colmap_lcl.data () != nullptr),
       std::runtime_error, "X, Y and R may not alias one another.");
  }

  if (A_lcl.numRows() == 0) {
    return;
  }

  int64_t numLocalRows = A_lcl.numRows ();
  int64_t myNnz = A_lcl.nnz();

  int team_size = -1;
  int vector_length = -1;
  int64_t rows_per_thread = -1;

  using execution_space = typename CrsMatrix<SC,LO,GO,NO>::execution_space;
  using policy_type = typename Kokkos::TeamPolicy<execution_space>;

  int64_t rows_per_team = KokkosSparse::Impl::spmv_launch_parameters<execution_space>(numLocalRows, myNnz, rows_per_thread, team_size, vector_length);
  int64_t worksets = (B_lcl.extent (0) + rows_per_team - 1) / rows_per_team;

  policy_type policy (1, 1);
  if (team_size < 0) {
    policy = policy_type (worksets, Kokkos::AUTO, vector_length);
  }
  else {
    policy = policy_type (worksets, team_size, vector_length);
  }

  bool is_vector = (X_colmap_lcl.extent(1) == 1);

  if(is_vector) {

    using functor_type = LocalResidualFunctor<local_matrix_device_type,local_view_device_type,const_local_view_device_type,offset_type,false,true,true>;
    functor_type func (A_lcl, X_colmap_lcl, B_lcl, R_lcl, rows_per_team, offsets, X_domainmap_lcl);
    Kokkos::parallel_for("residual-vector",policy,func);

    RCP<const import_type> importer = A.getGraph ()->getImporter ();
    X_colmap.endImport (X_domainmap, *importer, INSERT, true);

    Kokkos::fence("Tpetra::localResidualWithCommCompOverlap-1");

    using functor_type2 = OffRankUpdateFunctor<local_matrix_device_type,local_view_device_type,const_local_view_device_type,offset_type,false>;
    functor_type2 func2 (A_lcl, X_colmap_lcl, R_lcl, rows_per_team, offsets);
    Kokkos::parallel_for("residual-vector-offrank",policy,func2);

  }
  else {

    using functor_type = LocalResidualFunctor<local_matrix_device_type,local_view_device_type,const_local_view_device_type,offset_type,true,true,true>;
    functor_type func (A_lcl, X_colmap_lcl, B_lcl, R_lcl, rows_per_team, offsets, X_domainmap_lcl);
    Kokkos::parallel_for("residual-multivector",policy,func);

    RCP<const import_type> importer = A.getGraph ()->getImporter ();
    X_colmap.endImport (X_domainmap, *importer, INSERT, true);

    Kokkos::fence("Tpetra::localResidualWithCommCompOverlap-2");

    using functor_type2 = OffRankUpdateFunctor<local_matrix_device_type,local_view_device_type,const_local_view_device_type,offset_type,true>;
    functor_type2 func2 (A_lcl, X_colmap_lcl, R_lcl, rows_per_team, offsets);
    Kokkos::parallel_for("residual-vector-offrank",policy,func2);

  }
}


//! Computes R = B - A * X 
template<class SC, class LO, class GO, class NO>
void residual(const Operator<SC,LO,GO,NO> &   Aop,
              const MultiVector<SC,LO,GO,NO> & X_in,
              const MultiVector<SC,LO,GO,NO> & B_in,
              MultiVector<SC,LO,GO,NO> & R_in) {
  using Tpetra::Details::ProfilingRegion;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcpFromRef;

  const bool debug = ::Tpetra::Details::Behavior::debug ();
  const bool skipCopyAndPermuteIfPossible = ::Tpetra::Details::Behavior::skipCopyAndPermuteIfPossible ();
  const bool overlapCommunicationAndComputation = ::Tpetra::Details::Behavior::overlapCommunicationAndComputation ();
  if (overlapCommunicationAndComputation)
    TEUCHOS_ASSERT(skipCopyAndPermuteIfPossible);

  // Whether we are using restrictedMode in the import from domain to
  // column map. Restricted mode skips the copy and permutation of the
  // local part of X. We are using restrictedMode only when domain and
  // column map are locally fitted, i.e. when the local indices of
  // domain and column map match.
  bool restrictedMode = false;
  
  const CrsMatrix<SC,LO,GO,NO> * Apt  = dynamic_cast<const CrsMatrix<SC,LO,GO,NO>*>(&Aop);
  if(!Apt) {
    // If we're not a CrsMatrix, we can't do fusion, so just do apply+update
     SC one = Teuchos::ScalarTraits<SC>::one(), negone = -one, zero = Teuchos::ScalarTraits<SC>::zero();    
     Aop.apply(X_in,R_in,Teuchos::NO_TRANS, one, zero);
     R_in.update(one,B_in,negone);
     return;
  }
  const CrsMatrix<SC,LO,GO,NO> & A = *Apt;

  using import_type = typename CrsMatrix<SC,LO,GO,NO>::import_type;
  using export_type = typename CrsMatrix<SC,LO,GO,NO>::export_type;
  using MV = MultiVector<SC,LO,GO,NO>;
  using graph_type = Tpetra::CrsGraph<LO,GO,NO>;
  using offset_type = typename graph_type::offset_device_view_type;

  // We treat the case of a replicated MV output specially.
  const bool R_is_replicated =
    (! R_in.isDistributed () && A.getComm ()->getSize () != 1);

  // It's possible that R is a view of X or B.  
  // We don't try to to detect the more subtle cases (e.g., one is a
  // subview of the other, but their initial pointers differ).  We
  // only need to do this if this matrix's Import is trivial;
  // otherwise, we don't actually apply the operator from X into Y.
  
  RCP<const import_type> importer = A.getGraph ()->getImporter ();
  RCP<const export_type> exporter = A.getGraph ()->getExporter ();

  // Temporary MV for Import operation.  After the block of code
  // below, this will be an (Imported if necessary) column Map MV
  // ready to give to localApply(...).
  RCP<MV> X_colMap;
  if (importer.is_null ()) {
    if (! X_in.isConstantStride ()) {
      // Not all sparse mat-vec kernels can handle an input MV with
      // nonconstant stride correctly, so we have to copy it in that
      // case into a constant stride MV.  To make a constant stride
      // copy of X_in, we force creation of the column (== domain)
      // Map MV (if it hasn't already been created, else fetch the
      // cached copy).  This avoids creating a new MV each time.
      X_colMap = A.getColumnMapMultiVector (X_in, true);
      Tpetra::deep_copy (*X_colMap, X_in);
      // X_colMap = rcp_const_cast<const MV> (X_colMapNonConst);
    }
    else {
      // The domain and column Maps are the same, so do the local
      // multiply using the domain Map input MV X_in.
      X_colMap = rcp_const_cast<MV> (rcpFromRef (X_in) );
    }
  }
  else { // need to Import source (multi)vector
    ProfilingRegion regionImport ("Tpetra::CrsMatrix::residual: Import");
    // We're doing an Import anyway, which will copy the relevant
    // elements of the domain Map MV X_in into a separate column Map
    // MV.  Thus, we don't have to worry whether X_in is constant
    // stride.
    X_colMap = A.getColumnMapMultiVector (X_in);
    
    // Do we want to use restrictedMode?
    restrictedMode = skipCopyAndPermuteIfPossible && importer->isLocallyFitted();

    if (debug && restrictedMode) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (!importer->getTargetMap()->isLocallyFitted(*importer->getSourceMap()), std::runtime_error,
         "Source map and target map are not locally fitted, but Tpetra::residual thinks they are.");
    }

    // Import from the domain Map MV to the column Map MV.
    X_colMap->beginImport (X_in, *importer, INSERT, restrictedMode);
  }

  offset_type offsets;
  if (restrictedMode)
    A.getCrsGraph()->getLocalOffRankOffsets(offsets);

  // Get a vector for the R_rowMap output residual, handling the 
  // non-constant stride and exporter cases.  Since R gets clobbered
  // we don't need to worry about the data in it
  RCP<MV> R_rowMap;
  if(exporter.is_null()) {
    if (! R_in.isConstantStride ()) {
      R_rowMap = A.getRowMapMultiVector(R_in);
    }
    else {
      R_rowMap = rcpFromRef (R_in);
    }
  }
  else {
    R_rowMap = A.getRowMapMultiVector (R_in);    
  }
  
  // Get a vector for the B_rowMap output residual, handling the 
  // non-constant stride and exporter cases
  RCP<const MV> B_rowMap;
  if(exporter.is_null()) {
    if (! B_in.isConstantStride ()) {
      // Do an allocation here.  If we need to optimize this later, we can have the matrix 
      // cache this.
      RCP<MV> B_rowMapNonConst = rcp(new MV(A.getRowMap(),B_in.getNumVectors()));
      Tpetra::deep_copy (*B_rowMapNonConst, B_in);
      B_rowMap = rcp_const_cast<const MV> (B_rowMapNonConst);
    }
    else {
      B_rowMap = rcpFromRef (B_in);
    }
  }
  else {
    // Do an allocation here.  If we need to optimize this later, we can have the matrix 
    // cache this.
    ProfilingRegion regionExport ("Tpetra::CrsMatrix::residual: B Import");
    RCP<MV> B_rowMapNonConst = rcp(new MV(A.getRowMap(),B_in.getNumVectors()));
    B_rowMapNonConst->doImport(B_in, *exporter, ADD);
    B_rowMap = rcp_const_cast<const MV> (B_rowMapNonConst);
  }

  // If we have a nontrivial Export object, we must perform an
  // Export.  In that case, the local multiply result will go into
  // the row Map multivector.  We don't have to make a
  // constant-stride version of R_in in this case, because we had to
  // make a constant stride R_rowMap MV and do an Export anyway.
  if (! exporter.is_null ()) {
    if ( ! importer.is_null ())
      X_colMap->endImport (X_in, *importer, INSERT, restrictedMode);
    if (restrictedMode && !importer.is_null ())
      localResidual (A, *X_colMap, *B_rowMap, *R_rowMap, offsets, &X_in);
    else
      localResidual (A, *X_colMap, *B_rowMap, *R_rowMap, offsets);
    
    {
      ProfilingRegion regionExport ("Tpetra::CrsMatrix::residual: R Export");
      
      // Do the Export operation.
      R_in.doExport (*R_rowMap, *exporter, ADD);
    }
  }
  else { // Don't do an Export: row Map and range Map are the same.

    if (restrictedMode)
      if (overlapCommunicationAndComputation) {
        localResidualWithCommCompOverlap (A, *X_colMap, *B_rowMap, *R_rowMap, offsets, X_in);
      } else {
        X_colMap->endImport (X_in, *importer, INSERT, restrictedMode);
        localResidual (A, *X_colMap, *B_rowMap, *R_rowMap, offsets, &X_in);
      }
    else {
      if ( ! importer.is_null ())
        X_colMap->endImport (X_in, *importer, INSERT, restrictedMode);
      localResidual (A, *X_colMap, *B_rowMap, *R_rowMap, offsets);
    }

    //
    // If R_in does not have constant stride,
    // then we can't let the kernel write directly to R_in.
    // Instead, we have to use the cached row (== range)
    // Map MV as temporary storage.
    //
    if (! R_in.isConstantStride () ) {
      // We need to be sure to do a copy out in this case.
      Tpetra::deep_copy (R_in, *R_rowMap);
    }
  }

  // If the range Map is a locally replicated Map, sum up
  // contributions from each process. 
  if (R_is_replicated) {
    ProfilingRegion regionReduce ("Tpetra::CrsMatrix::residual: Reduce Y");
    R_in.reduce ();
  }
}





  } // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_RESIDUAL_HPP
