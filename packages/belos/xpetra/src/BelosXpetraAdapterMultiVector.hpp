// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_XPETRA_ADAPTER_MULTIVECTOR_HPP
#define BELOS_XPETRA_ADAPTER_MULTIVECTOR_HPP

#include "Xpetra_BlockedMultiVector_decl.hpp"
#include "Xpetra_MultiVectorFactory_decl.hpp"
#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Exceptions.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>

#include <Xpetra_TpetraMultiVector.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosMultiVecTraits.hpp>
#include <BelosOperatorTraits.hpp>
#include <BelosTypes.hpp>

#include <BelosTpetraAdapter.hpp>
#include <TpetraCore_config.h>

#ifdef HAVE_BELOS_TSQR
namespace BelosXpetraTsqrImpl {

template <class Scalar, class LO, class GO, class Node>
class XpetraStubTsqrAdaptor : public Teuchos::ParameterListAcceptorDefaultBase {
public:
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MV;
  typedef Scalar scalar_type;
  typedef LO ordinal_type;
  typedef Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>
      dense_matrix_type;
  typedef
      typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  XpetraStubTsqrAdaptor(
      const Teuchos::RCP<Teuchos::ParameterList> & /* plist */) {}

  XpetraStubTsqrAdaptor() {}

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const {
    return Teuchos::rcp(new Teuchos::ParameterList);
  }

  void
  setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & /* plist */) {}

  void factorExplicit(MV & /* A */, MV & /* Q */, dense_matrix_type & /* R */,
                      const bool /* forceNonnegativeDiagonal */ = false) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "Xpetra TSQR adaptor is only implemented "
                               "for the Tpetra case");
  }

  int revealRank(MV & /* Q */, dense_matrix_type & /* R */,
                 const magnitude_type & /* tol */) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "Xpetra TSQR adaptor is only implemented "
                               "for the Tpetra case");
  }
};

template <class Scalar, class LO, class GO, class Node>
class XpetraTpetraTsqrAdaptor
    : public Teuchos::ParameterListAcceptorDefaultBase {
public:
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MV;
  typedef Scalar scalar_type;
  typedef LO ordinal_type;
  typedef Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>
      dense_matrix_type;
  typedef
      typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  XpetraTpetraTsqrAdaptor(const Teuchos::RCP<Teuchos::ParameterList> &plist)
      : tpetraImpl_(plist) {}

  XpetraTpetraTsqrAdaptor() {}

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const {
    return tpetraImpl_.getValidParameters();
  }

  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> &plist) {
    tpetraImpl_.setParameterList(plist);
  }

  void factorExplicit(MV &A, MV &Q, dense_matrix_type &R,
                      const bool forceNonnegativeDiagonal = false) {
    if (A.getMap()->lib() == Xpetra::UseTpetra) {
      tpetraImpl_.factorExplicit(toTpetra(A), toTpetra(Q), R,
                                 forceNonnegativeDiagonal);
      return;
    }
    XPETRA_FACTORY_END;
  }

  int revealRank(MV &Q, dense_matrix_type &R, const magnitude_type &tol) {
    if (Q.getMap()->lib() == Xpetra::UseTpetra) {
      return tpetraImpl_.revealRank(toTpetra(Q), R, tol);
    }
    XPETRA_FACTORY_END;
  }

private:
  typedef ::Tpetra::TsqrAdaptor<::Tpetra::MultiVector<Scalar, LO, GO, Node>>
      tpetra_tsqr_adaptor_type;
  tpetra_tsqr_adaptor_type tpetraImpl_;
};

} // namespace BelosXpetraTsqrImpl
#endif // HAVE_BELOS_TSQR

namespace Belos {

using Teuchos::RCP;
using Teuchos::rcp;

////////////////////////////////////////////////////////////////////
//
// Implementation of the Belos::MultiVecTraits for Xpetra::MultiVector.
//
////////////////////////////////////////////////////////////////////

/*! \brief Template specialization of Belos::MultiVecTraits class using the
  Xpetra::MultiVector class. This interface will ensure that any
  Xpetra::MultiVector will be accepted by the Belos templated solvers.
*/
template <class Scalar, class LO, class GO, class Node>
class MultiVecTraits<Scalar, Xpetra::MultiVector<Scalar, LO, GO, Node>> {
private:
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MV;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Xpetra::TpetraMultiVector<Scalar, LO, GO, Node> TpetraMultiVector;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> Vector;
  typedef MultiVecTraits<Scalar, Tpetra::MultiVector<Scalar, LO, GO, Node>>
      MultiVecTraitsTpetra;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

  static RCP<const BlockedMultiVector> asBlocked(const MV &mv) {
    return Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(
        Teuchos::rcpFromRef(mv), false);
  }

  static RCP<BlockedMultiVector> asBlocked(MV &mv) {
    return Teuchos::rcp_dynamic_cast<BlockedMultiVector>(
        Teuchos::rcpFromRef(mv), false);
  }

  static size_t getNumBlocks(const BlockedMultiVector &bmv) {
    return bmv.getBlockedMap()->getNumMaps();
  }

public:
#ifdef HAVE_BELOS_XPETRA_TIMERS
  static RCP<Teuchos::Time> mvTimesMatAddMvTimer_, mvTransMvTimer_;
#endif

  static RCP<MV> Clone(const MV &mv, const int numvecs) {
    if (auto bmv = asBlocked(mv))
      return rcp(new BlockedMultiVector(bmv->getBlockedMap(), numvecs, true));

    return Xpetra::MultiVectorFactory<Scalar, LO, GO, Node>::Build(mv.getMap(),
                                                                   numvecs);
  }

  static RCP<MV> CloneCopy(const MV &mv) {
    if (auto bmv = asBlocked(mv)) {
      RCP<BlockedMultiVector> out = rcp(new BlockedMultiVector(
          bmv->getBlockedMap(), mv.getNumVectors(), true));
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> out_r = out->getMultiVector(r);
        RCP<MV> in_r = bmv->getMultiVector(r);
        MultiVecTraits<Scalar, MV>::Assign(*in_r, *out_r);
      }
      return out;
    }

    auto copy = Xpetra::MultiVectorFactory<Scalar, LO, GO, Node>::Build(
        mv.getMap(), mv.getNumVectors());
    *copy = mv;
    return copy;
  }

  static RCP<MV> CloneCopy(const MV &mv, const std::vector<int> &index) {
    if (auto bmv = asBlocked(mv)) {
      RCP<BlockedMultiVector> out =
          rcp(new BlockedMultiVector(bmv->getBlockedMap(), index.size(), true));
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> out_r = out->getMultiVector(r);
        RCP<MV> in_r = bmv->getMultiVector(r);
        RCP<MV> tmp = MultiVecTraits<Scalar, MV>::CloneCopy(*in_r, index);
        MultiVecTraits<Scalar, MV>::Assign(*tmp, *out_r);
      }
      return out;
    }

#ifdef HAVE_TPETRA_DEBUG
    const char fnName[] = "Belos::MultiVecTraits::CloneCopy(mv,index)";
    const size_t inNumVecs = mv.getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(
        index.size() > 0 && *std::min_element(index.begin(), index.end()) < 0,
        std::runtime_error, fnName << ": All indices must be nonnegative.");
    TEUCHOS_TEST_FOR_EXCEPTION(
        index.size() > 0 && static_cast<size_t>(*std::max_element(
                                index.begin(), index.end())) >= inNumVecs,
        std::runtime_error,
        fnName << ": All indices must be strictly less than the number of "
                  "columns "
               << inNumVecs << " of the input multivector mv.");
#endif // HAVE_TPETRA_DEBUG

    RCP<MV> X_copy = Xpetra::MultiVectorFactory<Scalar, LO, GO, Node>::Build(
        mv.getMap(), index.size());

    for (size_t j = 0; j < index.size(); ++j) {
      RCP<const Vector> src = mv.getVector(static_cast<size_t>(index[j]));
      RCP<Vector> dst = X_copy->getVectorNonConst(j);
      *dst = *src;
    }

    return X_copy;
  }

  static RCP<MV> CloneCopy(const MV &mv, const Teuchos::Range1D &index) {
    if (auto bmv = asBlocked(mv)) {
      RCP<BlockedMultiVector> out =
          rcp(new BlockedMultiVector(bmv->getBlockedMap(), index.size(), true));
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> out_r = out->getMultiVector(r);
        RCP<MV> in_r = bmv->getMultiVector(r);
        RCP<MV> tmp = MultiVecTraits<Scalar, MV>::CloneCopy(*in_r, index);
        MultiVecTraits<Scalar, MV>::Assign(*tmp, *out_r);
      }
      return out;
    }

#ifdef HAVE_TPETRA_DEBUG
    const char fnName[] = "Belos::MultiVecTraits::CloneCopy(mv,index)";
    const size_t inNumVecs = mv.getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::runtime_error,
                               fnName << ": Lower bound must be nonnegative.");
    TEUCHOS_TEST_FOR_EXCEPTION(
        index.size() > 0 && static_cast<size_t>(index.ubound()) >= inNumVecs,
        std::runtime_error,
        fnName << ": Upper bound must be strictly less than the number of "
                  "columns "
               << inNumVecs << " of the input multivector mv.");
#endif // HAVE_TPETRA_DEBUG

    RCP<MV> X_copy = Xpetra::MultiVectorFactory<Scalar, LO, GO, Node>::Build(
        mv.getMap(), index.size());

    for (int j = index.lbound(), k = 0; j <= index.ubound(); ++j, ++k) {
      RCP<const Vector> src = mv.getVector(static_cast<size_t>(j));
      RCP<Vector> dst = X_copy->getVectorNonConst(static_cast<size_t>(k));
      *dst = *src;
    }

    return X_copy;
  }

  static RCP<MV> CloneViewNonConst(MV &mv, const std::vector<int> &index) {
    if (auto bmv = asBlocked(mv)) {
      std::vector<RCP<MV>> blocks;
      blocks.reserve(getNumBlocks(*bmv));
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> in_r = bmv->getMultiVector(r);
        blocks.push_back(
            MultiVecTraits<Scalar, MV>::CloneViewNonConst(*in_r, index));
      }
      return rcp(new BlockedMultiVector(bmv->getBlockedMap(), blocks));
    }

    if (mv.getMap()->lib() == Xpetra::UseTpetra)
      return rcp(new TpetraMultiVector(
          MultiVecTraitsTpetra::CloneViewNonConst(toTpetra(mv), index)));

    XPETRA_FACTORY_END;
  }

  static RCP<MV> CloneViewNonConst(MV &mv, const Teuchos::Range1D &index) {
    if (auto bmv = asBlocked(mv)) {
      std::vector<RCP<MV>> blocks;
      blocks.reserve(getNumBlocks(*bmv));
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> in_r = bmv->getMultiVector(r);
        blocks.push_back(
            MultiVecTraits<Scalar, MV>::CloneViewNonConst(*in_r, index));
      }
      return rcp(new BlockedMultiVector(bmv->getBlockedMap(), blocks));
    }

    if (mv.getMap()->lib() == Xpetra::UseTpetra)
      return rcp(new TpetraMultiVector(
          MultiVecTraitsTpetra::CloneViewNonConst(toTpetra(mv), index)));

    XPETRA_FACTORY_END;
  }

  static RCP<const MV> CloneView(const MV &mv, const std::vector<int> &index) {
    if (auto bmv = asBlocked(mv)) {
      std::vector<RCP<MV>> blocks;
      blocks.reserve(getNumBlocks(*bmv));
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> in_r = bmv->getMultiVector(r);
        RCP<const MV> view_r =
            MultiVecTraits<Scalar, MV>::CloneView(*in_r, index);
        blocks.push_back(Teuchos::rcp_const_cast<MV>(view_r));
      }
      return rcp(new BlockedMultiVector(bmv->getBlockedMap(), blocks));
    }

    if (mv.getMap()->lib() == Xpetra::UseTpetra) {
      // TODO: double check if the const_cast is safe here.
      RCP<const Tpetra::MultiVector<Scalar, LO, GO, Node>> r =
          MultiVecTraitsTpetra::CloneView(toTpetra(mv), index);
      return rcp(new TpetraMultiVector(
          Teuchos::rcp_const_cast<Tpetra::MultiVector<Scalar, LO, GO, Node>>(
              r)));
    }

    XPETRA_FACTORY_END;
  }

  static RCP<const MV> CloneView(const MV &mv, const Teuchos::Range1D &index) {
    if (auto bmv = asBlocked(mv)) {
      std::vector<RCP<MV>> blocks;
      blocks.reserve(getNumBlocks(*bmv));
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> in_r = bmv->getMultiVector(r);
        RCP<const MV> view_r =
            MultiVecTraits<Scalar, MV>::CloneView(*in_r, index);
        blocks.push_back(Teuchos::rcp_const_cast<MV>(view_r));
      }
      return rcp(new BlockedMultiVector(bmv->getBlockedMap(), blocks));
    }

    if (mv.getMap()->lib() == Xpetra::UseTpetra) {
      // TODO: double check if the const_cast is safe here.
      RCP<const Tpetra::MultiVector<Scalar, LO, GO, Node>> r =
          MultiVecTraitsTpetra::CloneView(toTpetra(mv), index);
      return rcp(new TpetraMultiVector(
          Teuchos::rcp_const_cast<Tpetra::MultiVector<Scalar, LO, GO, Node>>(
              r)));
    }

    XPETRA_FACTORY_END;
  }

  static ptrdiff_t GetGlobalLength(const MV &mv) {
    return static_cast<ptrdiff_t>(mv.getGlobalLength());
  }

  static int GetNumberVecs(const MV &mv) {
    return static_cast<int>(mv.getNumVectors());
  }

  static bool HasConstantStride(const MV &mv) {
    if (auto bmv = asBlocked(mv)) {
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> in_r = bmv->getMultiVector(r);
        if (!MultiVecTraits<Scalar, MV>::HasConstantStride(*in_r))
          return false;
      }
      return true;
    }

    if (mv.getMap()->lib() == Xpetra::UseTpetra)
      return MultiVecTraitsTpetra::HasConstantStride(toTpetra(mv));

    XPETRA_FACTORY_END;
  }

  static void MvTimesMatAddMv(Scalar alpha, const MV &A,
                              const Teuchos::SerialDenseMatrix<int, Scalar> &B,
                              Scalar beta, MV &mv) {
#ifdef HAVE_BELOS_XPETRA_TIMERS
    Teuchos::TimeMonitor lcltimer(*mvTimesMatAddMvTimer_);
#endif

    if (auto bmv = asBlocked(mv)) {
      auto A_bmv = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(
          Teuchos::rcpFromRef(A), true);
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> mv_r = bmv->getMultiVector(r);
        RCP<MV> A_r = A_bmv->getMultiVector(r);
        MultiVecTraits<Scalar, MV>::MvTimesMatAddMv(alpha, *A_r, B, beta,
                                                    *mv_r);
      }
      return;
    }

    if (mv.getMap()->lib() == Xpetra::UseTpetra) {
      MultiVecTraitsTpetra::MvTimesMatAddMv(alpha, toTpetra(A), B, beta,
                                            toTpetra(mv));
      return;
    }

    XPETRA_FACTORY_END;
  }

  static void MvAddMv(Scalar alpha, const MV &A, Scalar beta, const MV &B,
                      MV &mv) {
    if (auto bmv = asBlocked(mv)) {
      auto A_bmv = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(
          Teuchos::rcpFromRef(A), true);
      auto B_bmv = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(
          Teuchos::rcpFromRef(B), true);
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> mv_r = bmv->getMultiVector(r);
        RCP<MV> A_r = A_bmv->getMultiVector(r);
        RCP<MV> B_r = B_bmv->getMultiVector(r);
        MultiVecTraits<Scalar, MV>::MvAddMv(alpha, *A_r, beta, *B_r, *mv_r);
      }
      return;
    }

    mv.update(alpha, A, beta, B, Teuchos::ScalarTraits<Scalar>::zero());
  }

  static void MvScale(MV &mv, Scalar alpha) {
    if (auto bmv = asBlocked(mv)) {
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> mv_r = bmv->getMultiVector(r);
        MultiVecTraits<Scalar, MV>::MvScale(*mv_r, alpha);
      }
      return;
    }

    mv.scale(alpha);
  }

  static void MvScale(MV &mv, const std::vector<Scalar> &alphas) {
    if (auto bmv = asBlocked(mv)) {
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> mv_r = bmv->getMultiVector(r);
        MultiVecTraits<Scalar, MV>::MvScale(*mv_r, alphas);
      }
      return;
    }

    mv.scale(Teuchos::ArrayView<const Scalar>(alphas.data(), alphas.size()));
  }

  static void MvTransMv(Scalar alpha, const MV &A, const MV &B,
                        Teuchos::SerialDenseMatrix<int, Scalar> &C) {
#ifdef HAVE_BELOS_XPETRA_TIMERS
    Teuchos::TimeMonitor lcltimer(*mvTransMvTimer_);
#endif

    if (auto A_bmv = asBlocked(A)) {
      auto B_bmv = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(
          Teuchos::rcpFromRef(B), true);
      C.putScalar(Teuchos::ScalarTraits<Scalar>::zero());
      Teuchos::SerialDenseMatrix<int, Scalar> Ctmp(C.numRows(), C.numCols());

      for (size_t r = 0; r < getNumBlocks(*A_bmv); ++r) {
        RCP<MV> A_r = A_bmv->getMultiVector(r);
        RCP<MV> B_r = B_bmv->getMultiVector(r);
        Ctmp.putScalar(Teuchos::ScalarTraits<Scalar>::zero());
        MultiVecTraits<Scalar, MV>::MvTransMv(alpha, *A_r, *B_r, Ctmp);
        for (int i = 0; i < C.numRows(); ++i)
          for (int j = 0; j < C.numCols(); ++j)
            C(i, j) += Ctmp(i, j);
      }
      return;
    }

    if (A.getMap()->lib() == Xpetra::UseTpetra) {
      MultiVecTraitsTpetra::MvTransMv(alpha, toTpetra(A), toTpetra(B), C);
      return;
    }

    XPETRA_FACTORY_END;
  }

  static void MvDot(const MV &A, const MV &B, std::vector<Scalar> &dots) {
    if (auto A_bmv = asBlocked(A)) {
      auto B_bmv = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(
          Teuchos::rcpFromRef(B), true);
      std::fill(dots.begin(), dots.end(),
                Teuchos::ScalarTraits<Scalar>::zero());
      std::vector<Scalar> tmp(dots.size());

      for (size_t r = 0; r < getNumBlocks(*A_bmv); ++r) {
        RCP<MV> A_r = A_bmv->getMultiVector(r);
        RCP<MV> B_r = B_bmv->getMultiVector(r);
        std::fill(tmp.begin(), tmp.end(),
                  Teuchos::ScalarTraits<Scalar>::zero());
        MultiVecTraits<Scalar, MV>::MvDot(*A_r, *B_r, tmp);
        for (size_t j = 0; j < dots.size(); ++j)
          dots[j] += tmp[j];
      }
      return;
    }

    mvDotImpl(A, B, dots);
  }

  static void
  MvNorm(const MV &mv,
         std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>
             &normvec,
         NormType type = TwoNorm) {
    if (auto bmv = asBlocked(mv)) {
      const magnitude_type zero = Teuchos::ScalarTraits<magnitude_type>::zero();
      std::fill(normvec.begin(), normvec.end(), zero);
      std::vector<magnitude_type> tmp(normvec.size(), zero);

      switch (type) {
      case OneNorm:
        for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
          RCP<MV> mv_r = bmv->getMultiVector(r);
          std::fill(tmp.begin(), tmp.end(), zero);
          MultiVecTraits<Scalar, MV>::MvNorm(*mv_r, tmp, OneNorm);
          for (size_t j = 0; j < normvec.size(); ++j)
            normvec[j] += tmp[j];
        }
        return;

      case TwoNorm:
        for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
          RCP<MV> mv_r = bmv->getMultiVector(r);
          std::fill(tmp.begin(), tmp.end(), zero);
          MultiVecTraits<Scalar, MV>::MvNorm(*mv_r, tmp, TwoNorm);
          for (size_t j = 0; j < normvec.size(); ++j)
            normvec[j] += tmp[j] * tmp[j];
        }
        for (size_t j = 0; j < normvec.size(); ++j)
          normvec[j] =
              Teuchos::ScalarTraits<magnitude_type>::squareroot(normvec[j]);
        return;

      case InfNorm:
        for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
          RCP<MV> mv_r = bmv->getMultiVector(r);
          std::fill(tmp.begin(), tmp.end(), zero);
          MultiVecTraits<Scalar, MV>::MvNorm(*mv_r, tmp, InfNorm);
          for (size_t j = 0; j < normvec.size(); ++j)
            if (tmp[j] > normvec[j])
              normvec[j] = tmp[j];
        }
        return;

      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
                                   "Belos::MultiVecTraits<Xpetra::MultiVector>:"
                                   ":MvNorm: invalid NormType.");
      }
    }

    Teuchos::ArrayView<magnitude_type> norms(normvec.data(), normvec.size());
    switch (type) {
    case OneNorm:
      mv.norm1(norms);
      break;
    case TwoNorm:
      mv.norm2(norms);
      break;
    case InfNorm:
      mv.normInf(norms);
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
                                 "Belos::MultiVecTraits<Xpetra::MultiVector>::"
                                 "MvNorm: invalid NormType.");
    }
  }

  static void SetBlock(const MV &A, const std::vector<int> &index, MV &mv) {
    if (auto bmv = asBlocked(mv)) {
      auto A_bmv = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(
          Teuchos::rcpFromRef(A), true);
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> mv_r = bmv->getMultiVector(r);
        RCP<MV> A_r = A_bmv->getMultiVector(r);
        MultiVecTraits<Scalar, MV>::SetBlock(*A_r, index, *mv_r);
      }
      return;
    }

    if (mv.getMap()->lib() == Xpetra::UseTpetra) {
      MultiVecTraitsTpetra::SetBlock(toTpetra(A), index, toTpetra(mv));
      return;
    }

    XPETRA_FACTORY_END;
  }

  static void SetBlock(const MV &A, const Teuchos::Range1D &index, MV &mv) {
    if (auto bmv = asBlocked(mv)) {
      auto A_bmv = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(
          Teuchos::rcpFromRef(A), true);
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> mv_r = bmv->getMultiVector(r);
        RCP<MV> A_r = A_bmv->getMultiVector(r);
        MultiVecTraits<Scalar, MV>::SetBlock(*A_r, index, *mv_r);
      }
      return;
    }

    if (mv.getMap()->lib() == Xpetra::UseTpetra) {
      MultiVecTraitsTpetra::SetBlock(toTpetra(A), index, toTpetra(mv));
      return;
    }

    XPETRA_FACTORY_END;
  }

  static void Assign(const MV &A, MV &mv) {
    if (auto bmv = asBlocked(mv)) {
      auto A_bmv = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(
          Teuchos::rcpFromRef(A), true);
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> mv_r = bmv->getMultiVector(r);
        RCP<MV> A_r = A_bmv->getMultiVector(r);
        MultiVecTraits<Scalar, MV>::Assign(*A_r, *mv_r);
      }
      return;
    }

    mv = A;
  }

  static void MvRandom(MV &mv) {
    if (auto bmv = asBlocked(mv)) {
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> mv_r = bmv->getMultiVector(r);
        MultiVecTraits<Scalar, MV>::MvRandom(*mv_r);
      }
      return;
    }

    mv.randomize();
  }

  static void MvInit(MV &mv,
                     Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero()) {
    if (auto bmv = asBlocked(mv)) {
      for (size_t r = 0; r < getNumBlocks(*bmv); ++r) {
        RCP<MV> mv_r = bmv->getMultiVector(r);
        MultiVecTraits<Scalar, MV>::MvInit(*mv_r, alpha);
      }
      return;
    }

    mv.putScalar(alpha);
  }

  static void MvPrint(const MV &mv, std::ostream &os) {
    Teuchos::FancyOStream fos(Teuchos::rcpFromRef(os));
    mv.describe(fos, Teuchos::Describable::verbLevel_default);
  }

private:
  static void mvDotImpl(const MV &A, const MV &B, std::vector<Scalar> &dots) {
    A.dot(B, Teuchos::ArrayView<Scalar>(dots.data(), dots.size()));
  }
};

} // namespace Belos

#endif
