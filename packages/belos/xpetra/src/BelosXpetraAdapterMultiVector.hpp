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

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Exceptions.hpp>
#include <Xpetra_MultiVector.hpp>

#include <Xpetra_TpetraMultiVector.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosTypes.hpp>
#include <BelosMultiVecTraits.hpp>
#include <BelosOperatorTraits.hpp>

#include <BelosTpetraAdapter.hpp>
#include <TpetraCore_config.h>

#ifdef HAVE_BELOS_TSQR
namespace BelosXpetraTsqrImpl {

  template<class Scalar, class LO, class GO, class Node>
  class XpetraStubTsqrAdaptor : public Teuchos::ParameterListAcceptorDefaultBase {
  public:
    typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MV;
    typedef Scalar scalar_type;
    typedef LO ordinal_type;
    typedef Teuchos::SerialDenseMatrix<ordinal_type, scalar_type> dense_matrix_type;
    typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

    XpetraStubTsqrAdaptor (const Teuchos::RCP<Teuchos::ParameterList>& /* plist */)
    {}

    XpetraStubTsqrAdaptor () {}

    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters () const
    {
      return Teuchos::rcp (new Teuchos::ParameterList);
    }

    void
    setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& /* plist */)
    {}

    void
    factorExplicit (MV& /* A */,
                    MV& /* Q */,
                    dense_matrix_type& /* R */,
                    const bool /* forceNonnegativeDiagonal */ = false)
    {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "Xpetra TSQR adaptor is only implemented "
         "for the Tpetra case");
    }

    int
    revealRank (MV& /* Q */,
                dense_matrix_type& /* R */,
                const magnitude_type& /* tol */)
    {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "Xpetra TSQR adaptor is only implemented "
         "for the Tpetra case");
    }
  };

  template<class Scalar, class LO, class GO, class Node>
  class XpetraTpetraTsqrAdaptor : public Teuchos::ParameterListAcceptorDefaultBase {
  public:
    typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MV;
    typedef Scalar scalar_type;
    typedef LO ordinal_type;
    typedef Teuchos::SerialDenseMatrix<ordinal_type, scalar_type> dense_matrix_type;
    typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

    XpetraTpetraTsqrAdaptor (const Teuchos::RCP<Teuchos::ParameterList>& plist)
      : tpetraImpl_ (plist)
    {}

    XpetraTpetraTsqrAdaptor () {}

    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters () const
    {
      return tpetraImpl_.getValidParameters ();
    }

    void
    setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist)
    {
      tpetraImpl_.setParameterList (plist);
    }

    void
    factorExplicit (MV& A,
                    MV& Q,
                    dense_matrix_type& R,
                    const bool forceNonnegativeDiagonal = false)
    {
      if (A.getMap()->lib() == Xpetra::UseTpetra) {
        tpetraImpl_.factorExplicit (toTpetra (A), toTpetra (Q), R, forceNonnegativeDiagonal);
        return;
      }
      XPETRA_FACTORY_END;
    }

    int
    revealRank (MV& Q,
                dense_matrix_type& R,
                const magnitude_type& tol)
    {
      if (Q.getMap()->lib() == Xpetra::UseTpetra) {
        return tpetraImpl_.revealRank (toTpetra (Q), R, tol);
      }
      XPETRA_FACTORY_END;
    }

  private:
    typedef ::Tpetra::TsqrAdaptor< ::Tpetra::MultiVector<Scalar,LO,GO,Node> > tpetra_tsqr_adaptor_type;
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

  /*! \brief Template specialization of Belos::MultiVecTraits class using the Xpetra::MultiVector class.
    This interface will ensure that any Xpetra::MultiVector will be accepted by the Belos
    templated solvers.
  */
  template<class Scalar, class LO, class GO, class Node>
  class MultiVecTraits<Scalar, Xpetra::MultiVector<Scalar,LO,GO,Node> > {
  private:
    typedef Xpetra::TpetraMultiVector<Scalar,LO,GO,Node>                   TpetraMultiVector;
    typedef MultiVecTraits<Scalar,Tpetra::MultiVector<Scalar,LO,GO,Node> > MultiVecTraitsTpetra;

  public:

#ifdef HAVE_BELOS_XPETRA_TIMERS
    static RCP<Teuchos::Time> mvTimesMatAddMvTimer_, mvTransMvTimer_;
#endif

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > Clone( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const int numvecs )
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::Clone(toTpetra(mv), numvecs)));

      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneCopy( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv))));

      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneCopy( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv), index)));

      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneCopy (const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const Teuchos::Range1D& index)
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv), index)));

      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneViewNonConst( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneViewNonConst(toTpetra(mv), index)));

      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneViewNonConst(Xpetra::MultiVector<Scalar,LO,GO,Node>& mv,
                      const Teuchos::Range1D& index)
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneViewNonConst(toTpetra(mv), index)));

      XPETRA_FACTORY_END;
    }

    static RCP<const Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneView(const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        //TODO: double check if the const_cast is safe here.
        RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> > r = MultiVecTraitsTpetra::CloneView(toTpetra(mv), index);
        return rcp(new TpetraMultiVector(Teuchos::rcp_const_cast<Tpetra::MultiVector<Scalar,LO,GO,Node> >(r)));
      }

      XPETRA_FACTORY_END;
    }

    static RCP<const Xpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneView (const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const Teuchos::Range1D& index)
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        //TODO: double check if the const_cast is safe here.
        RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> > r = MultiVecTraitsTpetra::CloneView(toTpetra(mv), index);
        return rcp(new TpetraMultiVector(Teuchos::rcp_const_cast<Tpetra::MultiVector<Scalar,LO,GO,Node> >(r)));
      }

      XPETRA_FACTORY_END;
    }

    static ptrdiff_t GetGlobalLength( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return MultiVecTraitsTpetra::GetGlobalLength(toTpetra(mv));

      XPETRA_FACTORY_END;
    }

    static int GetNumberVecs( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return MultiVecTraitsTpetra::GetNumberVecs(toTpetra(mv));

      XPETRA_FACTORY_END;
    }

    static bool HasConstantStride( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return  MultiVecTraitsTpetra::HasConstantStride(toTpetra(mv));

      XPETRA_FACTORY_END;
    }

    static void MvTimesMatAddMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A,
                                 const Teuchos::SerialDenseMatrix<int,Scalar>& B,
                                 Scalar beta, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_BELOS_XPETRA_TIMERS
      Teuchos::TimeMonitor lcltimer(*mvTimesMatAddMvTimer_);
#endif


      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvTimesMatAddMv(alpha, toTpetra(A), B, beta, toTpetra(mv));
        return;
      }

      XPETRA_FACTORY_END;
    }

    static void MvAddMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, Scalar beta, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvAddMv(alpha, toTpetra(A), beta, toTpetra(B), toTpetra(mv));
        return;
      }

      XPETRA_FACTORY_END;
    }

    static void MvScale ( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, Scalar alpha )
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvScale(toTpetra(mv), alpha);
        return;
      }

      XPETRA_FACTORY_END;
    }

    static void MvScale ( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<Scalar>& alphas )
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvScale(toTpetra(mv), alphas);
        return;
      }

      XPETRA_FACTORY_END;
    }

    static void MvTransMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, Teuchos::SerialDenseMatrix<int,Scalar>& C)
    {

#ifdef HAVE_BELOS_XPETRA_TIMERS
      Teuchos::TimeMonitor lcltimer(*mvTransMvTimer_);
#endif

      if (A.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvTransMv(alpha, toTpetra(A), toTpetra(B), C);
        return;
      }

      XPETRA_FACTORY_END;
    }

    static void MvDot( const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, std::vector<Scalar> &dots)
    {

      if (A.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvDot(toTpetra(A), toTpetra(B), dots);
        return;
      }

      XPETRA_FACTORY_END;
    }

    static void MvNorm(const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec, NormType type=TwoNorm)
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvNorm(toTpetra(mv), normvec, type);
        return;
      }

      XPETRA_FACTORY_END;
    }

    static void SetBlock( const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const std::vector<int>& index, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::SetBlock(toTpetra(A), index, toTpetra(mv));
        return;
      }

      XPETRA_FACTORY_END;
    }

    static void
    SetBlock (const Xpetra::MultiVector<Scalar,LO,GO,Node>& A,
              const Teuchos::Range1D& index,
              Xpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::SetBlock(toTpetra(A), index, toTpetra(mv));
        return;
      }

      XPETRA_FACTORY_END;
    }

    static void
    Assign (const Xpetra::MultiVector<Scalar,LO,GO,Node>& A,
            Xpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::Assign(toTpetra(A), toTpetra(mv));
        return;
      }

      XPETRA_FACTORY_END;
    }

    static void MvRandom( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvRandom(toTpetra(mv));
        return;
      }

      XPETRA_FACTORY_END;
    }

    static void MvInit( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() )
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvInit(toTpetra(mv), alpha);
        return;
      }

      XPETRA_FACTORY_END;
    }

    static void MvPrint( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::ostream& os )
    {

      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvPrint(toTpetra(mv), os);
        return;
      }

      XPETRA_FACTORY_END;
    }

#ifdef HAVE_BELOS_TSQR
    typedef BelosXpetraTsqrImpl::XpetraTpetraTsqrAdaptor<Scalar, LO, GO, Node> tsqr_adaptor_type;
#endif // HAVE_BELOS_TSQR
  };

} // end of Belos namespace

#endif
