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

#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraMultiVector.hpp>
#endif

#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraMultiVector.hpp>
#endif

#include <BelosConfigDefs.hpp>
#include <BelosTypes.hpp>
#include <BelosMultiVecTraits.hpp>
#include <BelosOperatorTraits.hpp>

#ifdef HAVE_XPETRA_EPETRA
#include <BelosEpetraAdapter.hpp>
#endif

#ifdef HAVE_XPETRA_TPETRA
#include <BelosTpetraAdapter.hpp>
#include <TpetraCore_config.h>
#endif

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

#ifdef HAVE_XPETRA_TPETRA
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
      XPETRA_FACTORY_ERROR_IF_EPETRA(A.getMap()->lib());
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
      XPETRA_FACTORY_ERROR_IF_EPETRA(Q.getMap()->lib());
      XPETRA_FACTORY_END;
    }

  private:
    typedef ::Tpetra::TsqrAdaptor< ::Tpetra::MultiVector<Scalar,LO,GO,Node> > tpetra_tsqr_adaptor_type;
    tpetra_tsqr_adaptor_type tpetraImpl_;
  };
#endif // HAVE_XPETRA_TPETRA

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
#ifdef HAVE_XPETRA_TPETRA
    typedef Xpetra::TpetraMultiVector<Scalar,LO,GO,Node>                   TpetraMultiVector;
    typedef MultiVecTraits<Scalar,Tpetra::MultiVector<Scalar,LO,GO,Node> > MultiVecTraitsTpetra;
#endif

  public:

#ifdef HAVE_BELOS_XPETRA_TIMERS
    static RCP<Teuchos::Time> mvTimesMatAddMvTimer_, mvTransMvTimer_;
#endif

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > Clone( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const int numvecs )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::Clone(toTpetra(mv), numvecs)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneCopy( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv))));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneCopy( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv), index)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneCopy (const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const Teuchos::Range1D& index)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv), index)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneViewNonConst( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneViewNonConst(toTpetra(mv), index)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneViewNonConst(Xpetra::MultiVector<Scalar,LO,GO,Node>& mv,
                      const Teuchos::Range1D& index)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneViewNonConst(toTpetra(mv), index)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static RCP<const Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneView(const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        //TODO: double check if the const_cast is safe here.
        RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> > r = MultiVecTraitsTpetra::CloneView(toTpetra(mv), index);
        return rcp(new TpetraMultiVector(Teuchos::rcp_const_cast<Tpetra::MultiVector<Scalar,LO,GO,Node> >(r)));
      }
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static RCP<const Xpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneView (const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const Teuchos::Range1D& index)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        //TODO: double check if the const_cast is safe here.
        RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> > r = MultiVecTraitsTpetra::CloneView(toTpetra(mv), index);
        return rcp(new TpetraMultiVector(Teuchos::rcp_const_cast<Tpetra::MultiVector<Scalar,LO,GO,Node> >(r)));
      }
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static ptrdiff_t GetGlobalLength( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return MultiVecTraitsTpetra::GetGlobalLength(toTpetra(mv));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static int GetNumberVecs( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return MultiVecTraitsTpetra::GetNumberVecs(toTpetra(mv));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static bool HasConstantStride( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return  MultiVecTraitsTpetra::HasConstantStride(toTpetra(mv));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvTimesMatAddMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A,
                                 const Teuchos::SerialDenseMatrix<int,Scalar>& B,
                                 Scalar beta, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_BELOS_XPETRA_TIMERS
      Teuchos::TimeMonitor lcltimer(*mvTimesMatAddMvTimer_);
#endif


#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvTimesMatAddMv(alpha, toTpetra(A), B, beta, toTpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvAddMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, Scalar beta, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvAddMv(alpha, toTpetra(A), beta, toTpetra(B), toTpetra(mv));
        return;
      }

#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvScale ( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, Scalar alpha )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvScale(toTpetra(mv), alpha);
        return;
      }
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvScale ( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<Scalar>& alphas )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvScale(toTpetra(mv), alphas);
        return;
      }
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvTransMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, Teuchos::SerialDenseMatrix<int,Scalar>& C)
    {

#ifdef HAVE_BELOS_XPETRA_TIMERS
      Teuchos::TimeMonitor lcltimer(*mvTransMvTimer_);
#endif

#ifdef HAVE_XPETRA_TPETRA
      if (A.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvTransMv(alpha, toTpetra(A), toTpetra(B), C);
        return;
      }
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(A.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvDot( const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, std::vector<Scalar> &dots)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (A.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvDot(toTpetra(A), toTpetra(B), dots);
        return;
      }
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(A.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvNorm(const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec, NormType type=TwoNorm)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvNorm(toTpetra(mv), normvec, type);
        return;
      }
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void SetBlock( const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const std::vector<int>& index, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::SetBlock(toTpetra(A), index, toTpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void
    SetBlock (const Xpetra::MultiVector<Scalar,LO,GO,Node>& A,
              const Teuchos::Range1D& index,
              Xpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::SetBlock(toTpetra(A), index, toTpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void
    Assign (const Xpetra::MultiVector<Scalar,LO,GO,Node>& A,
            Xpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::Assign(toTpetra(A), toTpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvRandom( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvRandom(toTpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvInit( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvInit(toTpetra(mv), alpha);
        return;
      }
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvPrint( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::ostream& os )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
        MultiVecTraitsTpetra::MvPrint(toTpetra(mv), os);
        return;
      }
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

#ifdef HAVE_BELOS_TSQR
#  ifdef HAVE_XPETRA_TPETRA
    typedef BelosXpetraTsqrImpl::XpetraTpetraTsqrAdaptor<Scalar, LO, GO, Node> tsqr_adaptor_type;
#  else
    typedef BelosXpetraTsqrImpl::XpetraStubTsqrAdaptor<Scalar, LO, GO, Node> tsqr_adaptor_type;
#  endif // HAVE_XPETRA_TPETRA
#endif // HAVE_BELOS_TSQR
  };

#ifdef HAVE_XPETRA_EPETRA
#ifndef  EPETRA_NO_32BIT_GLOBAL_INDICES
  template<>
  class MultiVecTraits<double, Xpetra::MultiVector<double,int,int,Xpetra::EpetraNode> > {
  private:
    typedef double              Scalar;
    typedef int                 LO;
    typedef int                 GO;
    typedef Xpetra::EpetraNode  Node;

#ifdef HAVE_XPETRA_TPETRA // TODO check whether Tpetra is instantiated on all template parameters!
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
    typedef Xpetra::TpetraMultiVector<Scalar,LO,GO,Node>                    TpetraMultiVector;
    typedef MultiVecTraits<Scalar, Tpetra::MultiVector<Scalar,LO,GO,Node> > MultiVecTraitsTpetra;
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
    typedef Xpetra::EpetraMultiVectorT<GO,Node>            EpetraMultiVector;
    typedef MultiVecTraits   <Scalar, Epetra_MultiVector>  MultiVecTraitsEpetra;
#endif

  public:

#ifdef HAVE_BELOS_XPETRA_TIMERS
    static RCP<Teuchos::Time> mvTimesMatAddMvTimer_, mvTransMvTimer_;
#endif

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > Clone( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const int numvecs )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::Clone(toTpetra(mv), numvecs)));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return rcp(new EpetraMultiVector(MultiVecTraitsEpetra::Clone(toEpetra(mv), numvecs)));
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneCopy( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv))));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return rcp(new EpetraMultiVector(MultiVecTraitsEpetra::CloneCopy(toEpetra(mv))));
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneCopy( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv), index)));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return rcp(new EpetraMultiVector(MultiVecTraitsEpetra::CloneCopy(toEpetra(mv), index)));
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneCopy (const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const Teuchos::Range1D& index)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv), index)));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return rcp(new EpetraMultiVector(MultiVecTraitsEpetra::CloneCopy(toEpetra(mv), index)));
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneViewNonConst( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneViewNonConst(toTpetra(mv), index)));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return rcp(new EpetraMultiVector(MultiVecTraitsEpetra::CloneViewNonConst(toEpetra(mv), index)));
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneViewNonConst(Xpetra::MultiVector<Scalar,LO,GO,Node>& mv,
                      const Teuchos::Range1D& index)
    {

#ifdef HAVE_MUELUA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneViewNonConst(toTpetra(mv), index)));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return rcp(new EpetraMultiVector(MultiVecTraitsEpetra::CloneViewNonConst(toEpetra(mv), index)));
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<const Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneView(const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        //TODO: double check if the const_cast is safe here.
        RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> > r = MultiVecTraitsTpetra::CloneView(toTpetra(mv), index);
        return rcp(new TpetraMultiVector(Teuchos::rcp_const_cast<Tpetra::MultiVector<Scalar,LO,GO,Node> >(r)));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        //TODO: double check if the const_cast is safe here.
        RCP<const Epetra_MultiVector > r = MultiVecTraitsEpetra::CloneView(toEpetra(mv), index);
        return rcp(new EpetraMultiVector(Teuchos::rcp_const_cast<Epetra_MultiVector>(r)));
      }
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<const Xpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneView (const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const Teuchos::Range1D& index)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        //TODO: double check if the const_cast is safe here.
        RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> > r = MultiVecTraitsTpetra::CloneView(toTpetra(mv), index);
        return rcp(new TpetraMultiVector(Teuchos::rcp_const_cast<Tpetra::MultiVector<Scalar,LO,GO,Node> >(r)));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
     if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        //TODO: double check if the const_cast is safe here.
        RCP<const Epetra_MultiVector > r = MultiVecTraitsEpetra::CloneView(toEpetra(mv), index);
        return rcp(new EpetraMultiVector(Teuchos::rcp_const_cast<Epetra_MultiVector>(r)));
      }
#endif

      XPETRA_FACTORY_END;
    }

    static ptrdiff_t GetGlobalLength( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        return MultiVecTraitsTpetra::GetGlobalLength(toTpetra(mv));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return MultiVecTraitsEpetra::GetGlobalLength(toEpetra(mv));
#endif

      XPETRA_FACTORY_END;
    }

    static int GetNumberVecs( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        return MultiVecTraitsTpetra::GetNumberVecs(toTpetra(mv));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return MultiVecTraitsEpetra::GetNumberVecs(toEpetra(mv));
#endif

      XPETRA_FACTORY_END;
    }

    static bool HasConstantStride( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        return  MultiVecTraitsTpetra::HasConstantStride(toTpetra(mv));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return  MultiVecTraitsEpetra::HasConstantStride(toEpetra(mv));
#endif

      XPETRA_FACTORY_END;
    }

    static void MvTimesMatAddMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A,
                                 const Teuchos::SerialDenseMatrix<int,Scalar>& B,
                                 Scalar beta, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_BELOS_XPETRA_TIMERS
      Teuchos::TimeMonitor lcltimer(*mvTimesMatAddMvTimer_);
#endif

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        MultiVecTraitsTpetra::MvTimesMatAddMv(alpha, toTpetra(A), B, beta, toTpetra(mv));
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvTimesMatAddMv(alpha, toEpetra(A), B, beta, toEpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvAddMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, Scalar beta, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        MultiVecTraitsTpetra::MvAddMv(alpha, toTpetra(A), beta, toTpetra(B), toTpetra(mv));
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvAddMv(alpha, toEpetra(A), beta, toEpetra(B), toEpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvScale ( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, Scalar alpha )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        MultiVecTraitsTpetra::MvScale(toTpetra(mv), alpha);
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvScale(toEpetra(mv), alpha);
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvScale ( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<Scalar>& alphas )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        MultiVecTraitsTpetra::MvScale(toTpetra(mv), alphas);
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvScale(toEpetra(mv), alphas);
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvTransMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, Teuchos::SerialDenseMatrix<int,Scalar>& C)
    {

#ifdef HAVE_BELOS_XPETRA_TIMERS
      Teuchos::TimeMonitor lcltimer(*mvTransMvTimer_);
#endif

#ifdef HAVE_XPETRA_TPETRA
      if (A.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        MultiVecTraitsTpetra::MvTransMv(alpha, toTpetra(A), toTpetra(B), C);
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (A.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvTransMv(alpha, toEpetra(A), toEpetra(B), C);
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvDot( const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, std::vector<Scalar> &dots)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (A.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        MultiVecTraitsTpetra::MvDot(toTpetra(A), toTpetra(B), dots);
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (A.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvDot(toEpetra(A), toEpetra(B), dots);
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvNorm(const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::vector<Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec, NormType type=TwoNorm)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        MultiVecTraitsTpetra::MvNorm(toTpetra(mv), normvec, type);
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvNorm(toEpetra(mv), normvec, type);
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void SetBlock( const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const std::vector<int>& index, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        MultiVecTraitsTpetra::SetBlock(toTpetra(A), index, toTpetra(mv));
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::SetBlock(toEpetra(A), index, toEpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void
    SetBlock (const Xpetra::MultiVector<Scalar,LO,GO,Node>& A,
              const Teuchos::Range1D& index,
              Xpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        MultiVecTraitsTpetra::SetBlock(toTpetra(A), index, toTpetra(mv));
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::SetBlock(toEpetra(A), index, toEpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void
    Assign (const Xpetra::MultiVector<Scalar,LO,GO,Node>& A,
            Xpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        MultiVecTraitsTpetra::Assign(toTpetra(A), toTpetra(mv));
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::Assign(toEpetra(A), toEpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvRandom( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        MultiVecTraitsTpetra::MvRandom(toTpetra(mv));
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvRandom(toEpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvInit( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        MultiVecTraitsTpetra::MvInit(toTpetra(mv), alpha);
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvInit(toEpetra(mv), alpha);
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvPrint( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::ostream& os )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        MultiVecTraitsTpetra::MvPrint(toTpetra(mv), os);
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvPrint(toEpetra(mv), os);
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

#ifdef HAVE_BELOS_TSQR
#  if defined(HAVE_XPETRA_TPETRA) && \
  ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
   (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
    typedef BelosXpetraTsqrImpl::XpetraTpetraTsqrAdaptor<Scalar, LO, GO, Node> tsqr_adaptor_type;
#  else
    typedef BelosXpetraTsqrImpl::XpetraStubTsqrAdaptor<Scalar, LO, GO, Node> tsqr_adaptor_type;
#  endif
#endif // HAVE_BELOS_TSQR
  };
#endif // #ifndef  EPETRA_NO_32BIT_GLOBAL_INDICES
#endif // HAVE_XPETRA_EPETRA


#ifdef HAVE_XPETRA_EPETRA
#ifndef  EPETRA_NO_64BIT_GLOBAL_INDICES
  template<>
  class MultiVecTraits<double, Xpetra::MultiVector<double,int,long long,Xpetra::EpetraNode> > {
  private:
    typedef double              Scalar;
    typedef int                 LO;
    typedef long long           GO;
    typedef Xpetra::EpetraNode  Node;

#ifdef HAVE_XPETRA_TPETRA // TODO check whether Tpetra is instantiated on all template parameters!
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
    typedef Xpetra::TpetraMultiVector<Scalar,LO,GO,Node>                    TpetraMultiVector;
    typedef MultiVecTraits<Scalar, Tpetra::MultiVector<Scalar,LO,GO,Node> > MultiVecTraitsTpetra;
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
    typedef Xpetra::EpetraMultiVectorT<GO,Node>            EpetraMultiVector;
    typedef MultiVecTraits   <Scalar, Epetra_MultiVector>  MultiVecTraitsEpetra;
#endif

  public:

#ifdef HAVE_BELOS_XPETRA_TIMERS
    static RCP<Teuchos::Time> mvTimesMatAddMvTimer_, mvTransMvTimer_;
#endif

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > Clone( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const int numvecs )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::Clone(toTpetra(mv), numvecs)));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return rcp(new EpetraMultiVector(MultiVecTraitsEpetra::Clone(toEpetra(mv), numvecs)));
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneCopy( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv))));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return rcp(new EpetraMultiVector(MultiVecTraitsEpetra::CloneCopy(toEpetra(mv))));
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneCopy( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv), index)));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return rcp(new EpetraMultiVector(MultiVecTraitsEpetra::CloneCopy(toEpetra(mv), index)));
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneCopy (const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const Teuchos::Range1D& index)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv), index)));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return rcp(new EpetraMultiVector(MultiVecTraitsEpetra::CloneCopy(toEpetra(mv), index)));
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneViewNonConst( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneViewNonConst(toTpetra(mv), index)));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return rcp(new EpetraMultiVector(MultiVecTraitsEpetra::CloneViewNonConst(toEpetra(mv), index)));
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneViewNonConst(Xpetra::MultiVector<Scalar,LO,GO,Node>& mv,
                      const Teuchos::Range1D& index)
    {

#ifdef HAVE_MUELUA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneViewNonConst(toTpetra(mv), index)));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return rcp(new EpetraMultiVector(MultiVecTraitsEpetra::CloneViewNonConst(toEpetra(mv), index)));
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<const Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneView(const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        //TODO: double check if the const_cast is safe here.
        RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> > r = MultiVecTraitsTpetra::CloneView(toTpetra(mv), index);
        return rcp(new TpetraMultiVector(Teuchos::rcp_const_cast<Tpetra::MultiVector<Scalar,LO,GO,Node> >(r)));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        //TODO: double check if the const_cast is safe here.
        RCP<const Epetra_MultiVector > r = MultiVecTraitsEpetra::CloneView(toEpetra(mv), index);
        return rcp(new EpetraMultiVector(Teuchos::rcp_const_cast<Epetra_MultiVector>(r)));
      }
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<const Xpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneView (const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const Teuchos::Range1D& index)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        //TODO: double check if the const_cast is safe here.
        RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> > r = MultiVecTraitsTpetra::CloneView(toTpetra(mv), index);
        return rcp(new TpetraMultiVector(Teuchos::rcp_const_cast<Tpetra::MultiVector<Scalar,LO,GO,Node> >(r)));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
     if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        //TODO: double check if the const_cast is safe here.
        RCP<const Epetra_MultiVector > r = MultiVecTraitsEpetra::CloneView(toEpetra(mv), index);
        return rcp(new EpetraMultiVector(Teuchos::rcp_const_cast<Epetra_MultiVector>(r)));
      }
#endif

      XPETRA_FACTORY_END;
    }

    static ptrdiff_t GetGlobalLength( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        return MultiVecTraitsTpetra::GetGlobalLength(toTpetra(mv));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return MultiVecTraitsEpetra::GetGlobalLength(toEpetra(mv));
#endif

      XPETRA_FACTORY_END;
    }

    static int GetNumberVecs( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        return MultiVecTraitsTpetra::GetNumberVecs(toTpetra(mv));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return MultiVecTraitsEpetra::GetNumberVecs(toEpetra(mv));
#endif

      XPETRA_FACTORY_END;
    }

    static bool HasConstantStride( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        return  MultiVecTraitsTpetra::HasConstantStride(toTpetra(mv));
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return  MultiVecTraitsEpetra::HasConstantStride(toEpetra(mv));
#endif

      XPETRA_FACTORY_END;
    }

    static void MvTimesMatAddMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A,
                                 const Teuchos::SerialDenseMatrix<int,Scalar>& B,
                                 Scalar beta, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_BELOS_XPETRA_TIMERS
      Teuchos::TimeMonitor lcltimer(*mvTimesMatAddMvTimer_);
#endif

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        MultiVecTraitsTpetra::MvTimesMatAddMv(alpha, toTpetra(A), B, beta, toTpetra(mv));
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvTimesMatAddMv(alpha, toEpetra(A), B, beta, toEpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvAddMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, Scalar beta, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        MultiVecTraitsTpetra::MvAddMv(alpha, toTpetra(A), beta, toTpetra(B), toTpetra(mv));
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvAddMv(alpha, toEpetra(A), beta, toEpetra(B), toEpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvScale ( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, Scalar alpha )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        MultiVecTraitsTpetra::MvScale(toTpetra(mv), alpha);
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvScale(toEpetra(mv), alpha);
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvScale ( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<Scalar>& alphas )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        MultiVecTraitsTpetra::MvScale(toTpetra(mv), alphas);
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvScale(toEpetra(mv), alphas);
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvTransMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, Teuchos::SerialDenseMatrix<int,Scalar>& C)
    {

#ifdef HAVE_BELOS_XPETRA_TIMERS
      Teuchos::TimeMonitor lcltimer(*mvTransMvTimer_);
#endif

#ifdef HAVE_XPETRA_TPETRA
      if (A.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        MultiVecTraitsTpetra::MvTransMv(alpha, toTpetra(A), toTpetra(B), C);
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (A.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvTransMv(alpha, toEpetra(A), toEpetra(B), C);
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvDot( const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, std::vector<Scalar> &dots)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (A.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        MultiVecTraitsTpetra::MvDot(toTpetra(A), toTpetra(B), dots);
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (A.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvDot(toEpetra(A), toEpetra(B), dots);
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvNorm(const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::vector<Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec, NormType type=TwoNorm)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        MultiVecTraitsTpetra::MvNorm(toTpetra(mv), normvec, type);
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvNorm(toEpetra(mv), normvec, type);
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void SetBlock( const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const std::vector<int>& index, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        MultiVecTraitsTpetra::SetBlock(toTpetra(A), index, toTpetra(mv));
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::SetBlock(toEpetra(A), index, toEpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void
    SetBlock (const Xpetra::MultiVector<Scalar,LO,GO,Node>& A,
              const Teuchos::Range1D& index,
              Xpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        MultiVecTraitsTpetra::SetBlock(toTpetra(A), index, toTpetra(mv));
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::SetBlock(toEpetra(A), index, toEpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void
    Assign (const Xpetra::MultiVector<Scalar,LO,GO,Node>& A,
            Xpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        MultiVecTraitsTpetra::Assign(toTpetra(A), toTpetra(mv));
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::Assign(toEpetra(A), toEpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvRandom( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        MultiVecTraitsTpetra::MvRandom(toTpetra(mv));
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvRandom(toEpetra(mv));
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvInit( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        MultiVecTraitsTpetra::MvInit(toTpetra(mv), alpha);
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvInit(toEpetra(mv), alpha);
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

    static void MvPrint( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::ostream& os )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        MultiVecTraitsTpetra::MvPrint(toTpetra(mv), os);
        return;
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra) {
        MultiVecTraitsEpetra::MvPrint(toEpetra(mv), os);
        return;
      }
#endif

      XPETRA_FACTORY_END;
    }

#ifdef HAVE_BELOS_TSQR
#  if defined(HAVE_XPETRA_TPETRA) && \
  ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
   (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
    typedef BelosXpetraTsqrImpl::XpetraTpetraTsqrAdaptor<Scalar, LO, GO, Node> tsqr_adaptor_type;
#  else
    typedef BelosXpetraTsqrImpl::XpetraStubTsqrAdaptor<Scalar, LO, GO, Node> tsqr_adaptor_type;
#  endif
#endif // HAVE_BELOS_TSQR

  };
#endif // #ifndef  EPETRA_NO_64BIT_GLOBAL_INDICES
#endif // HAVE_XPETRA_EPETRA


} // end of Belos namespace

#endif
