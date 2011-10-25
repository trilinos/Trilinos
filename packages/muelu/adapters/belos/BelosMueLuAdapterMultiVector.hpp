#ifndef BELOS_MUELU_MULTIVECTOR_ADAPTER_HPP
#define BELOS_MUELU_MULTIVECTOR_ADAPTER_HPP

#include <Xpetra_ConfigDefs.hpp>

/*! \file  BelosXpetraAdapter.hpp
    \brief Provides several interfaces between Belos virtual classes and Xpetra concrete classes.
*/

#include <Xpetra_MultiVector.hpp>
 
#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraMultiVector.hpp>
#endif

#include <BelosConfigDefs.hpp>
#include <BelosTypes.hpp>
#include <BelosMultiVecTraits.hpp>
#include <BelosOperatorTraits.hpp>

//TODO: missing headers Belos/Tpetra and Belos/Epetra adapters?

namespace Belos { // should be moved to Belos or Xpetra?

  using Teuchos::RCP;

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
  class MultiVecTraits<Scalar, Xpetra::MultiVector<Scalar,LO,GO,Node> >
  {

  private:
    typedef Xpetra::TpetraMultiVector<Scalar,LO,GO,Node> TpetraMultiVector;
    typedef MultiVecTraits<Scalar, Tpetra::MultiVector<Scalar,LO,GO,Node> > MultiVecTraitsTpetra;

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
      if (mv.getMap()->lib() == Xpetra::UseTpetra) 
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneView(toTpetra(mv), index)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static int GetVecLength( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {
 
#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) 
        return MultiVecTraitsTpetra::GetVecLength(toTpetra(mv));
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

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;

 
#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) 
        MultiVecTraitsTpetra::MvTimesMatAddMv(alpha, toTpetra(A), B, beta, toTpetra(mv));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvAddMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, Scalar beta, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {
 
#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) 
        MultiVecTraitsTpetra::MvAddMv(alpha, toTpetra(A), beta, toTpetra(B), toTpetra(mv));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvScale ( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, Scalar alpha )
    { 
 
#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) 
        MultiVecTraitsTpetra::MvScale(toTpetra(mv) ,alpha);
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvScale ( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<Scalar>& alphas )
    { 
 
#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) 
        MultiVecTraitsTpetra::MvScale(toTpetra(mv) ,alphas);
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvTransMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, Teuchos::SerialDenseMatrix<int,Scalar>& C)
    { 
 
#ifdef HAVE_BELOS_XPETRA_TIMERS
      Teuchos::TimeMonitor lcltimer(*mvTransMvTimer_);
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(A.getMap()->lib());
      XPETRA_FACTORY_END;
 
#ifdef HAVE_XPETRA_TPETRA
      if (A.getMap()->lib() == Xpetra::UseTpetra) 
        MultiVecTraitsTpetra::MvTransMv(alpha, toTpetra(A), toTpetra(B), C);
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(A.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvDot( const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, std::vector<Scalar> &dots)
    {
 
#ifdef HAVE_XPETRA_TPETRA
      if (A.getMap()->lib() == Xpetra::UseTpetra) 
        MultiVecTraitsTpetra::MvDot(toTpetra(A), toTpetra(B), dots);
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(A.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvNorm(const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec, NormType type=TwoNorm)
    { 
 
#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) 
        MultiVecTraitsTpetra::MvNorm(toTpetra(mv), normvec, type);
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void SetBlock( const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const std::vector<int>& index, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {
 
#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) 
        MultiVecTraitsTpetra::SetBlock(toTpetra(A), index, toTpetra(mv));
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
      if (mv.getMap()->lib() == Xpetra::UseTpetra) 
        MultiVecTraitsTpetra::SetBlock(toTpetra(A), index, toTpetra(mv));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void
    Assign (const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, 
            Xpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
 
#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) 
        MultiVecTraitsTpetra::Assign(toTpetra(A), toTpetra(mv));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvRandom( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {
 
#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) 
        MultiVecTraitsTpetra::MvRandom(toTpetra(mv));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvInit( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() )
    { 
 
#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) 
        MultiVecTraitsTpetra::MvInit(toTpetra(mv), alpha);
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

    static void MvPrint( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::ostream& os )
    { 
 
#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra) 
        MultiVecTraitsTpetra::MvPrint(toTpetra(mv), os);
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(mv.getMap()->lib());
      XPETRA_FACTORY_END;
    }

  };        

} // end of Belos namespace 

#endif
