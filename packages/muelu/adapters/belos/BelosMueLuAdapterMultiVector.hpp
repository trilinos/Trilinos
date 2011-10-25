#ifndef BELOS_MUELU_MULTIVECTOR_ADAPTER_HPP
#define BELOS_MUELU_MULTIVECTOR_ADAPTER_HPP

/*! \file BelosXpetraAdapter.hpp
  \brief Provides several interfaces between Belos virtual classes and Xpetra concrete classes.
*/

#include <Xpetra_MultiVector.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosTypes.hpp>
#include <BelosMultiVecTraits.hpp>
#include <BelosOperatorTraits.hpp>

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
    typedef MultiVecTraits<Scalar, Tpetra::MultiVector<Scalar,LO,GO,Node> > MultiVecTraitsTpetra;

  public:
#ifdef HAVE_BELOS_XPETRA_TIMERS
    static RCP<Teuchos::Time> mvTimesMatAddMvTimer_, mvTransMvTimer_;
#endif

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > Clone( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const int numvecs )
    { 
      return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::Clone(toTpetra(mv), numvecs)));
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneCopy( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {
      return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv))));
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneCopy( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    { 
      return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv), index)));
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > 
    CloneCopy (const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, 
               const Teuchos::Range1D& index)
    { 
      return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv), index)));
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneViewNonConst( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {
      return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneViewNonConst(toTpetra(mv), index)));
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > 
    CloneViewNonConst(Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, 
                      const Teuchos::Range1D& index)
    {
      return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneViewNonConst(toTpetra(mv), index)));
    }

    static RCP<const Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneView(const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {
      //TODO: double check if the const_cast is safe here.
      RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> > r = MultiVecTraitsTpetra::CloneView(toTpetra(mv), index);
      return rcp(new TpetraMultiVector(Teuchos::rcp_const_cast<Tpetra::MultiVector<Scalar,LO,GO,Node> >(r)));
    }

    static RCP<const Xpetra::MultiVector<Scalar,LO,GO,Node> > 
    CloneView (const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, 
               const Teuchos::Range1D& index)
    {
      return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneView(toTpetra(mv), index)));
    }

    static int GetVecLength( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {
      return MultiVecTraitsTpetra::GetVecLength(toTpetra(mv));
    }

    static int GetNumberVecs( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    { 
      return MultiVecTraitsTpetra::GetNumberVecs(toTpetra(mv));
    }

    static bool HasConstantStride( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    { 
      return MultiVecTraitsTpetra::HasConstantStride(toTpetra(mv));
    }

    static void MvTimesMatAddMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, 
                                 const Teuchos::SerialDenseMatrix<int,Scalar>& B, 
                                 Scalar beta, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {
#ifdef HAVE_BELOS_XPETRA_TIMERS
      Teuchos::TimeMonitor lcltimer(*mvTimesMatAddMvTimer_);
#endif

      MultiVecTraitsTpetra::MvTimesMatAddMv(alpha, toTpetra(A), B, beta, toTpetra(mv));
    }

    static void MvAddMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, Scalar beta, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

      MultiVecTraitsTpetra::MvAddMv(alpha, toTpetra(A), beta, toTpetra(B), toTpetra(mv));
    }

    static void MvScale ( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, Scalar alpha )
    { 
      MultiVecTraitsTpetra::MvScale(toTpetra(mv) ,alpha);
    }

    static void MvScale ( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<Scalar>& alphas )
    { 
      MultiVecTraitsTpetra::MvScale(toTpetra(mv) ,alphas);
    }

    static void MvTransMv( Scalar alpha, const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, Teuchos::SerialDenseMatrix<int,Scalar>& C)
    { 
#ifdef HAVE_BELOS_XPETRA_TIMERS
      Teuchos::TimeMonitor lcltimer(*mvTransMvTimer_);
#endif

      MultiVecTraitsTpetra::MvTransMv(alpha, toTpetra(A), toTpetra(B), C);
    }

    static void MvDot( const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const Xpetra::MultiVector<Scalar,LO,GO,Node>& B, std::vector<Scalar> &dots)
    {
      MultiVecTraitsTpetra::MvDot(toTpetra(A), toTpetra(B), dots);
    }

    static void MvNorm(const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec, NormType type=TwoNorm)
    { 
      MultiVecTraitsTpetra::MvNorm(toTpetra(mv), normvec, type);
    }

    static void SetBlock( const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, const std::vector<int>& index, Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {
      MultiVecTraitsTpetra::SetBlock(toTpetra(A), index, toTpetra(mv));
    }

    static void
    SetBlock (const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, 
              const Teuchos::Range1D& index, 
              Xpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      MultiVecTraitsTpetra::SetBlock(toTpetra(A), index, toTpetra(mv));
    }

    static void
    Assign (const Xpetra::MultiVector<Scalar,LO,GO,Node>& A, 
            Xpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      MultiVecTraitsTpetra::Assign(toTpetra(A), toTpetra(mv));
    }

    static void MvRandom( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {
      MultiVecTraitsTpetra::MvRandom(toTpetra(mv));
    }

    static void MvInit( Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() )
    { 
      MultiVecTraitsTpetra::MvInit(toTpetra(mv), alpha);
    }

    static void MvPrint( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::ostream& os )
    { 
      MultiVecTraitsTpetra::MvPrint(toTpetra(mv), os);
    }

  };        

} // end of Belos namespace 

#endif // end of file BELOS_MUELU_MULTIVECTOR_ADAPTER_HPP
