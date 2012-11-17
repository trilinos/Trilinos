// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef BELOS_XPETRA_ADAPTER_MULTIVECTOR_HPP
#define BELOS_XPETRA_ADAPTER_MULTIVECTOR_HPP

#include <Xpetra_ConfigDefs.hpp>

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
#endif

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
#ifdef HAVE_XPETRA_TPETRA
    typedef Xpetra::TpetraMultiVector<Scalar,LO,GO,Node> TpetraMultiVector;
    typedef MultiVecTraits<Scalar, Tpetra::MultiVector<Scalar,LO,GO,Node> > MultiVecTraitsTpetra;
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

  };

  template<>
  class MultiVecTraits<double, Xpetra::MultiVector<double,int,int,Kokkos::DefaultNode::DefaultNodeType> >
  {

  private:

    typedef double Scalar;
    typedef int LO;
    typedef int GO;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;

#ifdef HAVE_XPETRA_TPETRA
    typedef Xpetra::TpetraMultiVector<Scalar,LO,GO,Node>                    TpetraMultiVector;
    typedef MultiVecTraits<Scalar, Tpetra::MultiVector<Scalar,LO,GO,Node> > MultiVecTraitsTpetra;
#endif

#ifdef HAVE_XPETRA_EPETRA
    typedef Xpetra::EpetraMultiVector                  EpetraMultiVector;
    typedef MultiVecTraits<Scalar, Epetra_MultiVector> MultiVecTraitsEpetra;
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

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return rcp(new EpetraMultiVector(MultiVecTraitsEpetra::Clone(toEpetra(mv), numvecs)));
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Xpetra::MultiVector<Scalar,LO,GO,Node> > CloneCopy( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv))));
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
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv), index)));
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
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneCopy(toTpetra(mv), index)));
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
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneViewNonConst(toTpetra(mv), index)));
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

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return rcp(new TpetraMultiVector(MultiVecTraitsTpetra::CloneViewNonConst(toTpetra(mv), index)));
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
        //TODO: double check if the const_cast is safe here.
        RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> > r = MultiVecTraitsTpetra::CloneView(toTpetra(mv), index);
        return rcp(new TpetraMultiVector(Teuchos::rcp_const_cast<Tpetra::MultiVector<Scalar,LO,GO,Node> >(r)));
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
        //TODO: double check if the const_cast is safe here.
        RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> > r = MultiVecTraitsTpetra::CloneView(toTpetra(mv), index);
        return rcp(new TpetraMultiVector(Teuchos::rcp_const_cast<Tpetra::MultiVector<Scalar,LO,GO,Node> >(r)));
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

    static int GetVecLength( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return MultiVecTraitsTpetra::GetVecLength(toTpetra(mv));
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (mv.getMap()->lib() == Xpetra::UseEpetra)
        return MultiVecTraitsEpetra::GetVecLength(toEpetra(mv));
#endif

      XPETRA_FACTORY_END;
    }

    static int GetNumberVecs( const Xpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {

#ifdef HAVE_XPETRA_TPETRA
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return MultiVecTraitsTpetra::GetNumberVecs(toTpetra(mv));
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
      if (mv.getMap()->lib() == Xpetra::UseTpetra)
        return  MultiVecTraitsTpetra::HasConstantStride(toTpetra(mv));
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
        MultiVecTraitsTpetra::MvTimesMatAddMv(alpha, toTpetra(A), B, beta, toTpetra(mv));
        return;
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
        MultiVecTraitsTpetra::MvAddMv(alpha, toTpetra(A), beta, toTpetra(B), toTpetra(mv));
        return;
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
        MultiVecTraitsTpetra::MvScale(toTpetra(mv), alpha);
        return;
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
        MultiVecTraitsTpetra::MvScale(toTpetra(mv), alphas);
        return;
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
        MultiVecTraitsTpetra::MvTransMv(alpha, toTpetra(A), toTpetra(B), C);
        return;
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
        MultiVecTraitsTpetra::MvDot(toTpetra(A), toTpetra(B), dots);
        return;
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
        MultiVecTraitsTpetra::MvNorm(toTpetra(mv), normvec, type);
        return;
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
        MultiVecTraitsTpetra::SetBlock(toTpetra(A), index, toTpetra(mv));
        return;
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
        MultiVecTraitsTpetra::SetBlock(toTpetra(A), index, toTpetra(mv));
        return;
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
        MultiVecTraitsTpetra::Assign(toTpetra(A), toTpetra(mv));
        return;
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
        MultiVecTraitsTpetra::MvRandom(toTpetra(mv));
        return;
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
        MultiVecTraitsTpetra::MvInit(toTpetra(mv), alpha);
        return;
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
        MultiVecTraitsTpetra::MvPrint(toTpetra(mv), os);
        return;
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

  };

} // end of Belos namespace

#endif
