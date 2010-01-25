// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef ANASAZI_TPETRA_ADAPTER_HPP
#define ANASAZI_TPETRA_ADAPTER_HPP

/*! \file AnasaziTpetraAdapter.hpp
  \brief Specializaitons of Anasazi multi-vector and operator traits classes for the Tpetra MultiVector and Operator classes.
*/

// TODO: the assumption is made that the solver, multivector and operator are templated on the same scalar. this will need to be modified.

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Array.hpp>

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"

namespace Anasazi {

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::MultiVecTraits for Tpetra::MultiVector.
  //
  ////////////////////////////////////////////////////////////////////

  /*!  \brief Template specialization of Anasazi::MultiVecTraits class using the Tpetra::MultiVector class.

    This interface will ensure that any Tpetra::MultiVector will be accepted by the Anasazi
    templated solvers.  */
  template<class Scalar, class LO, class GO, class Node>
  class MultiVecTraits<Scalar, Tpetra::MultiVector<Scalar,LO,GO,Node> >
  {
  public:

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> > Clone( const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, const int numvecs )
    { 
      return Teuchos::rcp( new Tpetra::MultiVector<Scalar,LO,GO,Node>(mv.getMap(),numvecs));
    }

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> > CloneCopy( const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {
      return Teuchos::rcp( new Tpetra::MultiVector<Scalar,LO,GO,Node>( mv ) ); 
    }

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> > CloneCopy( const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    { 
      TEST_FOR_EXCEPTION(index.size() == 0,std::runtime_error,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::Clone(mv,index): numvecs must be greater than zero.");
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::runtime_error,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::Clone(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.getNumVectors(), std::runtime_error,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::Clone(mv,index): indices must be < mv.getNumVectors().");
#endif
      for (typename std::vector<int>::size_type j=1; j<index.size(); ++j) {
        if (index[j] != index[j-1]+1) {
          // not contiguous; short circuit
          Teuchos::Array<size_t> stinds(index.begin(), index.end());
          return mv.subCopy(stinds());
        }
      }
      // contiguous
      return mv.subCopy(Teuchos::Range1D(index.front(),index.back()));
    }

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> > CloneView( Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {
      TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): numvecs must be greater than zero.");
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.getNumVectors(), std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be < mv.getNumVectors().");
#endif
      for (typename std::vector<int>::size_type j=1; j<index.size(); ++j) {
        if (index[j] != index[j-1]+1) {
          // not contiguous; short circuit
          Teuchos::Array<size_t> stinds(index.begin(), index.end());
          return mv.subViewNonConst(stinds);
        }
      }
      // contiguous
      return mv.subViewNonConst(Teuchos::Range1D(index.front(),index.back()));
    }

    static Teuchos::RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> > CloneView( const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {
      TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): numvecs must be greater than zero.");
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.getNumVectors(), std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be < mv.getNumVectors().");
#endif
      for (typename std::vector<int>::size_type j=1; j<index.size(); ++j) {
        if (index[j] != index[j-1]+1) {
          // not contiguous; short circuit
          Teuchos::Array<size_t> stinds(index.begin(), index.end());
          return mv.subView(stinds);
        }
      }
      // contiguous
      return mv.subView(Teuchos::Range1D(index.front(),index.back()));
    }

    static int GetVecLength( const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    { return mv.getGlobalLength(); }

    static int GetNumberVecs( const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    { return mv.getNumVectors(); }

    static void MvTimesMatAddMv( Scalar alpha, const Tpetra::MultiVector<Scalar,LO,GO,Node>& A, 
                                 const Teuchos::SerialDenseMatrix<int,Scalar>& B, 
                                 Scalar beta, Tpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {
      // create local map
      Tpetra::Map<LO,GO,Node> LocalMap(B.numRows(), 0, A.getMap()->getComm(), Tpetra::LocallyReplicated);
      // encapsulate Teuchos::SerialDenseMatrix data in ArrayView
      Teuchos::ArrayView<const Scalar> Bvalues(B.values(),B.stride()*B.numCols());
      // create locally replicated MultiVector with a copy of this data
      Tpetra::MultiVector<Scalar,LO,GO,Node> B_mv(Teuchos::rcpFromRef(LocalMap),Bvalues,B.stride(),B.numCols());
      // multiply
      mv.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha, A, B_mv, beta);
    }

    static void MvAddMv( Scalar alpha, const Tpetra::MultiVector<Scalar,LO,GO,Node>& A, Scalar beta, const Tpetra::MultiVector<Scalar,LO,GO,Node>& B, Tpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {
      mv.update(alpha,A,beta,B,Teuchos::ScalarTraits<Scalar>::zero());
    }

    static void MvScale ( Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, Scalar alpha )
    { mv.scale(alpha); }

    static void MvScale ( Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<Scalar>& alphas )
    { mv.scale(alphas); }

    static void MvTransMv( Scalar alpha, const Tpetra::MultiVector<Scalar,LO,GO,Node>& A, const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, Teuchos::SerialDenseMatrix<int,Scalar>& B)
    { 
      // create local map
      Tpetra::Map<LO,GO,Node> LocalMap(B.numRows(), 0, A.getMap()->getComm(), Tpetra::LocallyReplicated);
      // create local multivector to hold the result
      Tpetra::MultiVector<Scalar,LO,GO,Node> B_mv(Teuchos::rcpFromRef(LocalMap),B.numCols(),true);
      // multiply result into local multivector
      B_mv.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,alpha,A,mv,Teuchos::ScalarTraits<Scalar>::zero());
      // create arrayview encapsulating the Teuchos::SerialDenseMatrix
      Teuchos::ArrayView<Scalar> av(B.values(),B.stride()*B.numCols());
      // extract a copy of the result into the array view (and therefore, the SerialDenseMatrix)
      B_mv.get1dCopy(av,B.stride());
    }

    static void MvDot( const Tpetra::MultiVector<Scalar,LO,GO,Node>& A, const Tpetra::MultiVector<Scalar,LO,GO,Node>& B, std::vector<Scalar> &dots)
    {
      TEST_FOR_EXCEPTION(A.getNumVectors() != B.getNumVectors(),std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvDot(A,B,dots): A and B must have the same number of vectors.");
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(dots.size() < (typename std::vector<int>::size_type)A.getNumVectors(),std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvDot(A,B,dots): dots must have room for all dot products.");
#endif
      Teuchos::ArrayView<Scalar> av(dots);
      A.dot(B,av(0,A.getNumVectors()));
    }

    static void MvNorm(const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec)
    { 
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(normvec.size() < (typename std::vector<int>::size_type)mv.getNumVectors(),std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvNorm(mv,normvec): normvec must have room for all norms.");
#endif
      Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> av(normvec);
      mv.norm2(av(0,mv.getNumVectors()));
    }

    static void SetBlock( const Tpetra::MultiVector<Scalar,LO,GO,Node>& A, const std::vector<int>& index, Tpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION((typename std::vector<int>::size_type)A.getNumVectors() < index.size(),std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::SetBlock(A,index,mv): index must be the same size as A.");
#endif
      Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> > mvsub = CloneView(mv,index);
      if ((typename std::vector<int>::size_type)A.getNumVectors() > index.size()) {
        Teuchos::RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> > Asub = A.subView(Teuchos::Range1D(0,index.size()-1));
        (*mvsub) = (*Asub);
      }
      else {
        (*mvsub) = A;
      }
      mvsub = Teuchos::null;
    }

    static void MvRandom( Tpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    { mv.randomize(); }

    static void MvInit( Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() )
    { mv.putScalar(alpha); }

    static void MvPrint( const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::ostream& os )
    { mv.print(os); }

  };        

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::OperatorTraits for Tpetra::Operator.
  //
  ////////////////////////////////////////////////////////////////////

  template <class Scalar, class LO, class GO, class Node> 
  class OperatorTraits < Scalar, Tpetra::MultiVector<Scalar,LO,GO,Node>, Tpetra::Operator<Scalar,LO,GO,Node> >
  {
  public:
    static void Apply ( const Tpetra::Operator<Scalar,LO,GO,Node> & Op, 
                        const Tpetra::MultiVector<Scalar,LO,GO,Node> & X,
                              Tpetra::MultiVector<Scalar,LO,GO,Node> & Y)
    { 
      Op.apply(X,Y,Teuchos::NO_TRANS);
    }
  };

} // end of Anasazi namespace 

#endif 
// end of file ANASAZI_TPETRA_ADAPTER_HPP
