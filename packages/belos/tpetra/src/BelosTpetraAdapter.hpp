// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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

#ifndef BELOS_TPETRA_ADAPTER_HPP
#define BELOS_TPETRA_ADAPTER_HPP

/*! \file BelosTpetraAdapter.hpp
    \brief Provides several interfaces between Belos virtual classes and Tpetra concrete classes.
*/

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_TestForException.hpp>

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosOperatorTraits.hpp"

namespace Belos {
 
  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::MultiVecTraits for Tpetra::MultiVector.
  //
  ////////////////////////////////////////////////////////////////////

  template<class Scalar, class LO, class GO>
  class MultiVecTraits<Scalar, Tpetra::MultiVector<Scalar,LO,GO> >
  {
  public:

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO> > Clone( const Tpetra::MultiVector<Scalar,LO,GO>& mv, const int numvecs )
    { 
      return Teuchos::rcp( new Tpetra::MultiVector<Scalar,LO,GO>(mv.getMap(),numvecs)); 
    }

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO> > CloneCopy( const Tpetra::MultiVector<Scalar,LO,GO>& mv )
    {
      return Teuchos::rcp( new Tpetra::MultiVector<Scalar,LO,GO>( mv ) ); 
    }

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO> > CloneCopy( const Tpetra::MultiVector<Scalar,LO,GO>& mv, const std::vector<int>& index )
    { 
      TEST_FOR_EXCEPTION(index.size() == 0,std::runtime_error,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneCopy(mv,index): numvecs must be greater than zero.");
#ifdef TPETRA_DEBUG
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::runtime_error,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneCopy(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.numVectors(), std::runtime_error,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneCopy(mv,index): indices must be < mv.numVectors().");
#endif
      for (typename std::vector<int>::size_type j=1; j<index.size(); ++j) {
        if (index[j] != index[j-1]+1) {
          // not contiguous; short circuit
          return mv.subCopy(Teuchos::arrayViewFromVector(index));
        }
      }
      // contiguous
      return mv.subCopy(Teuchos::Range1D(index.front(),index.back()));
    }

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO> > CloneView( Tpetra::MultiVector<Scalar,LO,GO>& mv, const std::vector<int>& index )
    {  
      TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): numvecs must be greater than zero.");
#ifdef TPETRA_DEBUG
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.numVectors(), std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be < mv.numVectors().");
#endif
      for (typename std::vector<int>::size_type j=1; j<index.size(); ++j) {
        if (index[j] != index[j-1]+1) {
          // not contiguous; short circuit
          return mv.subView(Teuchos::arrayViewFromVector(index));
        }
      }
      // contiguous
      return mv.subView(Teuchos::Range1D(index.front(),index.back()));
    }

    static Teuchos::RCP<const Tpetra::MultiVector<Scalar,LO,GO> > CloneView( const Tpetra::MultiVector<Scalar,LO,GO>& mv, const std::vector<int>& index )
    {
      TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): numvecs must be greater than zero.");
#ifdef TPETRA_DEBUG
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.numVectors(), std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be < mv.numVectors().");
#endif
      for (typename std::vector<int>::size_type j=1; j<index.size(); ++j) {
        if (index[j] != index[j-1]+1) {
          // not contiguous; short circuit
          return mv.subViewConst(Teuchos::arrayViewFromVector(index));
        }
      }
      // contiguous
      return mv.subViewConst(Teuchos::Range1D(index.front(),index.back()));
    }

    static int GetVecLength( const Tpetra::MultiVector<Scalar,LO,GO>& mv )
    { return mv.globalLength(); }

    static int GetNumberVecs( const Tpetra::MultiVector<Scalar,LO,GO>& mv )
    { return mv.numVectors(); }

    static void MvTimesMatAddMv( const Scalar alpha, const Tpetra::MultiVector<Scalar,LO,GO>& A, 
                                 const Teuchos::SerialDenseMatrix<int,Scalar>& B, 
                                 const Scalar beta, Tpetra::MultiVector<Scalar,LO,GO>& mv )
    {
      Tpetra::Map<LO,GO> LocalMap(B.numRows(), 0, A.getMap().getComm(),true);
      // TODO: this multivector should be a view of the data of B, not a copy
      Teuchos::ArrayView<const Scalar> Bvalues(B.values(),B.stride()*B.numCols());
      Tpetra::MultiVector<Scalar,LO,GO> B_mv(LocalMap,Bvalues,B.stride(),B.numCols());
      mv.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha, A, B_mv, beta);
    }

    static void MvAddMv( const Scalar alpha, const Tpetra::MultiVector<Scalar,LO,GO>& A, const Scalar beta, const Tpetra::MultiVector<Scalar,LO,GO>& B, Tpetra::MultiVector<Scalar,LO,GO>& mv )
    {
      mv.update(alpha,A,beta,B,Teuchos::ScalarTraits<Scalar>::zero());
    }

    static void MvScale ( Tpetra::MultiVector<Scalar,LO,GO>& mv, Scalar alpha )
    { mv.scale(alpha); }

    static void MvScale ( Tpetra::MultiVector<Scalar,LO,GO>& mv, const std::vector<Scalar>& alpha )
    { mv.scale(alpha); }

    static void MvTransMv( const Scalar alpha, const Tpetra::MultiVector<Scalar,LO,GO>& A, const Tpetra::MultiVector<Scalar,LO,GO>& mv, Teuchos::SerialDenseMatrix<int,Scalar>& B )
    { 
      Tpetra::Map<LO,GO> LocalMap(B.numRows(), 0, A.getMap().getComm(),true);
      // TODO: this multivector should be a view of the data of B, so we don't have to perform the copy afterwards
      Tpetra::MultiVector<Scalar,LO,GO> B_mv(LocalMap,B.numCols(),true);
      B_mv.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,alpha,A,mv,Teuchos::ScalarTraits<Scalar>::zero());
      Teuchos::ArrayView<Scalar> av(B.values(),B.stride()*B.numCols());
      B_mv.extractCopy1D(av,B.stride());
    }

    static void MvDot( const Tpetra::MultiVector<Scalar,LO,GO>& A, const Tpetra::MultiVector<Scalar,LO,GO>& B, std::vector<Scalar>& dots)
    { 
      TEST_FOR_EXCEPTION(A.numVectors() != B.numVectors(),std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvDot(A,B,dots): A and B must have the same number of vectors.");
#ifdef TPETRA_DEBUG
      TEST_FOR_EXCEPTION(dots.size() < (typename std::vector<int>::size_type)A.numVectors(),std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvDot(A,B,dots): dots must have room for all dot products.");
#endif
      Teuchos::ArrayView<Scalar> av(dots);
      A.dot(B,av(0,A.numVectors()));
    }

    static void MvNorm( const Tpetra::MultiVector<Scalar,LO,GO>& mv, std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec, NormType type=TwoNorm )
    { 
#ifdef TPETRA_DEBUG
      TEST_FOR_EXCEPTION(normvec.size() < (typename std::vector<int>::size_type)mv.numVectors(),std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvNorm(mv,normvec): normvec must have room for all norms.");
#endif
      Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> av(normvec);
      switch (type) {
        case OneNorm:
          mv.norm1(av(0,mv.numVectors()));
          break;
        case TwoNorm:
          mv.norm2(av(0,mv.numVectors()));
          break;
        case InfNorm:
          mv.normInf(av(0,mv.numVectors()));
          break;
      }
    }

    static void SetBlock( const Tpetra::MultiVector<Scalar,LO,GO>& A, const std::vector<int>& index, Tpetra::MultiVector<Scalar,LO,GO>& mv )
    {
#ifdef TPETRA_DEBUG
      TEST_FOR_EXCEPTION((typename std::vector<int>::size_type)A.numVectors() < index.size(),std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::SetBlock(A,index,mv): index must be the same size as A.");
#endif
      Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO> > mvsub = mv.subView(Teuchos::arrayViewFromVector(index));
      if ((typename std::vector<int>::size_type)A.numVectors() > index.size()) {
        Teuchos::RCP<const Tpetra::MultiVector<Scalar,LO,GO> > Asub = A.subViewConst(Teuchos::Range1D(0,index.size()-1));
        (*mvsub) = (*Asub);
      }
      else {
        (*mvsub) = A;
      }
      mvsub = Teuchos::null;
    }

    static void MvRandom( Tpetra::MultiVector<Scalar,LO,GO>& mv )
    { mv.random(); }

    static void MvInit( Tpetra::MultiVector<Scalar,LO,GO>& mv, Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() )
    { mv.putScalar(alpha); }

    static void MvPrint( const Tpetra::MultiVector<Scalar,LO,GO>& mv, std::ostream& os )
    { mv.print(os); }
  };        

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::OperatorTraits for Tpetra::Operator.
  //
  ////////////////////////////////////////////////////////////////////

  template <class Scalar,class LO, class GO> 
  class OperatorTraits < Scalar, Tpetra::MultiVector<Scalar,LO,GO>, Tpetra::Operator<Scalar,LO,GO> >
  {
  public:
    static void Apply ( const Tpetra::Operator<Scalar,LO,GO> & Op, 
                        const Tpetra::MultiVector<Scalar,LO,GO> & X,
                              Tpetra::MultiVector<Scalar,LO,GO> & Y,
                        ETrans trans=NOTRANS )
    { 
      TEST_FOR_EXCEPTION(trans != NOTRANS, std::logic_error, "Feature not yet implemented.");
      Op.apply(X,Y,Teuchos::NO_TRANS);
    }
  };

} // end of Belos namespace 

#endif 
// end of file BELOS_TPETRA_ADAPTER_HPP
