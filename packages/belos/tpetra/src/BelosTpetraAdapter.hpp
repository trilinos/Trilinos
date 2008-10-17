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

  template<class Scalar>
  class MultiVecTraits<Scalar, Tpetra::MultiVector<int,Scalar> >
  {
  public:

    static Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > Clone( const Tpetra::MultiVector<int,Scalar>& mv, const int numvecs )
    { return Teuchos::rcp( new Tpetra::MultiVector<int,Scalar>(mv.getMap(),numvecs)); }

    static Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > CloneCopy( const Tpetra::MultiVector<int,Scalar>& mv )
    {
      return Teuchos::rcp( new Tpetra::MultiVector<int,Scalar>( mv ) ); 
    }

    static Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > CloneCopy( const Tpetra::MultiVector<int,Scalar>& mv, const std::vector<int>& index )
    { 
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(index.size() == 0,std::runtime_error,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::Clone(mv,index): numvecs must be greater than zero.");
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::runtime_error,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::Clone(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.numVectors(), std::runtime_error,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::Clone(mv,index): indices must be < mv.numVectors().");
#endif
      // TODO: check to see if index is contiguous, if so, use Teuchos::Range1D
      return mv.subCopy(Teuchos::arrayViewFromVector(index));
    }

    static Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > CloneView( Tpetra::MultiVector<int,Scalar>& mv, const std::vector<int>& index )
    {  
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): numvecs must be greater than zero.");
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.numVectors(), std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be < mv.numVectors().");
#endif
      return mv.subView(Teuchos::arrayViewFromVector(index));
    }

    static Teuchos::RCP<const Tpetra::MultiVector<int,Scalar> > CloneView( const Tpetra::MultiVector<int,Scalar>& mv, const std::vector<int>& index )
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): numvecs must be greater than zero.");
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.numVectors(), std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be < mv.numVectors().");
#endif
      return mv.subViewConst(Teuchos::arrayViewFromVector(index));
    }

    static int GetVecLength( const Tpetra::MultiVector<int,Scalar>& mv )
    { mv.globalLength(); }

    static int GetNumberVecs( const Tpetra::MultiVector<int,Scalar>& mv )
    { return mv.numVectors(); }

    static void MvTimesMatAddMv( const Scalar alpha, const Tpetra::MultiVector<int,Scalar>& A, 
                                 const Teuchos::SerialDenseMatrix<int,Scalar>& B, 
                                 const Scalar beta, Tpetra::MultiVector<int,Scalar>& mv )
    {
      Tpetra::Map<int> LocalMap(B.numRows(), 0, *A.getMap().getPlatform(),true);
      // TODO: this multivector should be a view of the data of B, not a copy
      Teuchos::ArrayView<const Scalar> Bvalues(B.values(),B.stride()*B.numCols());
      Tpetra::MultiVector<int,Scalar> B_mv(LocalMap,Bvalues,B.stride(),B.numCols());
      mv.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha, A, B_mv, beta);
    }

    static void MvAddMv( const Scalar alpha, const Tpetra::MultiVector<int,Scalar>& A, const Scalar beta, const Tpetra::MultiVector<int,Scalar>& B, Tpetra::MultiVector<int,Scalar>& mv )
    {
      mv.update(alpha,A,beta,B,Teuchos::ScalarTraits<Scalar>::zero());
    }

    static void MvScale ( Tpetra::MultiVector<int,Scalar>& mv, Scalar alpha )
    { TEST_FOR_EXCEPT(true); }

    static void MvScale ( Tpetra::MultiVector<int,Scalar>& mv, const std::vector<Scalar>& alpha )
    { TEST_FOR_EXCEPT(true); }

    static void MvTransMv( const Scalar alpha, const Tpetra::MultiVector<int,Scalar>& A, const Tpetra::MultiVector<int,Scalar>& mv, Teuchos::SerialDenseMatrix<int,Scalar>& B )
    { 
      Tpetra::Map<int> LocalMap(B.numRows(), 0, *A.getMap().getPlatform(),true);
      // TODO: this multivector should be a view of the data of B, so we don't have to perform the copy afterwards
      Tpetra::MultiVector<int,Scalar> B_mv(LocalMap,B.numCols(),true);
      B_mv.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,alpha,A,mv,Teuchos::ScalarTraits<Scalar>::zero());
      // copy the data from B_mv to B
      if (B.stride() == B.numRows()) {
        // data in B_mv,B is packed, so we can perform a single copy
        int dummy;
        Teuchos::ArrayView<Scalar> av(B.values(),B.numRows()*B.numCols());
        B_mv.extractCopy(av,dummy);
      }
      else {
        // data in B is not-packed, must perform multiple copies
        // build an array of views into B
        Teuchos::Array<Teuchos::ArrayView<Scalar> > aoa(B.numCols(),Teuchos::null);
        for (int j=0; j<B.numCols(); ++j) {
          aoa[j] = Teuchos::arrayView(B[j],B.numRows());
        }
        // get Tpetra::MultiVector to put data at these addresses
        B_mv.extractCopy(aoa);
      }
    }

    static void MvDot( const Tpetra::MultiVector<int,Scalar>& A, const Tpetra::MultiVector<int,Scalar>& B, std::vector<Scalar>& dots)
    { 
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(A.numVectors() != B.numVectors(),std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvDot(A,B,dots): A and B must have the same number of vectors.");
      TEST_FOR_EXCEPTION(dots.size() < (unsigned int)A.numVectors(),std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvDot(A,B,dots): dots must have room for all dot products.");
#endif
      Teuchos::ArrayView<Scalar> av(dots);
      A.dot(B,av(0,A.numVectors()));
    }

    static void MvNorm( const Tpetra::MultiVector<int,Scalar>& mv, std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec, NormType type=TwoNorm )
    { 
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

    static void SetBlock( const Tpetra::MultiVector<int,Scalar>& A, const std::vector<int>& index, Tpetra::MultiVector<int,Scalar>& mv )
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION((typename std::vector<int>::size_type)A.numVectors() < index.size(),std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::SetBlock(A,index,mv): index must be the same size as A.");
#endif
      Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > mvsub = mv.subView(Teuchos::arrayViewFromVector(index));
      *mvsub = A;
      mvsub = Teuchos::null;
    }

    static void MvRandom( Tpetra::MultiVector<int,Scalar>& mv )
    { mv.random(); }

    static void MvInit( Tpetra::MultiVector<int,Scalar>& mv, Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() )
    { mv.putScalar(alpha); }

    static void MvPrint( const Tpetra::MultiVector<int,Scalar>& mv, std::ostream& os )
    { TEST_FOR_EXCEPT(true); }
  };        

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::OperatorTraits for Tpetra::Operator.
  //
  ////////////////////////////////////////////////////////////////////

  template <class Scalar> 
  class OperatorTraits < Scalar, Tpetra::MultiVector<int,Scalar>, Tpetra::Operator<int,Scalar> >
  {
  public:
    static void Apply ( const Tpetra::Operator<int,Scalar>& Op, 
                        const Tpetra::MultiVector<int,Scalar>& X,
                              Tpetra::MultiVector<int,Scalar>& Y,
                        ETrans trans=NOTRANS )
    { 
      TEST_FOR_EXCEPTION(trans != NOTRANS, std::logic_error, "Feature not yet implemented.");
      Op.apply(X,Y,Teuchos::NO_TRANS);
    }
  };

} // end of Belos namespace 

#endif 
// end of file BELOS_TPETRA_ADAPTER_HPP
