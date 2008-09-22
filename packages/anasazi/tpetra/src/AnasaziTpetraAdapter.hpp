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

/*! \file AnasaziTpetraAdapter.hpp
  \brief Specializaitons of Anasazi multi-vector and operator traits classes for the Tpetra MultiVector and Operator classes.
*/

#ifndef ANASAZI_TPETRA_ADAPTER_HPP
#define ANASAZI_TPETRA_ADAPTER_HPP

#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziConfigDefs.hpp"

#include "AnasaziTypes.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziOperator.hpp"

#include <Teuchos_ScalarTraits.hpp>
#include <Tpetra_Operator.hpp>
#include <Tpetra_MultiVector.hpp>

namespace Anasazi {

  /////////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::MultiVecTraits for Tpetra::MultiVector
  //
  /////////////////////////////////////////////////////////////////////////

  /*! 
    \brief Template specialization of Anasazi::MultiVecTraits class using the Tpetra::MultiVector class.

    This interface will ensure that any Tpetra::MultiVector will be accepted by the Anasazi
    templated solvers.  
  */

  template<typename Scalar>
  class MultiVecTraits< Scalar, Tpetra::MultiVector<int,Scalar> >
  {
  public:

    //! @name Creation methods
    //@{ 

    /*! \brief Creates a new empty Tpetra::MultiVector<int,Scalar> containing \c numvecs columns.
      
    \return Reference-counted pointer to the new Tpetra::MultiVector<int,Scalar>.
    */
    static Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > Clone( const Tpetra::MultiVector<int,Scalar>& mv, const int numvecs )
    { 
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(numvecs <= 0,std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::Clone(mv,numvecs): numvecs must be greater than zero.");
#endif
      return Teuchos::rcp( new Tpetra::MultiVector<int,Scalar>(mv.getMap(), numvecs) ); 
    }

    /*! \brief Creates a new Tpetra::MultiVector<int,Scalar> and copies contents of \c mv into the new vector (deep copy).
      
      \return Reference-counted pointer to the new Tpetra::MultiVector<int,Scalar>.
    */
    static Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > CloneCopy( const Tpetra::MultiVector<int,Scalar>& mv )
    { 
      return Teuchos::rcp( new Tpetra::MultiVector<int,Scalar>( mv ) ); 
    }

    /*! \brief Creates a new Tpetra::MultiVector<int,Scalar> and copies the selected contents of \c mv into the new vector (deep copy).  

      The copied vectors from \c mv are indicated by the \c indeX.size() indices in \c index.      
      \return Reference-counted pointer to the new Tpetra::MultiVector<int,Scalar>.
    */
    static Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > CloneCopy( const Tpetra::MultiVector<int,Scalar>& mv, const std::vector<int>& index )
    { 
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(index.size() == 0,std::runtime_error,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::Clone(mv,index): numvecs must be greater than zero.");
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::runtime_error,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::Clone(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.numVectors(), std::runtime_error,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::Clone(mv,index): indices must be < mv.numVectors().");
#endif
      // TODO: check to see if index is contiguous, if so, use Teuchos::Range1D
      return mv.subCopy(Teuchos::arrayViewFromVector(index));
    }

    /*! \brief Creates a new Tpetra::MultiVector<int,Scalar> that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new Tpetra::MultiVector<int,Scalar>.
    */      
    static Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > CloneView( Tpetra::MultiVector<int,Scalar>& mv, const std::vector<int>& index )
    { 
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): numvecs must be greater than zero.");
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.numVectors(), std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be < mv.numVectors().");
#endif
      return mv.subView(Teuchos::arrayViewFromVector(index));
    }

    /*! \brief Creates a new const Tpetra::MultiVector<int,Scalar> that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new const Tpetra::MultiVector<int,Scalar>.
    */      
    static Teuchos::RCP<const Tpetra::MultiVector<int,Scalar> > CloneView( const Tpetra::MultiVector<int,Scalar>& mv, const std::vector<int>& index )
    { 
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): numvecs must be greater than zero.");
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.numVectors(), std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be < mv.numVectors().");
#endif
      return mv.subViewConst(Teuchos::arrayViewFromVector(index));
    }

    //@}

    //! @name Attribute methods
    //@{ 

    //! Obtain the vector length of \c mv.
    static int GetVecLength( const Tpetra::MultiVector<int,Scalar>& mv )
    { return mv.globalLength(); }

    //! Obtain the number of vectors in \c mv
    static int GetNumberVecs( const Tpetra::MultiVector<int,Scalar>& mv )
    { return mv.numVectors(); }
    //@}

    //! @name Update methods
    //@{ 

    /*! \brief Update \c mv with \f$ \alpha AB + \beta mv \f$.
     */
    static void MvTimesMatAddMv( Scalar alpha, const Tpetra::MultiVector<int,Scalar>& A, 
                                 const Teuchos::SerialDenseMatrix<int,Scalar>& B, 
                                 Scalar beta, Tpetra::MultiVector<int,Scalar>& mv )
    { 
      Tpetra::Map<int> LocalMap(B.numRows(), 0, *A.getMap().getPlatform(),true);
      // TODO: this multivector should be a view of the data of B, not a copy
      Teuchos::ArrayView<const Scalar> Bvalues(B.values(),B.stride()*B.numCols());
      Tpetra::MultiVector<int,Scalar> B_mv(LocalMap,Bvalues,B.stride(),B.numCols());
      mv.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha, A, B_mv, beta);
    }

    /*! \brief Replace \c mv with \f$\alpha A + \beta B\f$.
     */
    static void MvAddMv( Scalar alpha, const Tpetra::MultiVector<int,Scalar>& A, Scalar beta, const Tpetra::MultiVector<int,Scalar>& B, Tpetra::MultiVector<int,Scalar>& mv )
    { 
      mv.update(alpha,A,beta,B,Teuchos::ScalarTraits<Scalar>::zero());
    }

    /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$ \alpha A^Tmv \f$.
    */
    static void MvTransMv( Scalar alpha, const Tpetra::MultiVector<int,Scalar>& A, const Tpetra::MultiVector<int,Scalar>& mv, Teuchos::SerialDenseMatrix<int,Scalar>& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
                          , ConjType conj = Anasazi::CONJ
#endif
                        )
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
    
    /*! \brief Compute a vector \c b where the components are the individual dot-products of the \c i-th columns of \c A and \c mv, i.e.\f$b[i] = A[i]^Tmv[i]\f$.
     */
    static void MvDot( const Tpetra::MultiVector<int,Scalar>& A, const Tpetra::MultiVector<int,Scalar>& B, std::vector<Scalar> &dots
#ifdef HAVE_ANASAZI_EXPERIMENTAL
                      , ConjType conj = Anasazi::CONJ
#endif
                      )
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(A.numVectors() != B.numVectors(),std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvDot(A,B,dots): A and B must have the same number of vectors.");
      TEST_FOR_EXCEPTION(dots.size() != (unsigned int)A.numVectors(),std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvDot(A,B,dots): dots must have room for all dot products.");
#endif
      A.dot(B,Teuchos::arrayViewFromVector(dots));
    }

    //@}
    //! @name Norm method
    //@{ 

    /*! \brief Compute the 2-norm of each individual vector of \c mv.  
      Upon return, \c normvec[i] holds the value of \f$||mv_i||_2\f$, the \c i-th column of \c mv.
    */
    static void MvNorm(const Tpetra::MultiVector<int,Scalar>& mv, std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms)
    { 
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION((unsigned int)mv.numVectors() != norms.size(),std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Tpetra::MultiVector>::MvNorm(mv,norms): norms must be the same size as mv.");
#endif
      mv.norm2(Teuchos::arrayViewFromVector(norms));
    }

    //@}
    
    //! @name Initialization methods
    //@{ 
    /*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.
     */
    static void SetBlock( const Tpetra::MultiVector<int,Scalar>& A, const std::vector<int>& index, Tpetra::MultiVector<int,Scalar>& mv )
    { 
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION((unsigned int)A.numVectors() != index.size(),std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Tpetra::MultiVector>::SetBlock(A,index,mv): index must be the same size as A.");
#endif
      Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > mvsub = mv.subView(Teuchos::arrayViewFromVector(index));
      *mvsub = A;
      mvsub = Teuchos::null;
    }

    /*! \brief Scale each element of the vectors in \c mv with \c alpha.
     */
    static void MvScale ( Tpetra::MultiVector<int,Scalar>& mv, double alpha ) 
    { mv.scale(alpha); }

    /*! \brief Scale each element of the \c i-th vector in \c mv with \c alpha[i].
     */
    static void MvScale ( Tpetra::MultiVector<int,Scalar>& mv, const std::vector<double>& alpha )
    { 
      TEST_FOR_EXCEPT(true);
    }

    /*! \brief Replace the vectors in \c mv with random vectors.
     */
    static void MvRandom( Tpetra::MultiVector<int,Scalar>& mv )
    { mv.random(); }

    /*! \brief Replace each element of the vectors in \c mv with \c alpha.
     */
    static void MvInit( Tpetra::MultiVector<int,Scalar>& mv, double alpha = Teuchos::ScalarTraits<double>::zero() )
    { mv.putScalar(alpha); }

    //@}

    //! @name Print method
    //@{ 

    /*! \brief Print the \c mv multi-vector to the \c os output stream.
     */
    static void MvPrint( const Tpetra::MultiVector<int,Scalar>& mv, std::ostream& os )
    { 
      mv.print(os);
    }

    //@}
  };        


  ///////////////////////////////////////////////////////////////////////// 
  //
  // Implementation of the Anasazi::OperatorTraits for Thyra::LinearOpBase
  //
  ///////////////////////////////////////////////////////////////////////// 

  /*! 
                         
    \brief Template specialization of Anasazi::OperatorTraits class using the
    Thyra::LinearOpBase virtual base class and Thyra::MultiVectorBase class.

    This interface will ensure that any LinearOpBase and MultiVectorBase
    implementations will be accepted by the Anasazi templated solvers.

  */
  template <class Scalar> 
  class OperatorTraits < Scalar, Tpetra::MultiVector<int,Scalar>, Tpetra::Operator<int,Scalar> >
  {
  public:
    
    /*! \brief This method takes the MultiVectorBase \c x and
      applies the LinearOpBase \c Op to it resulting in the MultiVectorBase \c y.
    */    
    static void Apply ( const Tpetra::Operator<int,Scalar>& Op, const Tpetra::MultiVector<int,Scalar>& x, Tpetra::MultiVector<int,Scalar>& y )
    { 
      Op.apply(x,y);
    }
    
  };
  
  
} // end of Anasazi namespace 

#endif 
// end of file ANASAZI_TPETRA_ADAPTER_HPP
