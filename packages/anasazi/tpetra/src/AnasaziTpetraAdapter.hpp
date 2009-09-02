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

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"


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

  template<class Scalar, class LO, class GO>
  class MultiVecTraits<Scalar, Tpetra::MultiVector<Scalar,LO,GO> >
  {
  public:

    //! @name Creation methods
    //@{ 

    /*! \brief Creates a new empty Tpetra::MultiVector<Scalar,LO,GO> containing \c numvecs columns.
      
    \return Reference-counted pointer to the new Tpetra::MultiVector<Scalar,LO,GO>.
    */
    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO> > Clone( const Tpetra::MultiVector<Scalar,LO,GO>& mv, const int numvecs )
    { 
#ifdef TPETRA_DEBUG
      TEST_FOR_EXCEPTION(numvecs <= 0,std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::Clone(mv,numvecs): numvecs must be greater than zero.");
#endif
      return Teuchos::rcp( new Tpetra::MultiVector<Scalar,LO,GO>(mv.getMap(), numvecs) ); 
    }

    /*! \brief Creates a new Tpetra::MultiVector<Scalar,LO,GO> and copies contents of \c mv into the new vector (deep copy).
      
      \return Reference-counted pointer to the new Tpetra::MultiVector<Scalar,LO,GO>.
    */
    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO> > CloneCopy( const Tpetra::MultiVector<Scalar,LO,GO>& mv )
    {
      return Teuchos::rcp( new Tpetra::MultiVector<Scalar,LO,GO>( mv ) ); 
    }

    /*! \brief Creates a new Tpetra::MultiVector<Scalar,LO,GO> and copies the selected contents of \c mv into the new vector (deep copy).  

      The copied vectors from \c mv are indicated by the \c indeX.size() indices in \c index.      
      \return Reference-counted pointer to the new Tpetra::MultiVector<Scalar,LO,GO>.
    */
    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO> > CloneCopy( const Tpetra::MultiVector<Scalar,LO,GO>& mv, const std::vector<int>& index )
    { 
      TEST_FOR_EXCEPTION(index.size() == 0,std::runtime_error,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::Clone(mv,index): numvecs must be greater than zero.");
#ifdef TPETRA_DEBUG
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::runtime_error,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::Clone(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.numVectors(), std::runtime_error,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::Clone(mv,index): indices must be < mv.numVectors().");
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

    /*! \brief Creates a new Tpetra::MultiVector<Scalar,LO,GO> that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new Tpetra::MultiVector<Scalar,LO,GO>.
    */      
    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO> > CloneView( Tpetra::MultiVector<Scalar,LO,GO>& mv, const std::vector<int>& index )
    {
      TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): numvecs must be greater than zero.");
#ifdef TPETRA_DEBUG
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.numVectors(), std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be < mv.numVectors().");
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

    /*! \brief Creates a new const Tpetra::MultiVector<Scalar,LO,GO> that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new const Tpetra::MultiVector<Scalar,LO,GO>.
    */      
    static Teuchos::RCP<const Tpetra::MultiVector<Scalar,LO,GO> > CloneView( const Tpetra::MultiVector<Scalar,LO,GO>& mv, const std::vector<int>& index )
    {
      TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): numvecs must be greater than zero.");
#ifdef TPETRA_DEBUG
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.numVectors(), std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index): indices must be < mv.numVectors().");
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

    //@}

    //! @name Attribute methods
    //@{ 

    //! Obtain the vector length of \c mv.
    static int GetVecLength( const Tpetra::MultiVector<Scalar,LO,GO>& mv )
    { return mv.globalLength(); }

    //! Obtain the number of vectors in \c mv
    static int GetNumberVecs( const Tpetra::MultiVector<Scalar,LO,GO>& mv )
    { return mv.numVectors(); }
    //@}

    //! @name Update methods
    //@{ 

    /*! \brief Update \c mv with \f$ \alpha AB + \beta mv \f$.
     */
    static void MvTimesMatAddMv(       Scalar alpha, const Tpetra::MultiVector<Scalar,LO,GO>& A, 
                                 const Teuchos::SerialDenseMatrix<int,Scalar>& B, 
                                 Scalar beta, Tpetra::MultiVector<Scalar,LO,GO>& mv )
    {
      Tpetra::Map<LO,GO> LocalMap(B.numRows(), 0, A.getMap().getComm(),true);
      // TODO: this multivector should be a view of the data of B, not a copy
      Teuchos::ArrayView<const Scalar> Bvalues(B.values(),B.stride()*B.numCols());
      Tpetra::MultiVector<Scalar,LO,GO> B_mv(LocalMap,Bvalues,B.stride(),B.numCols());
      mv.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha, A, B_mv, beta);
    }

    /*! \brief Replace \c mv with \f$\alpha A + \beta B\f$.
     */
    static void MvAddMv( Scalar alpha, const Tpetra::MultiVector<Scalar,LO,GO>& A, Scalar beta, const Tpetra::MultiVector<Scalar,LO,GO>& B, Tpetra::MultiVector<Scalar,LO,GO>& mv )
    { 
      mv.update(alpha,A,beta,B,Teuchos::ScalarTraits<Scalar>::zero());
    }

    /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$ \alpha A^Tmv \f$.
    */
    static void MvTransMv( Scalar alpha, const Tpetra::MultiVector<Scalar,LO,GO>& A, const Tpetra::MultiVector<Scalar,LO,GO>& mv, Teuchos::SerialDenseMatrix<int,Scalar>& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
                          , ConjType conj = Anasazi::CONJ
#endif
                        )
    { 
      Tpetra::Map<LO,GO> LocalMap(B.numRows(), 0, A.getMap().getComm(),true);
      // TODO: this multivector should be a view of the data of B, so we don't have to perform the copy afterwards
      Tpetra::MultiVector<Scalar,LO,GO> B_mv(LocalMap,B.numCols(),true);
      B_mv.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,alpha,A,mv,Teuchos::ScalarTraits<Scalar>::zero());
      Teuchos::ArrayView<Scalar> av(B.values(),B.stride()*B.numCols());
      B_mv.extractCopy1D(av,B.stride());
    }

    /*! \brief Compute a vector \c b where the components are the individual dot-products of the \c i-th columns of \c A and \c mv, i.e.\f$b[i] = A[i]^Tmv[i]\f$.
     */
    static void MvDot( const Tpetra::MultiVector<Scalar,LO,GO>& A, const Tpetra::MultiVector<Scalar,LO,GO>& B, std::vector<Scalar> &dots
#ifdef HAVE_ANASAZI_EXPERIMENTAL
                      , ConjType conj = Anasazi::CONJ
#endif
                      )
    {
#ifdef TPETRA_DEBUG
      TEST_FOR_EXCEPTION(dots.size() < (typename std::vector<int>::size_type)A.numVectors(),std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvDot(A,B,dots): dots must have room for all dot products.");
#endif
      Teuchos::ArrayView<Scalar> av(dots);
      A.dot(B,av(0,A.numVectors()));
    }

    //@}
    //! @name Norm method
    //@{ 

    /*! \brief Compute the 2-norm of each individual vector of \c mv.  
      Upon return, \c normvec[i] holds the value of \f$||mv_i||_2\f$, the \c i-th column of \c mv.
    */
    static void MvNorm(const Tpetra::MultiVector<Scalar,LO,GO>& mv, std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec)
    { 
#ifdef TPETRA_DEBUG
      TEST_FOR_EXCEPTION(normvec.size() < (typename std::vector<int>::size_type)mv.numVectors(),std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvNorm(mv,normvec): normvec must have room for all norms.");
#endif
      Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> av(normvec);
      mv.norm2(av(0,mv.numVectors()));
    }

    //@}
    
    //! @name Initialization methods
    //@{ 
    /*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.
     */
    static void SetBlock( const Tpetra::MultiVector<Scalar,LO,GO>& A, const std::vector<int>& index, Tpetra::MultiVector<Scalar,LO,GO>& mv )
    {
#ifdef TPETRA_DEBUG
      TEST_FOR_EXCEPTION((typename std::vector<int>::size_type)A.numVectors() < index.size(),std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::SetBlock(A,index,mv): index must be the same size as A.");
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

    /*! \brief Scale each element of the vectors in \c mv with \c alpha.
     */
    static void MvScale ( Tpetra::MultiVector<Scalar,LO,GO>& mv, Scalar alpha ) 
    { mv.scale(alpha); }

    /*! \brief Scale each element of the \c i-th vector in \c mv with \c alpha[i].
     */
    static void MvScale ( Tpetra::MultiVector<Scalar,LO,GO>& mv, const std::vector<Scalar>& alpha )
    { 
      mv.scale(alpha);
    }

    /*! \brief Replace the vectors in \c mv with random vectors.
     */
    static void MvRandom( Tpetra::MultiVector<Scalar,LO,GO>& mv )
    { mv.random(); }

    /*! \brief Replace each element of the vectors in \c mv with \c alpha.
     */
    static void MvInit( Tpetra::MultiVector<Scalar,LO,GO>& mv, Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() )
    { mv.putScalar(alpha); }

    //@}

    //! @name Print method
    //@{ 

    /*! \brief Print the \c mv multi-vector to the \c os output stream.
     */
    static void MvPrint( const Tpetra::MultiVector<Scalar,LO,GO>& mv, std::ostream& os )
    { mv.print(os); }

    //@}
  };        

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::OperatorTraits for Tpetra::Operator.
  //
  ////////////////////////////////////////////////////////////////////

  template <class Scalar, class LO, class GO> 
  class OperatorTraits < Scalar, Tpetra::MultiVector<Scalar,LO,GO>, Tpetra::Operator<Scalar,LO,GO> >
  {
  public:
    static void Apply ( const Tpetra::Operator<Scalar,LO,GO> & Op, 
                        const Tpetra::MultiVector<Scalar,LO,GO> & X,
                              Tpetra::MultiVector<Scalar,LO,GO> & Y)
    { 
      Op.apply(X,Y);
    }
  };

} // end of Anasazi namespace 

#endif 
// end of file ANASAZI_TPETRA_ADAPTER_HPP
