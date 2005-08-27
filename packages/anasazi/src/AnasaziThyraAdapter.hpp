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

/*! \file AnasaziThyraAdapter.hpp
  \brief Specializations of the Anasazi multi-vector and operator traits 
  classes using Thyra base classes MultiVectorBase and LinearOpBase.
*/

#ifndef ANASAZI_THYRA_ADAPTER_HPP
#define ANASAZI_THYRA_ADAPTER_HPP

#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziConfigDefs.hpp"

#include <Thyra_ExplicitMultiVectorView.hpp>
#include <Thyra_MultiVectorBaseDecl.hpp>

namespace Anasazi {
  
  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::MultiVecTraits for Epetra::MultiVector.
  //
  ////////////////////////////////////////////////////////////////////

  /*! \class MultiVecTraits< ScalarType, Thyra::MultiVectorBase<ScalarType> >
    \brief Template specialization of Anasazi::MultiVecTraits class using the Epetra_MultiVector class.

    This interface will ensure that any Epetra_MultiVector will be accepted by the Anasazi
    templated solvers.  

  */

  template<class ScalarType>
  class MultiVecTraits< ScalarType, Thyra::MultiVectorBase<ScalarType> >
  {
  private:
    typedef Thyra::MultiVectorBase<ScalarType> TMVB;
  public:

    //@{ \name Creation methods

    /*! \brief Creates a new empty MultiVectorBase containing \c numvecs columns.
      
    \return Reference-counted pointer to the new MultiVectorBase.
    */
    static Teuchos::RefCountPtr<TMVB> Clone( const TMVB& mv, const int numvecs )
    { 
      Teuchos::RefCountPtr<TMVB> c = Thyra::createMembers( mv.range(), numvecs ); 
      Thyra::assign(&(*c),Teuchos::ScalarTraits<ScalarType>::zero());
      return c;
    }

    /*! \brief Creates a new MultiVectorBase and copies contents of \c mv into the new vector (deep copy).
      
      \return Reference-counted pointer to the new MultiVectorBase.
    */
    static Teuchos::RefCountPtr<TMVB> CloneCopy( const TMVB& mv )
    { 
      int numvecs = mv.domain()->dim();
      // create the new multivector
      Teuchos::RefCountPtr< TMVB > cc = Thyra::createMembers( mv.range(), numvecs );
      // copy the data from the source multivector to the new multivector
      Thyra::assign(&(*cc), mv);
      return cc;
    }

    /*! \brief Creates a new MultiVectorBase and copies the selected contents of \c mv into the new vector (deep copy).  

      The copied vectors from \c mv are indicated by the \c indeX.size() indices in \c index.      
      \return Reference-counted pointer to the new MultiVectorBase.
    */
    static Teuchos::RefCountPtr<TMVB> CloneCopy( const TMVB& mv, const std::vector<int>& index )
    { 
      int numvecs = index.size();
      std::vector<int> index_plus(index);
      // Anasazi::MultiVecTraits uses [0,numvecs-1] indexing, 
      // while Thyra uses [1,numvecs] indexing
      for (int i=0; i<numvecs; i++) {
        index_plus[i]++;
      }
      // create the new multivector
      Teuchos::RefCountPtr< TMVB > cc = Thyra::createMembers( mv.range(), numvecs );
      // create a view to the relevant part of the source multivector
      Teuchos::RefCountPtr< const TMVB > relview = mv.subView( numvecs, &(index_plus[0]) );
      // copy the data from the relevant view to the new multivector
      Thyra::assign(&(*cc), *relview);
      return cc;
    }

    /*! \brief Creates a new MultiVectorBase that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new MultiVectorBase.
    */      
    static Teuchos::RefCountPtr<TMVB> CloneView( TMVB& mv, const std::vector<int>& index )
    {
      int numvecs = index.size();
      std::vector<int> index_plus(index);
      // Anasazi::MultiVecTraits uses [0,numvecs-1] indexing, 
      // while Thyra uses [1,numvecs] indexing
      for (int i=0; i<numvecs; i++) {
        index_plus[i]++;
      }
      // create a view to the relevant part of the source multivector
      Teuchos::RefCountPtr< TMVB > cc = mv.subView( numvecs, &(index_plus[0]) );
      return cc;
    }

    /*! \brief Creates a new const MultiVectorBase that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new const MultiVectorBase.
    */      
    static Teuchos::RefCountPtr<const TMVB> CloneView( const TMVB& mv, const std::vector<int>& index )
    {
      int numvecs = index.size();
      std::vector<int> index_plus(index);
      // Anasazi::MultiVecTraits uses [0,numvecs-1] indexing, 
      // while Thyra uses [1,numvecs] indexing
      for (int i=0; i<numvecs; i++) {
        index_plus[i]++;
      }
      // create a view to the relevant part of the source multivector
      Teuchos::RefCountPtr< const TMVB > cc = mv.subView( numvecs, &(index_plus[0]) );
      return cc;
    }

    //@}

    //@{ \name Attribute methods

    //! Obtain the vector length of \c mv.
    static int GetVecLength( const TMVB& mv )
    { return mv.range()->dim(); }

    //! Obtain the number of vectors in \c mv
    static int GetNumberVecs( const TMVB& mv )
    { return mv.domain()->dim(); }
    //@}

    //@{ \name Update methods

    /*! \brief Update \c mv with \f$ \alpha AB + \beta mv \f$.
     */
    static void MvTimesMatAddMv( const ScalarType alpha, const TMVB& A, 
				 const Teuchos::SerialDenseMatrix<int,ScalarType>& B, 
				 const ScalarType beta, TMVB& mv )
    {
      typedef Thyra::VectorBase<ScalarType> TVB;
      // create pointers to the columns of A
      int i, j;
      int m = A.domain()->dim();
      int n = mv.domain()->dim();
      std::vector<const TVB*> acols(m);
      for (i=0; i<m; i++) {
        // col() returns a RefCountPtr. Deref it and then grab the address.
        // remember, col is indexed like [1,m] (not [0,m-1])
        acols[i] = &(*A.col(i+1));
      }

      // create storage for the columns of B
      // we must have separate storage because we have to have somewhere
      // to store the product of alpha*B
      //
      // Note: we could just make a copy of B and then perform alpha*B. this
      // would use more memory (probably not an issue), while providing
      // negligible increases in efficiency. 
      std::vector<ScalarType> btemp(m);

      // now just perform the matrix multiplication, for each column of 
      // the answer.
      for (i=0; i<n; i++) {
        for (j=0; j<m; j++) {
          btemp[j] = alpha*B(j,i);
        }
        // remember, col is indexed like [1,n] (not [0,n-1])
        Thyra::linear_combination(m,&(btemp[0]),&(acols[0]),beta,&(*mv.col(i+1)) );
      }
    }

    /*! \brief Replace \c mv with \f$\alpha A + \beta B\f$.
     */
    static void MvAddMv( const ScalarType alpha, const TMVB& A, 
                         const ScalarType beta,  const TMVB& B, TMVB& mv )
    { 
      // Set mv = 0
      Thyra::assign(&mv,Teuchos::ScalarTraits<ScalarType>::zero());
      // Add alpha*A to mv
      Thyra::update(alpha,A,&mv);
      // Add beta*B to mv
      Thyra::update(beta ,B,&mv);
    }

    /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$ \alpha A^Tmv \f$.
    */
    static void MvTransMv( const ScalarType alpha, const TMVB& A, const TMVB& mv, 
                           Teuchos::SerialDenseMatrix<int,ScalarType>& B )
    { 
      // Create a multivector to hold the result (m by n)
      int m = A.domain()->dim();
      int n = mv.domain()->dim();
      Teuchos::RefCountPtr< TMVB > temp = Thyra::createMembers( A.domain(), n ); 
      // Note: using the conjugate transpose as opposed to a simple transpose
      A.applyTranspose(Thyra::CONJ_ELE,mv,&(*temp),alpha,Teuchos::ScalarTraits<ScalarType>::zero());
      // Move the data from the temporary multivector into the SerialDenseMatrix
      // To get the data out, we need a MutableMultiVectorView
      Thyra::ExplicitMutableMultiVectorView<ScalarType> cc(*temp);
      ScalarType *vals = cc.values();
      for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
          B(i,j) = vals[i+j*m];
        }
      }
    }
    
    /*! \brief Compute a vector \c b where the components are the individual dot-products of the \c i-th columns of \c A and \c mv, i.e.\f$b[i] = A[i]^Tmv[i]\f$.
     */
    static void MvDot( const TMVB& mv, const TMVB& A, std::vector<ScalarType>* b )
    { Thyra::dots(mv,A,&((*b)[0])); }

    //@}
    //@{ \name Norm method

    /*! \brief Compute the 2-norm of each individual vector of \c mv.  
      Upon return, \c normvec[i] holds the value of \f$||mv_i||_2\f$, the \c i-th column of \c mv.
    */
    static void MvNorm( const TMVB& mv, std::vector<ScalarType>* normvec )
    { Thyra::norms_2(mv,&((*normvec)[0])); }

    //@}

    //@{ \name Initialization methods
    /*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.
     */
    static void SetBlock( const TMVB& A, const std::vector<int>& index, TMVB& mv )
    { 
      // Extract the "numvecs" columns of mv indicated by the index vector.
      int numvecs = index.size();
      std::vector<int> index_plus(index), indexA(numvecs);
      int numAcols = A.domain()->dim();
      // Anasazi::MultiVecTraits uses [0,numvecs-1] indexing, 
      // while Thyra uses [1,numvecs] indexing
      for (int i=0; i<numvecs; i++) {
        index_plus[i]++;
        indexA[i] = i+1;
      }
      // Thyra::assign requires that both arguments have the same number of 
      // vectors. Enforce this, by shrinking one to match the other.
      if ( numAcols < numvecs ) {
        // A does not have enough columns to satisfy index_plus. Shrink
        // index_plus.
        numvecs = numAcols;
        index_plus.resize( numvecs );
      }
      else if ( numAcols > numvecs ) {
        numAcols = numvecs;
        indexA.resize( numAcols );
      }
      // create a view to the relevant part of the source multivector
      Teuchos::RefCountPtr< const TMVB > relsource = A.subView( numAcols, &(indexA[0]) );
      // create a view to the relevant part of the destination multivector
      Teuchos::RefCountPtr< TMVB > reldest = mv.subView( numvecs, &(index_plus[0]) );
      // copy the data to the destination multivector subview
      Thyra::assign(&(*reldest), *relsource);
    }

    /*! \brief Replace the vectors in \c mv with random vectors.
     */
    static void MvRandom( TMVB& mv )
    { 
      // Thyra::randomize generates via a uniform distribution on [l,u]
      // We will use this to generate on [-1,1]
      Thyra::randomize(-Teuchos::ScalarTraits<ScalarType>::one(),
                        Teuchos::ScalarTraits<ScalarType>::one(),
                       &mv);
    }

    /*! \brief Replace each element of the vectors in \c mv with \c alpha.
     */
    static void MvInit( TMVB& mv, ScalarType alpha = Teuchos::ScalarTraits<ScalarType>::zero() )
    { Thyra::assign(&mv,alpha); }

    //@}

    //@{ \name Print method

    /*! \brief Print the \c mv multi-vector to the \c os output stream.
     */
    static void MvPrint( const TMVB& mv, ostream& os )
    { }

    //@}
  };        

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::OperatorTraits for Epetra::Operator.
  //
  ////////////////////////////////////////////////////////////////////

  /*! \class OperatorTraits< ScalarType ,
                             Thyra::MultiVectorBase<ScalarType> ,  
                             Thyra::LinearOpBase<ScalarType> >
    \brief Template specialization of Anasazi::OperatorTraits class using the Epetra_Operator virtual base class and 
    Epetra_MultiVector class.

    This interface will ensure that any Epetra_Operator and Epetra_MultiVector will be accepted by the Anasazi
    templated solvers.

  */

  template <class ScalarType> 
  class OperatorTraits < ScalarType, Thyra::MultiVectorBase<ScalarType>, Thyra::LinearOpBase<ScalarType> >
  {
  private:
    typedef Thyra::MultiVectorBase<ScalarType> TMVB;
    typedef Thyra::LinearOpBase<ScalarType>    TLOB;
  public:
    
    /*! \brief This method takes the Epetra_MultiVector \c x and
      applies the Epetra_Operator \c Op to it resulting in the Epetra_MultiVector \c y.
    */    
    static ReturnType Apply ( const TLOB& Op, const TMVB& x, TMVB& y )
    { 
      Op.apply(Thyra::NONCONJ_ELE,x,&y);
      return Ok; 
    }
    
  };
  
} // end of Anasazi namespace 

#endif 
// end of file ANASAZI_EPETRA_ADAPTER_HPP
