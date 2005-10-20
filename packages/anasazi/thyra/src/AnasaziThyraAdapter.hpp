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
  classes using Thyra base classes LinearOpBase and MultiVectorBase.
*/

#ifndef ANASAZI_THYRA_ADAPTER_HPP
#define ANASAZI_THYRA_ADAPTER_HPP

#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziConfigDefs.hpp"

#include <Thyra_ExplicitMultiVectorView.hpp>
#include <Thyra_MultiVectorBaseDecl.hpp>

namespace Anasazi {
  
  ////////////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::MultiVecTraits for Thyra::MultiVectorBase
  //
  ////////////////////////////////////////////////////////////////////////////

  /*! \class MultiVecTraits< ScalarType, Thyra::MultiVectorBase<ScalarType> >
    \brief Template specialization of Anasazi::MultiVecTraits class using the
    Thyra::MultiVectorBase class.

    This interface will ensure that any implementation of MultiVectorBaseClass 
    will be accepted by the Anasazi templated solvers.  

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

      // We do not assume that the indices are sorted, nor do we check that
      // index.size() > 0. This code is fail-safe, in the sense that a zero
      // length index vector will pass the error on the Thyra.

      // Thyra has two ways to create an indexed View:
      // * contiguous (via a range of columns)
      // * indexed (via a vector of column indices)
      // The former is significantly more efficient than the latter, in terms of
      // computations performed with/against the created view.
      // We will therefore check to see if the given indices are contiguous, and
      // if so, we will use the contiguous view creation method.
      
      // Anasazi::MultiVecTraits uses 0-indexing, while Thyra uses 1-indexing. 
      // Increment the indices as we check for contiguity.
      int lb = index_plus[0]+1;
      bool contig = true;
      for (int i=0; i<numvecs; i++) {
        index_plus[i]++;
        if (lb+i != index_plus[i]) contig = false;
      }

      Teuchos::RefCountPtr< TMVB > cc;
      if (contig) {
        RangePack::Range1D rng(lb,lb+numvecs-1);
        // create a contiguous view to the relevant part of the source multivector
        cc = mv.subView(rng);
      }
      else {
        // create an indexed view to the relevant part of the source multivector
        cc = mv.subView( numvecs, &(index_plus[0]) );
      }
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

      // We do not assume that the indices are sorted, nor do we check that
      // index.size() > 0. This code is fail-safe, in the sense that a zero
      // length index vector will pass the error on the Thyra.

      // Thyra has two ways to create an indexed View:
      // * contiguous (via a range of columns)
      // * indexed (via a vector of column indices)
      // The former is significantly more efficient than the latter, in terms of
      // computations performed with/against the created view.
      // We will therefore check to see if the given indices are contiguous, and
      // if so, we will use the contiguous view creation method.
      
      // Anasazi::MultiVecTraits uses 0-indexing, while Thyra uses 1-indexing. 
      // Increment the indices as we check for contiguity.
      int lb = index_plus[0]+1;
      bool contig = true;
      for (int i=0; i<numvecs; i++) {
        index_plus[i]++;
        if (lb+i != index_plus[i]) contig = false;
      }

      Teuchos::RefCountPtr< const TMVB > cc;
      if (contig) {
        RangePack::Range1D rng(lb,lb+numvecs-1);
        // create a contiguous view to the relevant part of the source multivector
        cc = mv.subView(rng);
      }
      else {
        // create an indexed view to the relevant part of the source multivector
        cc = mv.subView( numvecs, &(index_plus[0]) );
      }
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
      int m = B.numRows();
      int n = B.numCols();
      // create a MultiVectorBase object to represent the matrix
      Teuchos::RefCountPtr< TMVB > mat = Thyra::createMembers( A.domain(), n );
      // get an explicit view which allows us to access the data
      Thyra::ExplicitMutableMultiVectorView<ScalarType> mat_view(*mat);
      // copy the data
      for (int j=0; j<n; j++) {
        for (int i=0; i<m; i++) {
          mat_view(i+1,j+1) = B(i,j);
        }
      }
      // perform the operation via A: mv <- alpha*A*mat + beta*mv
      A.apply(Thyra::NONCONJ_ELE,*mat,&mv,alpha,beta);
    }

    /*! \brief Replace \c mv with \f$\alpha A + \beta B\f$.
     */
    static void MvAddMv( const ScalarType alpha, const TMVB& A, 
                         const ScalarType beta,  const TMVB& B, TMVB& mv )
    { 
      ScalarType coef[2], zero = Teuchos::ScalarTraits<ScalarType>::zero();
      const TMVB * in[2];

      in[0] = &A;
      in[1] = &B;
      coef[0] = alpha;
      coef[1] = beta;

      Thyra::linear_combination(2,coef,in,zero,&mv);
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

    /*! \brief Compute a vector \c b where the components are the individual dot-products of the 
        \c i-th columns of \c A and \c mv, i.e.\f$b[i] = A[i]^Tmv[i]\f$.
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

  ///////////////////////////////////////////////////////////////////////// 
  //
  // Implementation of the Anasazi::OperatorTraits for Thyra::LinearOpBase
  //
  ///////////////////////////////////////////////////////////////////////// 

  /*! \class OperatorTraits< ScalarType ,
                             Thyra::MultiVectorBase<ScalarType> ,  
                             Thyra::LinearOpBase<ScalarType> >
    \brief Template specialization of Anasazi::OperatorTraits class using the
    Thyra::LinearOpBase virtual base class and Thyra::MultiVectorBase class.

    This interface will ensure that any LinearOpBase and MultiVectorBase
    implementations will be accepted by the Anasazi templated solvers.

  */

  template <class ScalarType> 
  class OperatorTraits < ScalarType, Thyra::MultiVectorBase<ScalarType>, Thyra::LinearOpBase<ScalarType> >
  {
  private:
    typedef Thyra::MultiVectorBase<ScalarType> TMVB;
    typedef Thyra::LinearOpBase<ScalarType>    TLOB;
  public:
    
    /*! \brief This method takes the MultiVectorBase \c x and
      applies the LinearOpBase \c Op to it resulting in the MultiVectorBase \c y.
    */    
    static ReturnType Apply ( const TLOB& Op, const TMVB& x, TMVB& y )
    { 
      Op.apply(Thyra::NONCONJ_ELE,x,&y);
      return Ok; 
    }
    
  };
  
} // end of Anasazi namespace 

#endif 
// end of file ANASAZI_THYRA_ADAPTER_HPP
