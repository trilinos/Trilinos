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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

/*! \file AnasaziSaddleContainer.hpp
 *  \brief Stores a set of vectors of the form [V; L] where V is a distributed multivector
 *  and L is a serialdense matrix.
 *
 *  Used in solving TraceMin's saddle point problem.
*/

#ifndef ANASAZI_SADDLE_CONTAINER_HPP
#define ANASAZI_SADDLE_CONTAINER_HPP

#include "AnasaziConfigDefs.hpp"
#include "Teuchos_VerboseObject.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace Anasazi {
namespace Experimental {

template <class ScalarType, class MV>
class SaddleContainer //: public Anasazi::SaddleContainer<ScalarType, MV>
{
  template <class Scalar_, class MV_, class OP_> friend class SaddleOperator;

private:
  typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVT;
  typedef Teuchos::SerialDenseMatrix<int,ScalarType> SerialDenseMatrix;
  const ScalarType ONE; 
  RCP<MV> X_;
  RCP<SerialDenseMatrix> Y_;

public:
  // Constructors
  SaddleContainer( ) : ONE(Teuchos::ScalarTraits<ScalarType>::one()) { };
  SaddleContainer( const RCP<MV> X, bool eye=false );

  // Things that are necessary for compilation
  // Returns a clone of the current vector
  SaddleContainer<ScalarType, MV> * Clone(const int nvecs) const;
  // Returns a duplicate of the current vector
  SaddleContainer<ScalarType, MV> * CloneCopy() const;
  // Returns a duplicate of the specified vectors
  SaddleContainer<ScalarType, MV> * CloneCopy(const std::vector< int > &index) const;
  int GetVecLength() const;
  int GetNumberVecs() const { return MVT::GetNumberVecs(*X_); };
  // Update *this with alpha * A * B + beta * (*this)
  void MvTimesMatAddMv(ScalarType alpha, const SaddleContainer<ScalarType,MV> &A, 
                       const Teuchos::SerialDenseMatrix<int, ScalarType> &B, 
                       ScalarType beta);
  // Replace *this with alpha * A + beta * B
  void MvAddMv(ScalarType alpha, const SaddleContainer<ScalarType,MV>& A, 
               ScalarType beta,  const SaddleContainer<ScalarType,MV>& B);
  // Scale the vectors by alpha
  void MvScale( ScalarType alpha );
  // Scale the i-th vector by alpha[i]
  void MvScale( const std::vector<ScalarType>& alpha );
  // Compute a dense matrix B through the matrix-matrix multiply alpha * A^H * (*this)
  void MvTransMv (ScalarType alpha, const SaddleContainer<ScalarType, MV>& A, 
                  Teuchos::SerialDenseMatrix< int, ScalarType >& B) const;
  // Compute a vector b where the components are the individual dot-products, i.e.b[i] = A[i]^H*this[i] where A[i] is the i-th column of A. 
  void MvDot (const SaddleContainer<ScalarType, MV>& A, std::vector<ScalarType> &b) const;
  // Compute the 2-norm of each individual vector
  void MvNorm ( std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &normvec) const;
  // Copy the vectors in A to a set of vectors in *this. The numvecs vectors in 
  // A are copied to a subset of vectors in *this indicated by the indices given 
  // in index.
  void SetBlock (const SaddleContainer<ScalarType, MV>& A, const std::vector<int> &index);
  // Deep copy.
  void Assign (const SaddleContainer<ScalarType, MV>&A);
  // Fill the vectors in *this with random numbers.
  void  MvRandom ();
  // Replace each element of the vectors in *this with alpha.
  void  MvInit (ScalarType alpha);
  // Prints the multivector to an output stream
  void MvPrint (std::ostream &os) const;
};



// Constructor
template <class ScalarType, class MV>
SaddleContainer<ScalarType, MV>::SaddleContainer( const RCP<MV> X, bool eye )
: ONE(Teuchos::ScalarTraits<ScalarType>::one()) 
{
  int nvecs = MVT::GetNumberVecs(*X);

  if(eye)
  {
    // Initialize X_ as all 0s
    X_ = MVT::Clone(*X, nvecs);
    MVT::MvInit(*X_);
    
    // Initialize Y to be I
    Y_ = rcp(new SerialDenseMatrix(nvecs,nvecs));
    for(int i=0; i < nvecs; i++)
      (*Y_)(i,i) = ONE;
  }
  else
  {
    // Point X_ to X
    X_ = X;

    // Initialize Y to be 0
    Y_ = rcp(new SerialDenseMatrix(nvecs,nvecs));
  }
}



// Returns a clone of the current vector
template <class ScalarType, class MV>
SaddleContainer<ScalarType, MV> * SaddleContainer<ScalarType, MV>::Clone(const int nvecs) const 
{  
  SaddleContainer<ScalarType, MV> * newSC = new SaddleContainer<ScalarType, MV>();
  
  newSC->X_ = MVT::Clone(*X_,nvecs);
  newSC->Y_ = rcp(new SerialDenseMatrix(Y_->numRows(),Y_->numCols()));

  return newSC;
}



// Returns a duplicate of the current vector
template <class ScalarType, class MV>
SaddleContainer<ScalarType, MV> * SaddleContainer<ScalarType, MV>::CloneCopy() const 
{
  SaddleContainer<ScalarType, MV> * newSC = new SaddleContainer<ScalarType, MV>();
  
  newSC->X_ = MVT::CloneCopy(*X_);
  newSC->Y_ = rcp(new SerialDenseMatrix(*Y_));
  
  return newSC;
}



// Returns a duplicate of the specified vectors
template <class ScalarType, class MV>
SaddleContainer<ScalarType, MV> * SaddleContainer<ScalarType, MV>::CloneCopy(const std::vector< int > &index) const
{
  SaddleContainer<ScalarType, MV> * newSC = new SaddleContainer<ScalarType, MV>();
  
  newSC->X_ = MVT::CloneCopy(*X_,index);
  
  int ncols = index.size();
  int nrows = Y_->numRows();
  newSC->Y_ = rcp(new SerialDenseMatrix(nrows,ncols));
  for(int c=0; c < ncols; c++)
  {
    for(int r=0; r < nrows; r++)
      (*newSC->Y_)(r,c) = (*Y_)(r,index[c]);
  }
  
  return newSC;
}



// Replace *this with alpha * A + beta * B
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvAddMv(ScalarType alpha, const SaddleContainer<ScalarType,MV>& A, 
                                              ScalarType beta,  const SaddleContainer<ScalarType,MV>& B)
{
  MVT::MvAddMv(alpha, *(A.X_), beta, *(B.X_), *X_);

  int ncolsA = A.Y_->numCols();
  int ncolsThis = Y_->numCols();
  int nrows = Y_->numRows();
    
  // check whether space needs to be reallocated
  if(ncolsA != ncolsThis)
    Y_->shapeUninitialized(nrows,ncolsA);

  // Y = alpha A
  Y_->assign(*A.Y_);
  if(alpha != ONE)
    Y_->scale(alpha);
  // Y += beta B
  if(beta == ONE)
    *Y_ += *B.Y_;
  else if(beta == -ONE)
    *Y_ -= *B.Y_;
  else
  {
    SerialDenseMatrix scaledB(*B.Y_);
    scaledB.scale(beta);
    *Y_ += *B.Y_;
  }
}



// Scale the vectors by alpha
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvScale( ScalarType alpha )
{
  MVT::MvScale(*X_, alpha);
  Y_->scale(alpha);
}



// Scale the i-th vector by alpha[i]
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvScale( const std::vector<ScalarType>& alpha )
{
  MVT::MvScale(*X_, alpha);

  int nrows = Y_->numRows();
  int ncols = Y_->numCols();
  
  for(int c=0; c<ncols; c++)
  {
    for(int r=0; r<nrows; r++)
      (*Y_)(r,c) *= alpha[c];
  }
}



// Compute a vector b where the components are the individual dot-products, i.e.b[i] = A[i]^H*this[i] where A[i] is the i-th column of A. 
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvDot (const SaddleContainer<ScalarType, MV>& A, std::vector<ScalarType> &b) const
{
  MVT::MvDot(*X_, *(A.X_), b);

  int nrows = Y_->numRows();
  int ncols = Y_->numCols();
  
  for(int c=0; c < ncols; c++)
  {
    for(int r=0; r < nrows; r++)
    {
      b[c] += ((*A.Y_)(r,c) * (*Y_)(r,c));
    }
  }
}



// Copy the vectors in A to a set of vectors in *this. The numvecs vectors in 
// A are copied to a subset of vectors in *this indicated by the indices given 
// in index.
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::SetBlock (const SaddleContainer<ScalarType, MV>& A, const std::vector<int> &index)
{
  MVT::SetBlock(*(A.X_), index, *X_);

  int nrows = Y_->numRows();
  
  int nvecs = index.size();
  for(int c=0; c<nvecs; c++)
  {
    for(int r=0; r<nrows; r++)
      (*Y_)(r,index[c]) = (*A.Y_)(r,c);
  }
}



// Deep copy.
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::Assign (const SaddleContainer<ScalarType, MV>&A)
{
  MVT::Assign(*(A.X_),*(X_));

  *Y_ = *A.Y_; // This is a well-defined operator for SerialDenseMatrix
}



// Replace each element of the vectors in *this with alpha.
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvInit (ScalarType alpha)
{
  MVT::MvInit(*X_,alpha);
  Y_->putScalar(alpha);
}



// Prints the multivector to an output stream
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvPrint (std::ostream &os) const
{
  int nrows = Y_->numRows();
  int ncols = Y_->numCols();

  os << "Object SaddleContainer" << std::endl;
  os << "X\n";
//  X_->describe(*(Teuchos::VerboseObjectBase::getDefaultOStream()),Teuchos::VERB_EXTREME);
//  os << "X\n" << *X_ << std::endl;
  
  os << "Y\n";
  for(int r=0; r<nrows; r++)
  {
    for(int c=0; c<ncols; c++)
      os << "Y[" << r << "][" << c << "]=" << (*Y_)[c][r] << std::endl;
  }
}

} // End namespace Experimental

template<class ScalarType, class MV >
class MultiVecTraits<ScalarType,Experimental::SaddleContainer<ScalarType, MV> >
{
public:
  static RCP<Experimental::SaddleContainer<ScalarType, MV> > Clone( const Experimental::SaddleContainer<ScalarType, MV>& mv, const int numvecs )
    { return rcp( const_cast<Experimental::SaddleContainer<ScalarType, MV>&>(mv).Clone(numvecs) ); }

  static RCP<Experimental::SaddleContainer<ScalarType, MV> > CloneCopy( const Experimental::SaddleContainer<ScalarType, MV>& mv )
    { return rcp( const_cast<Experimental::SaddleContainer<ScalarType, MV>&>(mv).CloneCopy() ); }

  static RCP<Experimental::SaddleContainer<ScalarType, MV> > CloneCopy( const Experimental::SaddleContainer<ScalarType, MV>& mv, const std::vector<int>& index )
    { return rcp( const_cast<Experimental::SaddleContainer<ScalarType, MV>&>(mv).CloneCopy(index) ); }
 
  static int GetVecLength( const Experimental::SaddleContainer<ScalarType, MV>& mv )
    { return mv.GetVecLength(); }

  static int GetNumberVecs( const Experimental::SaddleContainer<ScalarType, MV>& mv )
    { return mv.GetNumberVecs(); }

  static void MvTimesMatAddMv( ScalarType alpha, const Experimental::SaddleContainer<ScalarType, MV>& A, 
                               const Teuchos::SerialDenseMatrix<int,ScalarType>& B, 
                               ScalarType beta, Experimental::SaddleContainer<ScalarType, MV>& mv )
    { mv.MvTimesMatAddMv(alpha, A, B, beta); }

  static void MvAddMv( ScalarType alpha, const Experimental::SaddleContainer<ScalarType, MV>& A, ScalarType beta, const Experimental::SaddleContainer<ScalarType, MV>& B, Experimental::SaddleContainer<ScalarType, MV>& mv )
    { mv.MvAddMv(alpha, A, beta, B); }

  static void MvTransMv( ScalarType alpha, const Experimental::SaddleContainer<ScalarType, MV>& A, const Experimental::SaddleContainer<ScalarType, MV>& mv, Teuchos::SerialDenseMatrix<int,ScalarType>& B)
    { mv.MvTransMv(alpha, A, B); }

  static void MvDot( const Experimental::SaddleContainer<ScalarType, MV>& mv, const Experimental::SaddleContainer<ScalarType, MV>& A, std::vector<ScalarType> & b)
    { mv.MvDot( A, b); }

  static void MvScale ( Experimental::SaddleContainer<ScalarType, MV>& mv, ScalarType alpha )
    { mv.MvScale( alpha ); }

  static void MvScale ( Experimental::SaddleContainer<ScalarType, MV>& mv, const std::vector<ScalarType>& alpha )
    { mv.MvScale( alpha ); }
 
  static void MvNorm( const Experimental::SaddleContainer<ScalarType, MV>& mv, std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> & normvec)
    { mv.MvNorm(normvec); }

  static void SetBlock( const Experimental::SaddleContainer<ScalarType, MV>& A, const std::vector<int>& index, Experimental::SaddleContainer<ScalarType, MV>& mv )
    { mv.SetBlock(A, index); }

  static void Assign( const Experimental::SaddleContainer<ScalarType, MV>& A, Experimental::SaddleContainer<ScalarType, MV>& mv )
    { mv.Assign(A); }

  static void MvRandom( Experimental::SaddleContainer<ScalarType, MV>& mv )
    { mv.MvRandom(); }

  static void MvInit( Experimental::SaddleContainer<ScalarType, MV>& mv, ScalarType alpha = Teuchos::ScalarTraits<ScalarType>::zero() )
    { mv.MvInit(alpha); }

  static void MvPrint( const Experimental::SaddleContainer<ScalarType, MV>& mv, std::ostream& os )
    { mv.MvPrint(os); }
};

} // end namespace Anasazi

#endif
