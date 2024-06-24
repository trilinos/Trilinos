// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#ifdef HAVE_ANASAZI_BELOS
#include "BelosConfigDefs.hpp"
#include "BelosMultiVecTraits.hpp"
#endif

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
  const ScalarType ONE, ZERO;
  RCP<MV> upper_;
  RCP<SerialDenseMatrix> lowerRaw_;
  std::vector<int> indices_;

  RCP<SerialDenseMatrix> getLower() const;
  void setLower(const RCP<SerialDenseMatrix> lower);

public:
  // Constructors
  SaddleContainer( ) : ONE(Teuchos::ScalarTraits<ScalarType>::one()), ZERO(Teuchos::ScalarTraits<ScalarType>::zero()) { };
  SaddleContainer( const RCP<MV> X, bool eye=false );

  // Things that are necessary for compilation
  // Returns a clone of the current vector
  SaddleContainer<ScalarType, MV> * Clone(const int nvecs) const;
  // Returns a duplicate of the current vector
  SaddleContainer<ScalarType, MV> * CloneCopy() const;
  // Returns a duplicate of the specified vectors
  SaddleContainer<ScalarType, MV> * CloneCopy(const std::vector<int> &index) const;
  SaddleContainer<ScalarType, MV> * CloneViewNonConst(const std::vector<int> &index) const;
  SaddleContainer<ScalarType, MV> * CloneViewNonConst(const Teuchos::Range1D& index) const;
  const SaddleContainer<ScalarType, MV> * CloneView(const std::vector<int> &index) const;
  const SaddleContainer<ScalarType, MV> * CloneView(const Teuchos::Range1D& index) const;
  ptrdiff_t GetGlobalLength() const { return MVT::GetGlobalLength(*upper_) + lowerRaw_->numRows(); };
  int GetNumberVecs() const { return MVT::GetNumberVecs(*upper_); };
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



// THIS IS NEW!
template <class ScalarType, class MV>
RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > SaddleContainer<ScalarType, MV>::getLower() const
{
  if(indices_.empty())
  {
    return lowerRaw_;
  }

  int nrows = lowerRaw_->numRows();
  int ncols = indices_.size();

  RCP<SerialDenseMatrix> lower = rcp(new SerialDenseMatrix(nrows,ncols,false));

  for(int r=0; r<nrows; r++)
  {
    for(int c=0; c<ncols; c++)
    {
      (*lower)(r,c) = (*lowerRaw_)(r,indices_[c]);
    }
  }

  return lower;
}



// THIS IS NEW!
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::setLower(const RCP<SerialDenseMatrix> lower)
{
  // If the indices are empty, lower points to lowerRaw
  if(indices_.empty())
  {
    return;
  }

  int nrows = lowerRaw_->numRows();
  int ncols = indices_.size();

  for(int r=0; r<nrows; r++)
  {
    for(int c=0; c<ncols; c++)
    {
      (*lowerRaw_)(r,indices_[c]) = (*lower)(r,c);
    }
  }
}



// Constructor
template <class ScalarType, class MV>
SaddleContainer<ScalarType, MV>::SaddleContainer( const RCP<MV> X, bool eye )
: ONE(Teuchos::ScalarTraits<ScalarType>::one()), ZERO(Teuchos::ScalarTraits<ScalarType>::zero())
{
  int nvecs = MVT::GetNumberVecs(*X);

  if(eye)
  {
    // Initialize upper_ as all 0s
    upper_ = MVT::Clone(*X, nvecs);
    MVT::MvInit(*upper_);

    // Initialize Y to be I
    lowerRaw_ = rcp(new SerialDenseMatrix(nvecs,nvecs));
    for(int i=0; i < nvecs; i++)
      (*lowerRaw_)(i,i) = ONE;
  }
  else
  {
    // Point upper_ to X
    upper_ = X;

    // Initialize Y to be 0
    lowerRaw_ = rcp(new SerialDenseMatrix(nvecs,nvecs));
  }
}



// Returns a clone of the current vector
template <class ScalarType, class MV>
SaddleContainer<ScalarType, MV> * SaddleContainer<ScalarType, MV>::Clone(const int nvecs) const
{
  SaddleContainer<ScalarType, MV> * newSC = new SaddleContainer<ScalarType, MV>();

  newSC->upper_ = MVT::Clone(*upper_,nvecs);
  newSC->lowerRaw_ = rcp(new SerialDenseMatrix(lowerRaw_->numRows(),nvecs));

  return newSC;
}



// Returns a duplicate of the current vector
template <class ScalarType, class MV>
SaddleContainer<ScalarType, MV> * SaddleContainer<ScalarType, MV>::CloneCopy() const
{
  SaddleContainer<ScalarType, MV> * newSC = new SaddleContainer<ScalarType, MV>();

  newSC->upper_ = MVT::CloneCopy(*upper_);
  newSC->lowerRaw_ = getLower();

  return newSC;
}



// Returns a duplicate of the specified vectors
template <class ScalarType, class MV>
SaddleContainer<ScalarType, MV> * SaddleContainer<ScalarType, MV>::CloneCopy(const std::vector< int > &index) const
{
  SaddleContainer<ScalarType, MV> * newSC = new SaddleContainer<ScalarType, MV>();

  newSC->upper_ = MVT::CloneCopy(*upper_,index);

  int ncols = index.size();
  int nrows = lowerRaw_->numRows();
  RCP<SerialDenseMatrix> lower = getLower();
  newSC->lowerRaw_ = rcp(new SerialDenseMatrix(nrows,ncols));
  for(int c=0; c < ncols; c++)
  {
    for(int r=0; r < nrows; r++)
      (*newSC->lowerRaw_)(r,c) = (*lower)(r,index[c]);
  }

  return newSC;
}



// THIS IS NEW!
template <class ScalarType, class MV>
SaddleContainer<ScalarType, MV> * SaddleContainer<ScalarType, MV>::CloneViewNonConst(const std::vector<int> &index) const
{
  SaddleContainer<ScalarType, MV> * newSC = new SaddleContainer<ScalarType, MV>();

  newSC->upper_ = MVT::CloneViewNonConst(*upper_,index);

  newSC->lowerRaw_ = lowerRaw_;

  if(!indices_.empty())
  {
    newSC->indices_.resize(index.size());
    for(unsigned int i=0; i<index.size(); i++)
    {
      newSC->indices_[i] = indices_[index[i]];
    }
  }
  else
  {
    newSC->indices_ = index;
  }

  return newSC;
}


template <class ScalarType, class MV>
SaddleContainer<ScalarType, MV> * SaddleContainer<ScalarType, MV>::CloneViewNonConst(const Teuchos::Range1D& index) const
{
  SaddleContainer<ScalarType, MV> * newSC = new SaddleContainer<ScalarType, MV>();

  newSC->upper_ = MVT::CloneViewNonConst(*upper_,index);

  newSC->lowerRaw_ = lowerRaw_;

  newSC->indices_.resize(index.size());
  for(unsigned int i=0; i<index.size(); i++)
  {
    newSC->indices_[i] = indices_[index.lbound()+i];
  }

  return newSC;
}



// THIS IS NEW!
template <class ScalarType, class MV>
const SaddleContainer<ScalarType, MV> * SaddleContainer<ScalarType, MV>::CloneView(const std::vector<int> &index) const
{
  SaddleContainer<ScalarType, MV> * newSC = new SaddleContainer<ScalarType, MV>();

  newSC->upper_ = MVT::CloneViewNonConst(*upper_,index);

  newSC->lowerRaw_ = lowerRaw_;

  if(!indices_.empty())
  {
    newSC->indices_.resize(index.size());
    for(unsigned int i=0; i<index.size(); i++)
    {
      newSC->indices_[i] = indices_[index[i]];
    }
  }
  else
  {
    newSC->indices_ = index;
  }

  return newSC;
}


template <class ScalarType, class MV>
const SaddleContainer<ScalarType, MV> * SaddleContainer<ScalarType, MV>::CloneView(const Teuchos::Range1D& index) const
{
  SaddleContainer<ScalarType, MV> * newSC = new SaddleContainer<ScalarType, MV>();

  newSC->upper_ = MVT::CloneViewNonConst(*upper_,index);

  newSC->lowerRaw_ = lowerRaw_;

  newSC->indices_.resize(index.size());
  for(unsigned int i=0; i<index.size(); i++)
  {
    newSC->indices_[i] = indices_[index.lbound()+i];
  }

  return newSC;
}



// Update *this with alpha * A * B + beta * (*this)
// THIS IS NEW!
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvTimesMatAddMv(ScalarType alpha, const SaddleContainer<ScalarType,MV> &A,
                                                      const Teuchos::SerialDenseMatrix<int, ScalarType> &B,
                                                      ScalarType beta)
{
  MVT::MvTimesMatAddMv(alpha,*(A.upper_),B,beta,*upper_);
  RCP<SerialDenseMatrix> lower = getLower();
  RCP<SerialDenseMatrix> Alower = A.getLower();
  lower->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,alpha,*Alower,B,beta);
  setLower(lower);
}



// Replace *this with alpha * A + beta * B
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvAddMv(ScalarType alpha, const SaddleContainer<ScalarType,MV>& A,
                                              ScalarType beta,  const SaddleContainer<ScalarType,MV>& B)
{
  MVT::MvAddMv(alpha, *(A.upper_), beta, *(B.upper_), *upper_);

  RCP<SerialDenseMatrix> lower = getLower();
  RCP<SerialDenseMatrix> Alower = A.getLower();
  RCP<SerialDenseMatrix> Blower = B.getLower();

  //int ncolsA = Alower->numCols(); // unused
  //int ncolsThis = lower->numCols(); // unused
  //int nrows = lower->numRows(); // unused

  // Y = alpha A
  lower->assign(*Alower);
  if(alpha != ONE)
    lower->scale(alpha);
  // Y += beta B
  if(beta == ONE)
    *lower += *Blower;
  else if(beta == -ONE)
    *lower -= *Blower;
  else if(beta != ZERO)
  {
    SerialDenseMatrix scaledB(*Blower);
    scaledB.scale(beta);
    *lower += *Blower;
  }

  setLower(lower);
}



// Scale the vectors by alpha
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvScale( ScalarType alpha )
{
  MVT::MvScale(*upper_, alpha);


  RCP<SerialDenseMatrix> lower = getLower();
  lower->scale(alpha);
  setLower(lower);
}



// Scale the i-th vector by alpha[i]
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvScale( const std::vector<ScalarType>& alpha )
{
  MVT::MvScale(*upper_, alpha);

  RCP<SerialDenseMatrix> lower = getLower();

  int nrows = lower->numRows();
  int ncols = lower->numCols();

  for(int c=0; c<ncols; c++)
  {
    for(int r=0; r<nrows; r++)
      (*lower)(r,c) *= alpha[c];
  }

  setLower(lower);
}



// Compute a dense matrix B through the matrix-matrix multiply alpha * A^H * (*this)
// THIS IS NEW!
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvTransMv (ScalarType alpha, const SaddleContainer<ScalarType, MV>& A,
                                                 Teuchos::SerialDenseMatrix< int, ScalarType >& B) const
{
  MVT::MvTransMv(alpha,*(A.upper_),*upper_,B);
  RCP<SerialDenseMatrix> lower = getLower();
  RCP<SerialDenseMatrix> Alower = A.getLower();
  B.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,alpha,*(Alower),*lower,ONE);
}



// Compute a vector b where the components are the individual dot-products, i.e.b[i] = A[i]^H*this[i] where A[i] is the i-th column of A.
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvDot (const SaddleContainer<ScalarType, MV>& A, std::vector<ScalarType> &b) const
{
  MVT::MvDot(*upper_, *(A.upper_), b);

  RCP<SerialDenseMatrix> lower = getLower();
  RCP<SerialDenseMatrix> Alower = A.getLower();

  int nrows = lower->numRows();
  int ncols = lower->numCols();

  for(int c=0; c < ncols; c++)
  {
    for(int r=0; r < nrows; r++)
    {
      b[c] += ((*Alower)(r,c) * (*lower)(r,c));
    }
  }
}



// Compute the 2-norm of each individual vector
// THIS IS NEW!
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvNorm ( std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &normvec) const
{
  // TODO: Make this better
  MvDot(*this,normvec);
  for(unsigned int i=0; i<normvec.size(); i++)
    normvec[i] = sqrt(normvec[i]);
}



// Copy the vectors in A to a set of vectors in *this. The numvecs vectors in
// A are copied to a subset of vectors in *this indicated by the indices given
// in index.
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::SetBlock (const SaddleContainer<ScalarType, MV>& A, const std::vector<int> &index)
{
  MVT::SetBlock(*(A.upper_), index, *upper_);

  RCP<SerialDenseMatrix> lower = getLower();
  RCP<SerialDenseMatrix> Alower = A.getLower();

  int nrows = lower->numRows();

  int nvecs = index.size();
  for(int c=0; c<nvecs; c++)
  {
    for(int r=0; r<nrows; r++)
      (*lower)(r,index[c]) = (*Alower)(r,c);
  }

  setLower(lower);
}



// Deep copy.
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::Assign (const SaddleContainer<ScalarType, MV>&A)
{
  MVT::Assign(*(A.upper_),*(upper_));

  RCP<SerialDenseMatrix> lower = getLower();
  RCP<SerialDenseMatrix> Alower = A.getLower();

  *lower = *Alower; // This is a well-defined operator for SerialDenseMatrix

  setLower(lower);
}



// Fill the vectors in *this with random numbers.
// THIS IS NEW!
template <class ScalarType, class MV>
void  SaddleContainer<ScalarType, MV>::MvRandom ()
{
  MVT::MvRandom(*upper_);

  RCP<SerialDenseMatrix> lower = getLower();
  lower->random();
  setLower(lower);
}



// Replace each element of the vectors in *this with alpha.
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvInit (ScalarType alpha)
{
  MVT::MvInit(*upper_,alpha);

  RCP<SerialDenseMatrix> lower = getLower();
  lower->putScalar(alpha);
  setLower(lower);
}



// Prints the multivector to an output stream
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvPrint (std::ostream &os) const
{
  RCP<SerialDenseMatrix> lower = getLower();
  //int nrows = lower->numRows(); // unused
  //int ncols = lower->numCols(); // unused

  os << "Object SaddleContainer" << std::endl;
  os << "X\n";
  upper_->describe(*(Teuchos::VerboseObjectBase::getDefaultOStream()),Teuchos::VERB_EXTREME);
//  os << "X\n" << *upper_ << std::endl;

  os << "Y\n" << *lower << std::endl;
}

} // End namespace Experimental

template<class ScalarType, class MV >
class MultiVecTraits<ScalarType,Experimental::SaddleContainer<ScalarType, MV> >
{
typedef Experimental::SaddleContainer<ScalarType,MV>  SC;

public:
  static RCP<SC > Clone( const SC& mv, const int numvecs )
    { return rcp( const_cast<SC&>(mv).Clone(numvecs) ); }

  static RCP<SC > CloneCopy( const SC& mv )
    { return rcp( const_cast<SC&>(mv).CloneCopy() ); }

  static RCP<SC > CloneCopy( const SC& mv, const std::vector<int>& index )
    { return rcp( const_cast<SC&>(mv).CloneCopy(index) ); }

  static ptrdiff_t GetGlobalLength( const SC& mv )
    { return mv.GetGlobalLength(); }

  static int GetNumberVecs( const SC& mv )
    { return mv.GetNumberVecs(); }

  static void MvTimesMatAddMv( ScalarType alpha, const SC& A,
                               const Teuchos::SerialDenseMatrix<int,ScalarType>& B,
                               ScalarType beta, SC& mv )
    { mv.MvTimesMatAddMv(alpha, A, B, beta); }

  static void MvAddMv( ScalarType alpha, const SC& A, ScalarType beta, const SC& B, SC& mv )
    { mv.MvAddMv(alpha, A, beta, B); }

  static void MvTransMv( ScalarType alpha, const SC& A, const SC& mv, Teuchos::SerialDenseMatrix<int,ScalarType>& B)
    { mv.MvTransMv(alpha, A, B); }

  static void MvDot( const SC& mv, const SC& A, std::vector<ScalarType> & b)
    { mv.MvDot( A, b); }

  static void MvScale ( SC& mv, ScalarType alpha )
    { mv.MvScale( alpha ); }

  static void MvScale ( SC& mv, const std::vector<ScalarType>& alpha )
    { mv.MvScale( alpha ); }

  static void MvNorm( const SC& mv, std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> & normvec)
    { mv.MvNorm(normvec); }

  static void SetBlock( const SC& A, const std::vector<int>& index, SC& mv )
    { mv.SetBlock(A, index); }

  static void Assign( const SC& A, SC& mv )
    { mv.Assign(A); }

  static void MvRandom( SC& mv )
    { mv.MvRandom(); }

  static void MvInit( SC& mv, ScalarType alpha = Teuchos::ScalarTraits<ScalarType>::zero() )
    { mv.MvInit(alpha); }

  static void MvPrint( const SC& mv, std::ostream& os )
    { mv.MvPrint(os); }
};

} // end namespace Anasazi

#ifdef HAVE_ANASAZI_BELOS
namespace Belos
{

template<class ScalarType, class MV >
class MultiVecTraits< ScalarType, Anasazi::Experimental::SaddleContainer<ScalarType,MV> >
{
typedef Anasazi::Experimental::SaddleContainer<ScalarType,MV>  SC;
public:
  static RCP<SC > Clone( const SC& mv, const int numvecs )
    { return rcp( const_cast<SC&>(mv).Clone(numvecs) ); }

  static RCP<SC > CloneCopy( const SC& mv )
    { return rcp( const_cast<SC&>(mv).CloneCopy() ); }

  static RCP<SC > CloneCopy( const SC& mv, const std::vector<int>& index )
    { return rcp( const_cast<SC&>(mv).CloneCopy(index) ); }

  static RCP<SC> CloneViewNonConst( SC& mv, const std::vector<int>& index )
    { return rcp( mv.CloneViewNonConst(index) ); }

  static RCP<SC> CloneViewNonConst( SC& mv, const Teuchos::Range1D& index )
    { return rcp( mv.CloneViewNonConst(index) ); }

  static RCP<const SC> CloneView( const SC& mv, const std::vector<int> & index )
    { return rcp( mv.CloneView(index) ); }

  static RCP<const SC> CloneView( const SC& mv, const Teuchos::Range1D& index )
    { return rcp( mv.CloneView(index) ); }

  static ptrdiff_t GetGlobalLength( const SC& mv )
    { return mv.GetGlobalLength(); }

  static int GetNumberVecs( const SC& mv )
    { return mv.GetNumberVecs(); }

  static void MvTimesMatAddMv( ScalarType alpha, const SC& A,
                               const Teuchos::SerialDenseMatrix<int,ScalarType>& B,
                               ScalarType beta, SC& mv )
    { mv.MvTimesMatAddMv(alpha, A, B, beta); }

  static void MvAddMv( ScalarType alpha, const SC& A, ScalarType beta, const SC& B, SC& mv )
    { mv.MvAddMv(alpha, A, beta, B); }

  static void MvTransMv( ScalarType alpha, const SC& A, const SC& mv, Teuchos::SerialDenseMatrix<int,ScalarType>& B)
    { mv.MvTransMv(alpha, A, B); }

  static void MvDot( const SC& mv, const SC& A, std::vector<ScalarType> & b)
    { mv.MvDot( A, b); }

  static void MvScale ( SC& mv, ScalarType alpha )
    { mv.MvScale( alpha ); }

  static void MvScale ( SC& mv, const std::vector<ScalarType>& alpha )
    { mv.MvScale( alpha ); }

  // TODO: MAKE SURE TYPE == TWONORM
  static void MvNorm( const SC& mv, std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> & normvec, NormType type=TwoNorm)
    { mv.MvNorm(normvec); }

  static void SetBlock( const SC& A, const std::vector<int>& index, SC& mv )
    { mv.SetBlock(A, index); }

  static void Assign( const SC& A, SC& mv )
    { mv.Assign(A); }

  static void MvRandom( SC& mv )
    { mv.MvRandom(); }

  static void MvInit( SC& mv, ScalarType alpha = Teuchos::ScalarTraits<ScalarType>::zero() )
    { mv.MvInit(alpha); }

  static void MvPrint( const SC& mv, std::ostream& os )
    { mv.MvPrint(os); }

#ifdef HAVE_BELOS_TSQR
  /// \typedef tsqr_adaptor_type
  /// \brief TsqrAdaptor specialization for the multivector type MV.
  ///
  /// By default, we provide a "stub" implementation.  It has the
  /// right methods and typedefs, but its constructors and methods
  /// all throw std::logic_error.  If you plan to use TSQR in Belos
  /// (e.g., through TsqrOrthoManager or TsqrMatOrthoManager), and
  /// if your multivector type MV is neither Epetra_MultiVector nor
  /// Tpetra::MultiVector, you must implement a functional TSQR
  /// adapter.  Please refer to Epetra::TsqrAdapter (for
  /// Epetra_MultiVector) or Tpetra::TsqrAdaptor (for
  /// Tpetra::MultiVector) for examples.
  typedef Belos::details::StubTsqrAdapter<SC> tsqr_adaptor_type;
#endif // HAVE_BELOS_TSQR
};

} // end namespace Belos
#endif

#endif
