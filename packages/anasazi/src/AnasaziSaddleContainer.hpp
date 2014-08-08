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

using Teuchos::RCP;

namespace Anasazi {
namespace Experimental {

template <class ScalarType, class MV>
class SaddleContainer //: public Anasazi::SaddleContainer<ScalarType, MV>
{
public:
	// Constructors
	SaddleContainer( ) { };
	SaddleContainer( const Teuchos::RCP<MV> X, bool eye=false );

	// Things that are necessary for compilation
	// Returns a clone of the current vector
	SaddleContainer<ScalarType, MV> * Clone(const int nvecs) const;
	// Returns a duplicate of the current vector
	SaddleContainer<ScalarType, MV> * CloneCopy() const;
	// Returns a duplicate of the specified vectors
	SaddleContainer<ScalarType, MV> * CloneCopy(const std::vector< int > &index) const;
	// Returns a view of current vector (shallow copy)
	SaddleContainer<ScalarType, MV> * CloneViewNonConst(const std::vector< int > &index);
	// Returns a view of current vector (shallow copy)
	SaddleContainer<ScalarType, MV> * CloneView(const std::vector< int > &index) const;
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
	// Deep copy.  Although the documentation doesn't say it's required, it totally is.  Great.
	void Assign (const SaddleContainer<ScalarType, MV>&A);
	// Fill the vectors in *this with random numbers.
	void  MvRandom ();
	// Replace each element of the vectors in *this with alpha.
	void  MvInit (ScalarType alpha);
	// Prints the multivector to an output stream
	void MvPrint (std::ostream &os) const;
	
	Teuchos::RCP<MV> X_;
	//Teuchos::RCP<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > Y_;
	ScalarType ** Y_;
	int nyrows, nycols;

private:
	typedef Anasazi::MultiVecTraits<ScalarType,MV> MVT;
};



// Constructor
template <class ScalarType, class MV>
SaddleContainer<ScalarType, MV>::SaddleContainer( const Teuchos::RCP<MV> X, bool eye )
{
	int nvecs = MVT::GetNumberVecs(*X);

	if(eye)
	{
		// Initialize X_ as all 0s
		X_ = MVT::Clone(*X, nvecs);
		MVT::MvInit(*X_);
		
		// Initialize Y to be I
		Y_ = (ScalarType**) malloc(nvecs*sizeof(ScalarType*));
		for(int c=0; c < nvecs; c++)
		{
			Y_[c] = (ScalarType*) malloc(nvecs*sizeof(ScalarType));
			for(int r=0; r < nvecs; r++)
			{
				if(r == c)
					Y_[c][r] = 1.;
				else
					Y_[c][r] = 0.;
			}
		}
		nyrows = nvecs;
		nycols = nvecs;
	}
	else
	{
		// Point X_ to X
		X_ = X;

		// Initialize Y to be 0
		Y_ = (ScalarType**) malloc(nvecs*sizeof(ScalarType*));
		for(int c=0; c < nvecs; c++)
		{
			Y_[c] = (ScalarType*) malloc(nvecs*sizeof(ScalarType));
			for(int r=0; r < nvecs; r++)
				Y_[c][r] = 0.;
		}
		nyrows = nvecs;
		nycols = nvecs;
	}
}



// Returns a clone of the current vector
template <class ScalarType, class MV>
SaddleContainer<ScalarType, MV> * SaddleContainer<ScalarType, MV>::Clone(const int nvecs) const 
{	
	SaddleContainer<ScalarType, MV> * newSC = new SaddleContainer<ScalarType, MV>();
	
	newSC->X_ = MVT::Clone(*X_,nvecs);
		
	newSC->Y_ = (ScalarType**) malloc(nvecs*sizeof(ScalarType*));
	for(int c=0; c < nvecs; c++)
		(newSC->Y_)[c] = (ScalarType*) malloc(nyrows*sizeof(ScalarType));
	newSC->nyrows = nyrows;
	newSC->nycols = nvecs;	

	return newSC;
}



// Returns a duplicate of the current vector
template <class ScalarType, class MV>
SaddleContainer<ScalarType, MV> * SaddleContainer<ScalarType, MV>::CloneCopy() const 
{
	SaddleContainer<ScalarType, MV> * newSC = new SaddleContainer<ScalarType, MV>();
	
	newSC->X_ = MVT::CloneCopy(*X_);

	newSC->Y_ = (ScalarType**) malloc(nycols*sizeof(ScalarType*));
	for(int c=0; c < nycols; c++)
	{
		(newSC->Y_)[c] = (ScalarType*) malloc(nyrows*sizeof(ScalarType));
		for(int r=0; r < nyrows; r++)
			(newSC->Y_)[c][r] = Y_[c][r];
	}
	newSC->nyrows = nyrows;
	newSC->nycols = nycols;
	
	return newSC;
}



// Returns a duplicate of the specified vectors
template <class ScalarType, class MV>
SaddleContainer<ScalarType, MV> * SaddleContainer<ScalarType, MV>::CloneCopy(const std::vector< int > &index) const
{
	SaddleContainer<ScalarType, MV> * newSC = new SaddleContainer<ScalarType, MV>();
	
	newSC->X_ = MVT::CloneCopy(*X_,index);
	
	int ncols = index.size();
	newSC->Y_ = (ScalarType**) malloc(ncols*sizeof(ScalarType*));
	for(int c=0; c < ncols; c++)
	{
		(newSC->Y_)[c] = (ScalarType*) malloc(nyrows*sizeof(ScalarType));
		for(int r=0; r < nyrows; r++)
			(newSC->Y_)[c][r] = Y_[index[c]][r];
	}
	newSC->nyrows = nyrows;
	newSC->nycols = ncols;
	
	return newSC;
}



// Returns a view of current vector (shallow copy)
template <class ScalarType, class MV>
SaddleContainer<ScalarType, MV> * SaddleContainer<ScalarType, MV>::CloneViewNonConst(const std::vector< int > &index)
{	
	SaddleContainer<ScalarType, MV> * newSC = new SaddleContainer<ScalarType, MV>();
	
	newSC->X_ = MVT::CloneViewNonConst(*X_,index);
	
	int ncols = index.size();
	newSC->Y_ = (ScalarType**) malloc(ncols*sizeof(ScalarType*));
	for(int c=0; c < ncols; c++)
	{
		(newSC->Y_)[c] = Y_[index[c]];
	}
	newSC->nyrows = nyrows;
	newSC->nycols = ncols;
	
	return newSC;
}



// Returns a view of current vector (shallow copy)
template <class ScalarType, class MV>
SaddleContainer<ScalarType, MV> * SaddleContainer<ScalarType, MV>::CloneView(const std::vector< int > &index) const
{	
	SaddleContainer<ScalarType, MV> * newSC = new SaddleContainer<ScalarType, MV>();
	
	newSC->X_ = MVT::CloneViewNonConst(*X_,index);
	
	int ncols = index.size();
	newSC->Y_ = (ScalarType**) malloc(ncols*sizeof(ScalarType*));
	for(int c=0; c < ncols; c++)
	{
		(newSC->Y_)[c] = Y_[index[c]];
	}
	newSC->nyrows = nyrows;
	newSC->nycols = ncols;
	
	return newSC;
}



template <class ScalarType, class MV>
int SaddleContainer<ScalarType, MV>::GetVecLength() const
{
	return (MVT::GetVecLength(*X_)+nyrows);
}



// Update *this with alpha * A * B + beta * (*this)
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvTimesMatAddMv(ScalarType alpha, const SaddleContainer<ScalarType,MV> &A, 
                                                      const Teuchos::SerialDenseMatrix<int, ScalarType> &B, 
													  ScalarType beta)
{
	MVT::MvTimesMatAddMv(alpha, *(A.X_), B, beta, *X_);
	
	// Copy over Y
	Teuchos::SerialDenseMatrix<int,ScalarType> locY(nyrows,nycols), locAY(A.nyrows,A.nycols);
	
	for(int c=0; c<nycols; c++)
		for(int r=0; r<nyrows; r++)
			locY(r,c) = Y_[c][r];
			
	for(int c=0; c<A.nycols; c++)
		for(int r=0; r<A.nyrows; r++)
			locAY(r,c) = A.Y_[c][r];

	// Do the multiplication with Y
	locY.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha, locAY, B, beta);
	
	// Copy Y back
	for(int c=0; c<nycols; c++)
		for(int r=0; r<nyrows; r++)
			Y_[c][r] = locY(r,c);	
}



// Replace *this with alpha * A + beta * B
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvAddMv(ScalarType alpha, const SaddleContainer<ScalarType,MV>& A, 
                                              ScalarType beta,  const SaddleContainer<ScalarType,MV>& B)
{
	MVT::MvAddMv(alpha, *(A.X_), beta, *(B.X_), *X_);
		
	// check whether additional space is needed
	if(A.nycols > nycols)
	{
		Y_ = (ScalarType**) realloc(Y_,A.nycols*sizeof(ScalarType*));
		nycols = A.nycols;
	}
	// check whether space can be freed
	else if(A.nycols < nycols)
	{
		for(int i=A.nycols; i < nycols; i++)
			free(Y_[i]);
		Y_ = (ScalarType**) realloc(Y_,A.nycols*sizeof(ScalarType*));
		nycols = A.nycols;
	}

	// do the addition
	for(int c=0; c < nycols; c++)
	{
		for(int r=0; r < nyrows; r++)
		{
			Y_[c][r] = alpha*A.Y_[c][r] + beta*B.Y_[c][r];
		}
	}
}



// Scale the vectors by alpha
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvScale( ScalarType alpha )
{
	MVT::MvScale(*X_, alpha);
	
	for(int c=0; c<nycols; c++)
	{
		for(int r=0; r<nyrows; r++)
			Y_[c][r] *= alpha;
	}
}



// Scale the i-th vector by alpha[i]
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvScale( const std::vector<ScalarType>& alpha )
{
	MVT::MvScale(*X_, alpha);
	
	for(int c=0; c<nycols; c++)
	{
		for(int r=0; r<nyrows; r++)
			Y_[c][r] *= alpha[c];
	}
}



// Compute a dense matrix B through the matrix-matrix multiply alpha * A^H * (*this)
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvTransMv (ScalarType alpha, const SaddleContainer<ScalarType, MV>& A, 
                                                 Teuchos::SerialDenseMatrix< int, ScalarType >& B) const
{
	MVT::MvTransMv(alpha, *(A.X_), *X_, B);
	
	// Copy over Y
	Teuchos::SerialDenseMatrix<int,ScalarType> locY(nyrows,nycols), locAY(A.nyrows,A.nycols);
	
	for(int c=0; c<nycols; c++)
		for(int r=0; r<nyrows; r++)
			locY(r,c) = Y_[c][r];
			
	for(int c=0; c<A.nycols; c++)
		for(int r=0; r<A.nyrows; r++)
			locAY(r,c) = A.Y_[c][r];
	
	// Do the multiplication
	B.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, alpha, locAY, locY, 1.);
}



// Compute a vector b where the components are the individual dot-products, i.e.b[i] = A[i]^H*this[i] where A[i] is the i-th column of A. 
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvDot (const SaddleContainer<ScalarType, MV>& A, std::vector<ScalarType> &b) const
{
  std::cout << "in mvdot\n";
	MVT::MvDot(*X_, *(A.X_), b);
  std::cout << "nycols: " << nycols << std::endl;
  std::cout << "nyrows: " << nyrows << std::endl;
	
	for(int c=0; c < nycols; c++)
	{
		for(int r=0; r < nyrows; r++)
		{
			b[c] += (A.Y_[c][r] * Y_[c][r]);
      std::cout << "r: " << r << " and c: " << c << " give us b[c] = " << b[c] << std::endl;
		}
	}

  std::cout << "b.size(): " << b.size() << std::endl;
}



// Compute the 2-norm of each individual vector
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvNorm ( std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &normvec) const
{	
/*	ScalarType temp;

	// sum of absolute values
	if(type == Anasazi::OneNorm)
	{
		MVT::MvNorm(*X_, normvec, type);
	
		for(int c=0; c < nycols; c++)
		{
			for(int r=0; r < nyrows; r++)
			{
				temp = Y_[c][r];
				normvec[c] += abs(temp);
			}
		}
	
	// square root of dot product
	else if(type == Anasazi::TwoNorm)
	{*/
		MvDot(*this,normvec);
		
		int nvecs = normvec.size();
		for(int i=0; i<nvecs; i++)
			normvec[i] = sqrt(normvec[i]);
/*	}
	// max value
	else if(type == Anasazi::InfNorm)
	{
		MVT::MvNorm(*X_, normvec, type);
		
		for(int c=0; c < nycols; c++)
		{
			for(int r=0; r < nyrows; r++)
			{
				temp = Y_[c][r];
				normvec[c] = max(normvec[c], abs(temp));
			}
		}	
	}
	else
		std::cout << "Error: that is not a valid norm type\n";*/
}



// Copy the vectors in A to a set of vectors in *this. The numvecs vectors in 
// A are copied to a subset of vectors in *this indicated by the indices given 
// in index.
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::SetBlock (const SaddleContainer<ScalarType, MV>& A, const std::vector<int> &index)
{
	MVT::SetBlock(*X_, index, *(A.X_));
	
	int nvecs = index.size();
	for(int c=0; c<nvecs; c++)
	{
		for(int r=0; r<nyrows; r++)
			Y_[index[c]][r] = A.Y_[c][r];
	}
}



// Deep copy.  Although the documentation doesn't say it's required, it totally is.  Great.
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::Assign (const SaddleContainer<ScalarType, MV>&A)
{
	// hopefully this performs a deep copy...
	(*X_) = (*(A.X_));
	
	// check whether additional space is needed
	if(A.nycols > nycols)
	{
		Y_ = (ScalarType**) realloc(Y_,A.nycols*sizeof(ScalarType*));
		nycols = A.nycols;
	}
	// check whether space can be freed
	else if(A.nycols < nycols)
	{
		for(int i=A.nycols; i < nycols; i++)
			free(Y_[i]);
		Y_ = (ScalarType**) realloc(Y_,A.nycols*sizeof(ScalarType*));
		nycols = A.nycols;
	}
	
	// copy over the values
	for(int c=0; c<nycols; c++)
	{
		for(int r=0; r<nyrows; r++)
			Y_[c][r] = A.Y_[c][r];
	}
}



// Fill the vectors in *this with random numbers.
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvRandom ()
{
	X_->Random();
	
	// get random values
	Teuchos::SerialDenseMatrix<int,ScalarType> locY(nyrows,nycols);
	locY.random();
	
	// put the data back in Y
	for(int c=0; c<nycols; c++)
	{
		for(int r=0; r<nyrows; r++)
			Y_[c][r] = locY(r,c);
	}
}



// Replace each element of the vectors in *this with alpha.
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvInit (ScalarType alpha)
{
  MVT::MvInit(*X_,alpha);
	
	for(int c=0; c<nycols; c++)
	{
		for(int r=0; r<nyrows; r++)
			Y_[c][r] = alpha;
	}
}



// Prints the multivector to an output stream
template <class ScalarType, class MV>
void SaddleContainer<ScalarType, MV>::MvPrint (std::ostream &os) const
{
	os << "Object SaddleContainer" << std::endl;
	os << "X\n" << *X_ << std::endl;
	
	os << "Y\n";
	for(int r=0; r<nyrows; r++)
	{
		for(int c=0; c<nycols; c++)
			os << "Y[" << r << "][" << c << "]=" << Y_[c][r] << std::endl;
	}
}

} // End namespace Experimental

template<class ScalarType, class MV >
class MultiVecTraits<ScalarType,Experimental::SaddleContainer<ScalarType, MV> >
{
public:
	static Teuchos::RCP<Experimental::SaddleContainer<ScalarType, MV> > Clone( const Experimental::SaddleContainer<ScalarType, MV>& mv, const int numvecs )
		{ return Teuchos::rcp( const_cast<Experimental::SaddleContainer<ScalarType, MV>&>(mv).Clone(numvecs) ); }

	static Teuchos::RCP<Experimental::SaddleContainer<ScalarType, MV> > CloneCopy( const Experimental::SaddleContainer<ScalarType, MV>& mv )
		{ return Teuchos::rcp( const_cast<Experimental::SaddleContainer<ScalarType, MV>&>(mv).CloneCopy() ); }

	static Teuchos::RCP<Experimental::SaddleContainer<ScalarType, MV> > CloneCopy( const Experimental::SaddleContainer<ScalarType, MV>& mv, const std::vector<int>& index )
		{ return Teuchos::rcp( const_cast<Experimental::SaddleContainer<ScalarType, MV>&>(mv).CloneCopy(index) ); }

	static Teuchos::RCP<Experimental::SaddleContainer<ScalarType, MV> > CloneViewNonConst( Experimental::SaddleContainer<ScalarType, MV>& mv, const std::vector<int>& index )
		{ return Teuchos::rcp( mv.CloneViewNonConst(index) ); }

	static Teuchos::RCP<const Experimental::SaddleContainer<ScalarType, MV> > CloneView( const Experimental::SaddleContainer<ScalarType, MV>& mv, const std::vector<int>& index )
		{ return Teuchos::rcp( const_cast<Experimental::SaddleContainer<ScalarType, MV>&>(mv).CloneView(index) ); }
 
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
