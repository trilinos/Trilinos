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

/*! \file AnasaziSaddleOperator.hpp
 *  \brief An operator of the form [A Y; Y' 0] where A is a sparse matrix and Y a multivector.
 *
 *  Used by TraceMin to solve the saddle point problem.
*/

#ifndef ANASAZI_SADDLE_OPERATOR_HPP
#define ANASAZI_SADDLE_OPERATOR_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziSaddleContainer.hpp"
#include "AnasaziTraceMinRitzOp.hpp"

#include "Teuchos_SerialDenseSolver.hpp"

using Teuchos::RCP;

enum PrecType {NO_PREC, BD_PREC, UT_PREC};

namespace Anasazi {
namespace Experimental {

template <class ScalarType, class MV, class OP>
class SaddleOperator : public TraceMinOp<ScalarType,SaddleContainer<ScalarType,MV>,OP>
{
public:
	// Default constructor
	SaddleOperator( ) { };
	SaddleOperator( const Teuchos::RCP<OP> A, const Teuchos::RCP<const MV> B, PrecType pt=NO_PREC );

	// Applies the saddle point operator to a "multivector"
	void Apply(const SaddleContainer<ScalarType,MV>& X, SaddleContainer<ScalarType,MV>& Y) const;

  void removeIndices(const std::vector<int>& indicesToRemove) { A_->removeIndices(indicesToRemove); };

private:
	// A is the 1-1 block, and B the 1-2 block
	Teuchos::RCP<OP> A_;
	Teuchos::RCP<const MV> B_;
	Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > Schur_;
	PrecType pt_;
	
	typedef Anasazi::MultiVecTraits<ScalarType,MV> MVT;
};



// Default constructor
template <class ScalarType, class MV, class OP>
SaddleOperator<ScalarType, MV, OP>::SaddleOperator( const Teuchos::RCP<OP> A, const Teuchos::RCP<const MV> B, PrecType pt )
{
	// Get a pointer to A and B
	A_ = A;
	B_ = B;
	pt_ = pt;

	if(pt == BD_PREC || pt == UT_PREC)
	{
		// Form the Schur complement
		int nvecs = MVT::GetNumberVecs(*B);
		Teuchos::RCP<MV> AinvB = MVT::Clone(*B,nvecs);
		Schur_ = rcp(new Teuchos::SerialDenseMatrix<int,ScalarType>(nvecs,nvecs));

		A_->Apply(*B_,*AinvB);

		MVT::MvTransMv(1., *B_, *AinvB, *Schur_);
	}
}

// Applies the saddle point operator to a "multivector"
template <class ScalarType, class MV, class OP>
void SaddleOperator<ScalarType, MV, OP>::Apply(const SaddleContainer<ScalarType,MV>& X, SaddleContainer<ScalarType,MV>& Y) const
{
	if(pt_ == NO_PREC)
	{
		// trans does literally nothing, because the operator is symmetric
		// Y.bottom = B'X.top
		MVT::MvTransMv(1., *B_, *(X.X_), *Y.Y_);
	
		// Y.top = A*X.top+B*X.bottom
		A_->Apply(*(X.X_), *(Y.X_));
		MVT::MvTimesMatAddMv(1., *B_, *X.Y_, 1., *(Y.X_));
	}
	else if(pt_ == BD_PREC)
	{
		Teuchos::SerialDenseSolver<int,ScalarType> MySolver;

		// Solve A Y.X = X.X
		A_->Apply(*(X.X_),*(Y.X_));

		// So, let me tell you a funny story about how the SerialDenseSolver destroys the original matrix...
		Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > localSchur = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int,ScalarType>(*Schur_));

		// Solve the small system
		MySolver.setMatrix(localSchur);
		MySolver.setVectors(Y.Y_, X.Y_);
		MySolver.solve();
	}
	else if(pt_ == UT_PREC)
	{
		std::cout << "block upper triangular preconditioning is not implemented yet.  Sorry for any inconvenience!\n";
	}
	else
	{
		std::cout << "Not a valid preconditioner type\n";
	}
}

} // End namespace Experimental

template<class ScalarType, class MV, class OP>
class OperatorTraits<ScalarType, Experimental::SaddleContainer<ScalarType,MV>, Experimental::SaddleOperator<ScalarType,MV,OP> >
{
public:
  static void Apply( const Experimental::SaddleOperator<ScalarType,MV,OP>& Op, 
                     const Experimental::SaddleContainer<ScalarType,MV>& x, 
                     Experimental::SaddleContainer<ScalarType,MV>& y)
		{ Op.Apply( x, y); };
};

} // end namespace Anasazi

#endif
