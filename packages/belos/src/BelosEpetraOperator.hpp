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

// File BelosEpetraOperator.hpp :   This file implements a BelosEpetraOperator that derives from
// 					the EpetraOperator class, so Belos can be integrated into
//					other codes as an abstract operator.
//
#ifndef BELOS_EPETRA_OPERATOR_HPP
#define BELOS_EPETRA_OPERATOR_HPP

#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "BelosPetraInterface.hpp"

#include "BelosBlockGmres.hpp"
#include "BelosBlockCG.hpp"

/*! \class Belos::EpetraOperator
 	\brief This class  derives from	the Epetra_Operator class, so Belos can be 
	integrated into other codes as an abstract operator.
*/

namespace Belos {

///////////////////////////////////////////////////////////////
//--------template class BelosEpetraOperator--------------------
//
// This class will allow Belos to be called as an Epetra_Operator.
// Thus, it can use itself as a preconditioner if need be.  It can
// also be used as the inner iteration of Anasazi :)
//
///////////////////////////////////////////////////////////////

template <class TYPE>
class EpetraOperator : public virtual Epetra_Operator {
public:
	//@{ \name Constructor / Destructor

	//! Constructor
	EpetraOperator( PetraMat<TYPE>& mat, PetraPrec<TYPE>& precond,
			    const string solver="BlockGMRES", const TYPE tol=1.0e-6, const int maxits=25, 
			    const int block=1, const int debuglevel=0, bool verbose=false );

	//! Destructor
	~EpetraOperator();
	//@}

	//@{ \name Attribute methods.

	//! Set whether the operator or its inverse should be applied. [ This option is not implemented ]
	int SetUseTranspose( bool UseTranspose ) { return(-1); };
	//@}

	//@{ \name Operator application methods.

	//! Apply the operator.
	int Apply( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const;

	//! Apply the operator's inverse.
	int ApplyInverse( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const;
	//@}

	//@{ \name Norm methods.

	//! Compute the infinity norm of the operator. [ This option is not implemented ]
	double NormInf() const { return(0.0); };
	//@}

	//@{ \name Attribute access functions

	//! Return the label of the operator.
	char* Label() const { return(Solver); };

	//! Return whether the operator is using the transpose.
	bool UseTranspose() const { return(false); };

	//! Return whether the infinity norm is available for this operator.
	bool HasNormInf() const { return(false); };

	//! Return the communicator for this operator.
	const Epetra_Comm& Comm() const { return(MyComm); };

	//! Return the domain map for this operator.
	const Epetra_Map& OperatorDomainMap() const { return(DomainMap); };

	//! Return the range map for this operator.
	const Epetra_Map& OperatorRangeMap() const { return(RangeMap); };	
	//@}	   
private:
	PetraMat<TYPE>& Mat;
	PetraPrec<TYPE>& Prec;
	const Epetra_Comm& MyComm;
	const Epetra_Map& DomainMap;
	const Epetra_Map& RangeMap;
	const TYPE Tol;
	const int Maxits, BlkSz, DbgLvl;
	char* Solver;
	bool Vb;
};
//--------------------------------------------------------------
//
// implementation of the EpetraOperator class.
//
// Constructor.
//
template <class TYPE>
EpetraOperator<TYPE>::EpetraOperator( PetraMat<TYPE>& mat, PetraPrec<TYPE>& precond,
			    	const string solver, const TYPE tol, const int maxits, const int block, 
				const int debuglevel, bool verbose )
	: Mat(mat), Prec(precond), Solver(0), Tol(tol),
	  Maxits(maxits), BlkSz(block), DbgLvl(debuglevel), Vb(verbose),
	  MyComm(mat.GetMatrix().Comm()), DomainMap(mat.GetMatrix().OperatorDomainMap()),
	  RangeMap(mat.GetMatrix().OperatorRangeMap())
{
        // Copy string to character array.  
        // Not using conversion routine copy() because it's not supported by RW on Janus. (HKT 11/13/2003) 
        Solver = new char[solver.length()+1];
        for (int i=0; i<solver.length()+1; i++) {
            Solver[i] = solver[i];
        } 
	Solver[solver.length()] = 0;
}

template<class TYPE>
EpetraOperator<TYPE>::~EpetraOperator()
{
}

template<class TYPE>
int EpetraOperator<TYPE>::Apply( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const
{
	TYPE zero = 0.0, one = 1.0;
	PetraVec<TYPE> vec_X(X), vec_Y(Y);
	int NumRHS = vec_X.GetNumberVecs();
	//
	// Create solver and solve problem.  This is inefficient, an instance of the solver should
	// exist already and just be reset with a new RHS.
	//
	if (strcmp(Solver,"BlockGMRES")==0) {
		BlockGmres<TYPE> MyBlockGmres( Mat, Prec, vec_X, NumRHS, Tol, Maxits, BlkSz, Vb );
		MyBlockGmres.SetDebugLevel(DbgLvl);
		MyBlockGmres.Solve(Vb);
		MyBlockGmres.GetSolutions( vec_Y );
		MyBlockGmres.TrueResiduals(Vb);
	}
	if (strcmp(Solver,"BlockCG")==0) {
		BlockCG<TYPE> MyBlockCG( Mat, Prec, vec_X, NumRHS, Tol, Maxits, BlkSz, Vb );
		MyBlockCG.SetDebugLevel(DbgLvl);
		MyBlockCG.Solve(Vb);
		MyBlockCG.GetSolutions( vec_Y );
		MyBlockCG.TrueResiduals(Vb);
	}
	// Copy solution into output vector Y.
     	Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&vec_Y); assert(vec_y!=NULL);

	// Y = vec_y
	int info = Y.Update( one, *vec_y, zero );		
	assert(info==0);

	// Assume a good return right now since Belos doesn't have return types yet.
	return(0);
}

template<class TYPE>
int EpetraOperator<TYPE>::ApplyInverse( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const
{
	TYPE zero = 0.0, one = 1.0;
	PetraVec<TYPE> vec_X(X), vec_Y(Y);
	int NumRHS = vec_X.GetNumberVecs();
	//
	// Create solver and solve problem.  This is inefficient, an instance of the solver should
	// exist already and just be reset with a new RHS.
	//
	if (strcmp(Solver,"BlockGMRES")==0) {
		BlockGmres<TYPE> MyBlockGmres( Mat, Prec, vec_X, NumRHS, Tol, Maxits, BlkSz, Vb );
		MyBlockGmres.SetDebugLevel(DbgLvl);
		MyBlockGmres.Solve(Vb);
		MyBlockGmres.GetSolutions( vec_Y );
		MyBlockGmres.TrueResiduals(Vb);
	}
	if (strcmp(Solver,"BlockCG")==0) {
		BlockCG<TYPE> MyBlockCG( Mat, Prec, vec_X, NumRHS, Tol, Maxits, BlkSz, Vb );
		MyBlockCG.SetDebugLevel(DbgLvl);
		MyBlockCG.Solve(Vb);
		MyBlockCG.GetSolutions( vec_Y );
		MyBlockCG.TrueResiduals(Vb);
	}
	// Copy solution into output vector Y.
     	Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&vec_Y); assert(vec_y!=NULL);

	// Y = vec_y
	int info = Y.Update( one, *vec_y, zero );		
	assert(info==0);

	// Assume a good return right now since Belos doesn't have return types yet.
	return(0);
}

} // end of Belos namespace

// end of file BELOS_EPETRA_OPERATOR_HPP
#endif 

