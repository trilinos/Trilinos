// File BelosEpetraOperator.hpp : This file implements a BelosEpetraOperator that derives from
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

///////////////////////////////////////////////////////////////
//--------template class BelosEpetraOperator--------------------
//
// This class will allow Belos to be called as an Epetra_Operator.
// Thus, it can use itself as a preconditioner if need be.  It can
// also be used as the inner iteration of Anasazi :)
//
///////////////////////////////////////////////////////////////

template <class TYPE>
class BelosEpetraOperator : public virtual Epetra_Operator {
public:
	// Constructor / Destructor
	BelosEpetraOperator( BelosPetraMat<TYPE>& mat, BelosPetraPrec<TYPE>& precond,
			    const string solver="BlockGMRES", const TYPE tol=1.0e-6, const int maxits=25, 
			    const int block=1, const int debuglevel=0, bool verbose=false );

	~BelosEpetraOperator();

	// Implement Epetra_Operator functionality

	// Attribute set methods.
	int SetUseTranspose( bool UseTranspose ) { return(-1); };

	// Mathematical Functions.
	int Apply( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const;
	int ApplyInverse( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const;
	double NormInf() const { return(0.0); };

	// Atribute access functions
	char* Label() const { return(Solver); };
	bool UseTranspose() const { return(false); };
	bool HasNormInf() const { return(false); };
	const Epetra_Comm& Comm() const { return(MyComm); };
	const Epetra_Map& OperatorDomainMap() const { return(DomainMap); };
	const Epetra_Map& OperatorRangeMap() const { return(RangeMap); };	
	   
private:
	BelosPetraMat<TYPE>& Mat;
	BelosPetraPrec<TYPE>& Prec;
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
// implementation of the BelosEpetraOperator class.
//
// Constructor.
//
template <class TYPE>
BelosEpetraOperator<TYPE>::BelosEpetraOperator( BelosPetraMat<TYPE>& mat, BelosPetraPrec<TYPE>& precond,
			    	const string solver, const TYPE tol, const int maxits, const int block, 
				const int debuglevel, bool verbose )
	: Mat(mat), Prec(precond), Solver(0), Tol(tol),
	  Maxits(maxits), BlkSz(block), DbgLvl(debuglevel), Vb(verbose),
	  MyComm(mat.GetMatrix().Comm()), DomainMap(mat.GetMatrix().OperatorDomainMap()),
	  RangeMap(mat.GetMatrix().OperatorRangeMap())
{
	Solver = new char[solver.length()+1];
	solver.copy(Solver,10,0);
	Solver[solver.length()] = 0;
}

template<class TYPE>
BelosEpetraOperator<TYPE>::~BelosEpetraOperator()
{
}

template<class TYPE>
int BelosEpetraOperator<TYPE>::Apply( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const
{
	TYPE zero = 0.0, one = 1.0;
	BelosPetraVec<TYPE> vec_X(X), vec_Y(Y);
	int NumRHS = vec_X.GetNumberVecs();
	//
	// Create solver and solve problem.  This is inefficient, an instance of the solver should
	// exist already and just be reset with a new RHS.
	//
	if (strcmp(Solver,"BlockGMRES")==0) {
		Belos::BlockGmres<TYPE> MyBlockGmres( Mat, Prec, vec_X, NumRHS, Tol, Maxits, BlkSz, Vb );
		MyBlockGmres.SetDebugLevel(DbgLvl);
		MyBlockGmres.Solve(Vb);
		MyBlockGmres.GetSolutions( vec_Y );
		MyBlockGmres.PrintResids(Vb);
	}
	if (strcmp(Solver,"BlockCG")==0) {
		Belos::BlockCG<TYPE> MyBlockCG( Mat, Prec, vec_X, NumRHS, Tol, Maxits, BlkSz, Vb );
		MyBlockCG.SetDebugLevel(DbgLvl);
		MyBlockCG.Solve(Vb);
		MyBlockCG.GetSolutions( vec_Y );
		MyBlockCG.PrintResids(Vb);
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
int BelosEpetraOperator<TYPE>::ApplyInverse( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const
{
	TYPE zero = 0.0, one = 1.0;
	BelosPetraVec<TYPE> vec_X(X), vec_Y(Y);
	int NumRHS = vec_X.GetNumberVecs();
	//
	// Create solver and solve problem.  This is inefficient, an instance of the solver should
	// exist already and just be reset with a new RHS.
	//
	if (strcmp(Solver,"BlockGMRES")==0) {
		Belos::BlockGmres<TYPE> MyBlockGmres( Mat, Prec, vec_X, NumRHS, Tol, Maxits, BlkSz, Vb );
		MyBlockGmres.SetDebugLevel(DbgLvl);
		MyBlockGmres.Solve(Vb);
		MyBlockGmres.GetSolutions( vec_Y );
		MyBlockGmres.PrintResids(Vb);
	}
	if (strcmp(Solver,"BlockCG")==0) {
		Belos::BlockCG<TYPE> MyBlockCG( Mat, Prec, vec_X, NumRHS, Tol, Maxits, BlkSz, Vb );
		MyBlockCG.SetDebugLevel(DbgLvl);
		MyBlockCG.Solve(Vb);
		MyBlockCG.GetSolutions( vec_Y );
		MyBlockCG.PrintResids(Vb);
	}
	// Copy solution into output vector Y.
     	Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&vec_Y); assert(vec_y!=NULL);

	// Y = vec_y
	int info = Y.Update( one, *vec_y, zero );		
	assert(info==0);

	// Assume a good return right now since Belos doesn't have return types yet.
	return(0);
}

// end of file BELOS_EPETRA_OPERATOR_HPP
#endif 

