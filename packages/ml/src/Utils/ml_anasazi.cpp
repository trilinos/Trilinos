#include "ml_common.h"

#if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_ANASAZI)

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Teuchos_ParameterList.hpp"

#include "AnasaziPetraInterface.hpp"
#include "AnasaziBlockArnoldi.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziReturnType.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziMatrix.hpp"
#include "AnasaziOperator.hpp"

using namespace Teuchos;

#include "ml_anasazi.h"

using namespace Anasazi;

template <class TYPE> 
class MLMat : public virtual Matrix<TYPE> {
public:
  MLMat(const Epetra_Operator&, const Epetra_Vector &, const bool);
  ~MLMat();
  ReturnType ApplyMatrix ( const MultiVec<TYPE>& x, 
			   MultiVec<TYPE>& y ) const;
  void SetApplyMatrix(const bool ApplyMatrix) 
  {
    ApplyMatrix_ = ApplyMatrix;
  }
  
private:
  const Epetra_Operator & Mat_;
  const Epetra_Vector & Diagonal_;
  bool ApplyMatrix_;
};

template <class TYPE>
MLMat<TYPE>::MLMat(const Epetra_Operator & Matrix,
		   const Epetra_Vector & Diagonal,
		   const bool ApplyMatrix) 
  : Mat_(Matrix),
    Diagonal_(Diagonal),
    ApplyMatrix_(ApplyMatrix)
{
}

template <class TYPE>
MLMat<TYPE>::~MLMat() 
{
}

template <class TYPE>
ReturnType MLMat<TYPE>::ApplyMatrix ( const MultiVec<TYPE>& x, 
				      MultiVec<TYPE>& y ) const 
{
  int info=0;
  MultiVec<TYPE> & temp_x = const_cast<MultiVec<TYPE> &>(x);
  Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
  
  assert( vec_x!=NULL && vec_y!=NULL );
  //
  // Need to cast away constness because the member function Apply
  // is not declared const.
  //
  info=const_cast<Epetra_Operator&>(Mat_).Apply( *vec_x, *vec_y );
  if (info==0) {

    for( int j=0 ; j<vec_y->NumVectors() ; ++j ) 
      for( int i=0 ; i<vec_y->Map().NumMyElements() ; ++i ) (*vec_y)[j][i] *= Diagonal_[i];

    if( ApplyMatrix_ == false ) {
      vec_y->Update (1.0,*vec_x,-1.0);
    }
    
    return Ok; 
  } else { 
    return Failed; 
  }	
}

using namespace Teuchos;

int ML_Anasazi_Interface(const Epetra_RowMatrix * RowMatrix, int & NullSpaceDim,
			 double * & NullSpacePtr, ParameterList & List) 
{

  
  // copy diagonal in an epetra_vector, being interested in I-D^{-1} A
  Epetra_Vector * RowMatrixDiagonal = new Epetra_Vector(RowMatrix->RowMatrixRowMap());
  EPETRA_CHK_ERR( RowMatrix->ExtractDiagonalCopy(*RowMatrixDiagonal) );

  int NumMyElements = RowMatrixDiagonal->Map().NumMyElements();
  
  for( int i=0 ; i<NumMyElements ; ++i )
    if( (*RowMatrixDiagonal)[i] != 0.0 ) (*RowMatrixDiagonal)[i] = 1.0/(*RowMatrixDiagonal)[i];
    else  (*RowMatrixDiagonal)[i] = 0.0;
  
  int block = 1;
  int length = 20;
  int nev = NullSpaceDim;
  double tol = 1.0e-10;
  string which="LR";
  int restarts = 100;
  //int step = 1;
  int step = restarts*length*block;

  Anasazi::PetraVec<double> ivec(RowMatrix->RowMatrixRowMap(), block);
  ivec.MvRandom();

  /* ********************************************************************** */
  /* Phase 1: estimate lambda max                                           */
  /* ********************************************************************** */
  
  // Call the ctor that calls the petra ctor for a matrix
  MLMat<double> Amat1(*RowMatrix,*RowMatrixDiagonal,true);	
  Anasazi::Eigenproblem<double> MyProblem1(&Amat1, &ivec);

  // Initialize the Block Arnoldi solver
  Anasazi::BlockArnoldi<double> MyBlockArnoldi1(MyProblem1, tol, nev, length, block, 
					       which, step, restarts);
	
  // Inform the solver that the problem is symmetric
  //MyBlockArnoldi1.setSymmetric(true);
  MyBlockArnoldi1.setDebugLevel(0);

  //  MyBlockArnoldi.iterate(5);
  
  // Solve the problem to the specified tolerances or length
  MyBlockArnoldi1.solve();
  
  // Obtain results directly
  double * resids = MyBlockArnoldi1.getResiduals();
  double * evalr = MyBlockArnoldi1.getEvals(); 
  double * evali = MyBlockArnoldi1.getiEvals();

  // Output results to screen
  MyBlockArnoldi1.currentStatus();

  /* ********************************************************************** */
  /* phase 2: estimate the null space                                       */
  /* - scale diagonal                                                       */
  /* - compute lambda max of I - D^{-1} A /lamdba_max                       */
  /* ********************************************************************** */

  double InvLambdaMax = 1.0/(1.5*sqrt(pow(evalr[0],2) + pow(evali[0],2)));
  
  for( int i=0 ; i<NumMyElements ; ++i )
    (*RowMatrixDiagonal)[i] *= InvLambdaMax;
  
  // Call the ctor that calls the petra ctor for a matrix
  MLMat<double> Amat2(*RowMatrix,*RowMatrixDiagonal,false);	
  Anasazi::Eigenproblem<double> MyProblem2(&Amat2, &ivec);

  // Initialize the Block Arnoldi solver
  Anasazi::BlockArnoldi<double> MyBlockArnoldi2(MyProblem2, tol, nev, length, block, 
					       which, step, restarts);
	
  // Inform the solver that the problem is symmetric
  //MyBlockArnoldi2.setSymmetric(true);
  MyBlockArnoldi2.setDebugLevel(0);

  MyBlockArnoldi2.solve();

  evalr = MyBlockArnoldi2.getEvals(); 
  evali = MyBlockArnoldi2.getiEvals();

  for( int i=0 ; i<nev ; ++i ) {
    cout << "eval[" << i << "] = " << evalr[i] << " + i " << evali[i] << endl;
  }  
  
  // Retrieve eigenvectors
  Anasazi::PetraVec<double> evecr(RowMatrix->RowMatrixRowMap(), nev);
  MyBlockArnoldi2.getEvecs( evecr );
  Anasazi::PetraVec<double> eveci(RowMatrix->RowMatrixRowMap(), nev);
  MyBlockArnoldi2.getiEvecs( eveci );

  /* ********************************************************************** */
  /* Scale vectors (real and imag part) so that max is 1 for each of them   */
  /* ********************************************************************** */

  double MaxRealValue[nev], MaxImagValue[nev];
  evecr.MaxValue(MaxRealValue);
  eveci.MaxValue(MaxImagValue);

  for( int i=0 ; i<nev ; ++i ) {
    if( MaxRealValue[i] != 0.0 ) evecr(i)->Scale(1.0/MaxRealValue[i]);
    if( MaxImagValue[i] != 0.0 ) eveci(i)->Scale(1.0/MaxImagValue[i]);
  }

  // Output results to screen
  MyBlockArnoldi2.currentStatus();

  delete RowMatrixDiagonal;

  NullSpacePtr = new double[NumMyElements * nev];

  for( int j=0 ; j<nev ; ++j ) 
    for( int i=0 ; i<NumMyElements ; ++i ) 
      NullSpacePtr[j*NumMyElements+i] = evecr[j][i];
  
  return 0;
  
}

#endif
