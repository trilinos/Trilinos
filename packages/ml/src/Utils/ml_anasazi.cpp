#include "ml_common.h"
#include "ml_include.h"

#if defined(ML_WITH_EPETRA) && defined(HAVE_ML_ANASAZI)

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

#include "AnasaziPetraInterface.hpp"
#include "AnasaziBlockArnoldi.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziReturnType.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziMatrix.hpp"
#include "AnasaziOperator.hpp"

#include "ml_epetra_utils.h"
#include "Epetra_CrsMatrix.h"
#include "ml_anasazi.h"

using namespace Anasazi;

// ================================================ ====== ==== ==== == =

enum MLMatOp { A_MATRIX, I_MINUS_A_MATRIX, A_PLUS_AT_MATRIX, A_MINUS_AT_MATRIX };

// ================================================ ====== ==== ==== == =

template <class TYPE> 
class MLMat : public virtual Matrix<TYPE> {
public:
  MLMat(const Epetra_RowMatrix &, const MLMatOp MatOp, const bool UseDiagScaling );
  ~MLMat();
  ReturnType ApplyMatrix ( const MultiVec<TYPE>& x, 
			   MultiVec<TYPE>& y ) const;
  void SetMatOp(const MLMatOp MatOp) 
  {
    MatOp_ = MatOp;
  }

  void SetDiagScaling(const bool UseDiagScaling) 
  {
    UseDiagScaling_ = UseDiagScaling;
  }
  
private:
  const Epetra_RowMatrix & Mat_;
  MLMatOp MatOp_;
  const bool UseDiagScaling_;
  double Scale_;
  Epetra_MultiVector * tmp_;
  Epetra_Vector * Diagonal_;
  
};

// ================================================ ====== ==== ==== == =

template <class TYPE>
MLMat<TYPE>::MLMat(const Epetra_RowMatrix & Matrix,
		   const MLMatOp MatOp, const bool UseDiagScaling ) 
  : Mat_(Matrix),
    MatOp_(MatOp),
    UseDiagScaling_(UseDiagScaling),
    Scale_(1.0),
    tmp_(0)
{
  if( UseDiagScaling_ ) {
    
    Diagonal_ = new Epetra_Vector(Matrix.RowMatrixRowMap());
    Matrix.ExtractDiagonalCopy(*Diagonal_);
    
    int NumMyElements = Matrix.RowMatrixRowMap().NumMyElements();
    for( int i=0 ; i<NumMyElements ; ++i ) {
      if( (*Diagonal_)[i] != 0.0 ) (*Diagonal_)[i] = 1.0/(*Diagonal_)[i];
      else  (*Diagonal_)[i] = 0.0;
    }

  }
  
}

// ================================================ ====== ==== ==== == =

template <class TYPE>
MLMat<TYPE>::~MLMat() 
{
  if( UseDiagScaling_ ) delete Diagonal_;
  if( tmp_ )  delete tmp_;
}

// ================================================ ====== ==== ==== == =

template <class TYPE>
ReturnType MLMat<TYPE>::ApplyMatrix ( const MultiVec<TYPE>& x, 
				      MultiVec<TYPE>& y ) const 
{

  int info = 0;
  MultiVec<TYPE> & temp_x = const_cast<MultiVec<TYPE> &>(x);
  Epetra_MultiVector * vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  Epetra_MultiVector * vec_y = dynamic_cast<Epetra_MultiVector* >(&y);

  if( vec_x==0 || vec_y==0 ) return Failed;

  info = const_cast<Epetra_RowMatrix &>(Mat_).Apply( *vec_x, *vec_y );
  vec_y->Scale(Scale_);

  if( info ) return Failed;

  if( tmp_ == 0 ) {
    const_cast<MLMat<TYPE> *>(this)->tmp_ = new Epetra_MultiVector(*vec_x);
  }
  
  if( MatOp_ == A_PLUS_AT_MATRIX || MatOp_ == A_MINUS_AT_MATRIX ) {
    const Epetra_RowMatrix & RowMatrix = dynamic_cast<const Epetra_RowMatrix &>(Mat_);
    info = RowMatrix.Multiply(true,*vec_x,*tmp_);
    tmp_->Scale(Scale_);
    if( info ) return Failed;
  }
  if( MatOp_ == A_PLUS_AT_MATRIX ) vec_y->Update(1.0,*tmp_,1.0);
  else if( MatOp_ == A_MINUS_AT_MATRIX ) vec_y->Update(-1.0,*tmp_,1.0);
  
  // diagonal scaling
  if( UseDiagScaling_ == true ) {
    for( int j=0 ; j<vec_y->NumVectors() ; ++j ) {
      for( int i=0 ; i<vec_y->Map().NumMyElements() ; ++i ) (*vec_y)[j][i] *= (*Diagonal_)[i];
    }    
  }
  
  if( MatOp_ == I_MINUS_A_MATRIX ) {
    vec_y->Update (1.0,*vec_x,-1.0);
  }

  return Ok; 
}

#ifdef HAVE_ML_TEUCHOS
#include "Teuchos_ParameterList.hpp"

namespace ML_Anasazi {

// ================================================ ====== ==== ==== == =

int Interface(const Epetra_RowMatrix * RowMatrix, Epetra_MultiVector & EigenVectors,
		 	  double RealEigenvalues[], double ImagEigenvalues[],
			  Teuchos::ParameterList & List) 
{

  int MyPID = RowMatrix->Comm().MyPID();
  
  /* ********************************************************************** */
  /* Retrive parameters' choices                                            */
  /* ********************************************************************** */

  MLMatOp MatOp;
  string MatOpStr = "A";
  MatOpStr = List.get("eigen-analysis: matrix operation", MatOpStr);
  if( MatOpStr == "A" ) MatOp = A_MATRIX;
  else if( MatOpStr == "I-A" ) MatOp = I_MINUS_A_MATRIX;
  else if( MatOpStr == "A+A^T" ) MatOp = A_PLUS_AT_MATRIX;
  else if( MatOpStr == "A-A^T" ) MatOp = A_MINUS_AT_MATRIX;

  bool UseDiagScaling = true;
  UseDiagScaling = List.get("eigen-analysis: use diagonal scaling", UseDiagScaling);

  int length = List.get("eigen-analysis: length", 20);
  double tol = List.get("eigen-analysis: tolerance", 1.0e-5);
  string which = List.get("eigen-analysis: action", "LM");
  int restarts = List.get("eigen-analysis: restart", 100);
  bool isSymmetric = List.get("eigen-analysis: symmetric problem", false);

  int output = List.get("eigen-analysis: output", 5);
  
  if( output > 5 && MyPID == 0 ) {
    if( MatOp == A_MATRIX ) cout << "ML_Anasazi : Computing eigenvalues of A" << endl;
    if( MatOp == I_MINUS_A_MATRIX ) cout << "ML_Anasazi : Computing eigenvalues of I - A" << endl;
    if( MatOp == A_PLUS_AT_MATRIX ) cout << "ML_Anasazi : Computing eigenvalues of A + A^T" << endl;
    if( MatOp == A_MINUS_AT_MATRIX ) cout << "ML_Anasazi : Computing eigenvalues of A - A^T" << endl;
    if( UseDiagScaling ) cout << "ML_Anasazi : where A is scaled by D^{-1}" << endl;
    if( isSymmetric ) cout << "ML_Anasazi : Problem is symmetric" << endl;
    cout << "ML_Anasazi : Tolerance = " << tol << endl;
    cout << "ML_Anasazi : Required Action = " << which << endl;
	
  }

  bool PrintCurrentStatus =  List.get("eigen-analysis: print current status", false);

  /* ********************************************************************** */
  /* view mode for PetraVectors                                             */
  /* ********************************************************************** */
  
  int NumBlocks = EigenVectors.NumVectors();

  int * indices = new int[NumBlocks];
  for( int i=0 ; i<NumBlocks ; ++i ) indices[i] = i;
  
  Anasazi::PetraVec<double> Vectors(View,EigenVectors,indices,NumBlocks);
  
  delete [] indices;

  /* ********************************************************************** */
  /* Perform required action                                                */
  /* ********************************************************************** */

  int step = restarts*length*NumBlocks;
  
  // Call the ctor that calls the petra ctor for a matrix
  MLMat<double> Amat(*RowMatrix,MatOp,UseDiagScaling);
  
  Anasazi::Eigenproblem<double> MyProblem(&Amat, &Vectors);

  // Initialize the Block Arnoldi solver
  Anasazi::BlockArnoldi<double> MyBlockArnoldi1(MyProblem, tol, NumBlocks, length, NumBlocks, 
						which, step, restarts);
	
  // Inform the solver that the problem is symmetric
  if( isSymmetric ) MyBlockArnoldi1.setSymmetric(true);
  MyBlockArnoldi1.setDebugLevel(0);

  //  MyBlockArnoldi.iterate(5);
  
  // Solve the problem to the specified tolerances or length
  MyBlockArnoldi1.solve();
  
  // Obtain results directly

  double * evalr = MyBlockArnoldi1.getEvals(); 
  double * evali = MyBlockArnoldi1.getiEvals();

  for( int i=0 ; i<NumBlocks ; ++i ) {
    RealEigenvalues[i] = evalr[i];
    ImagEigenvalues[i] = evali[i];
  }

  delete [] evalr;
  delete [] evali;
  
  // Output results to screen
  if( PrintCurrentStatus > 5 && MyPID == 0 ) MyBlockArnoldi1.currentStatus();

  MyBlockArnoldi1.getEvecs(Vectors);
  
  /* ********************************************************************** */
  /* Scale vectors (real and imag part) so that max is 1 for each of them   */
  /* ********************************************************************** */

  if( List.get("eigen-analysis: normalize eigenvectors",false) ) {
    
    double * MaxVals = new double[NumBlocks];
    double * MinVals = new double[NumBlocks];
    double * Vals = new double[NumBlocks];
    EigenVectors.MaxValue(MaxVals);
    EigenVectors.MinValue(MinVals);

    for( int i=0 ; i<NumBlocks ; ++i ) {
      Vals[i] = EPETRA_MAX(abs(MaxVals[i]),abs(MinVals[i]));
      if( Vals[i] == abs(MaxVals[i]) && MaxVals[i]<0 ) Vals[i] *= -1.0;
      if( Vals[i] == abs(MinVals[i]) && MinVals[i]<0 ) Vals[i] *= -1.0;
      Vals[i] = 1.0/Vals[i];
      EigenVectors(i)->Scale(Vals[i]);

    }

    delete [] MaxVals;
    delete [] MinVals;
    delete [] Vals;
  }

  return 0;
  
}

// ================================================ ====== ==== ==== == =

int GetFieldOfValuesBox(const Epetra_RowMatrix * RowMatrix, 
				    double & MaxReal, double & MaxImag,
				    Teuchos::ParameterList & List ) 
{

  Epetra_Time Time(RowMatrix->Comm());
  
  int MyPID = RowMatrix->Comm().MyPID();

  bool UseDiagScaling = true;
  UseDiagScaling = List.get("field-of-values: use diagonal scaling", UseDiagScaling);
  
  int length = List.get("field-of-values: length", 20);
  double tol = List.get("field-of-values: tolerance", 1.0e-5);
  int restarts = List.get("field-of-values: restart", 100);

  int output = List.get("output", 5);
  
  if( output > 5 && MyPID == 0 ) {
    cout << "ML_Anasazi : Estimate box containing the field of values" << endl;
    cout << "ML_Anasazi : Tolerance = " << tol << endl;	
    cout << "ML_Anasazi : Computing eigenvalues of A + A^T" << endl;
    if( UseDiagScaling ) cout << "ML_Anasazi : where A is scaled by D^{-1}" << endl;
  }

  bool PrintCurrentStatus =  List.get("field-of-values: print current status", false);

  bool isSymmetric = List.get("eigen-analysis: symmetric problem", false);
  
  /* ********************************************************************** */
  /* First compute A + A^T to get the real bound                            */
  /* ********************************************************************** */
  
  Anasazi::PetraVec<double> Vectors(RowMatrix->RowMatrixRowMap(), 1);
  Vectors.MvRandom();
  
  int step = restarts*length*1;
  
  // Call the ctor that calls the petra ctor for a matrix
  MLMat<double> Amat(*RowMatrix,A_PLUS_AT_MATRIX,UseDiagScaling);
  
  Anasazi::Eigenproblem<double> MyProblem(&Amat, &Vectors);

  // Initialize the Block Arnoldi solver
  Anasazi::BlockArnoldi<double> MyBlockArnoldi1(MyProblem, tol, 1, length, 1, 
						"LM", step, restarts);

  // Inform the solver that the problem is symmetric
  MyBlockArnoldi1.setSymmetric(true);
  MyBlockArnoldi1.setDebugLevel(0);

  //  MyBlockArnoldi.iterate(5);
  
  // Solve the problem to the specified tolerances or length
  MyBlockArnoldi1.solve();
  
  // Obtain results directly

  double * evalr = MyBlockArnoldi1.getEvals(); 

  MaxReal = evalr[0] / 2;

  double * residuals  = MyBlockArnoldi1.getResiduals();
  if( output > 5 && MyPID == 0 ) {
    cout << "ML_Anasazi : Ritz Residual for A^T + A = " << residuals[0] << endl;
  }

  delete [] evalr;
  delete [] residuals;  
  
  if( PrintCurrentStatus && MyPID == 0 ) MyBlockArnoldi1.currentStatus();

  /* ********************************************************************** */
  /* First compute A - A^T to get the real bound                            */
  /* ********************************************************************** */
  
  if( isSymmetric == false ) {
    
    if( output > 5 && MyPID == 0 ) {
      cout << "ML_Anasazi : Computing eigenvalues of A - A^T" << endl;
    }
    
    Vectors.MvRandom();
  
    // Call the ctor that calls the petra ctor for a matrix
    MLMat<double> Amat2(*RowMatrix,A_MINUS_AT_MATRIX,UseDiagScaling);
    
    Anasazi::Eigenproblem<double> MyProblem2(&Amat2, &Vectors);
    
    // Initialize the Block Arnoldi solver
    Anasazi::BlockArnoldi<double> MyBlockArnoldi2(MyProblem2, tol, 1, length, 1, 
						  "LM", step, restarts);
    
    // Inform the solver that the problem is symmetric
    MyBlockArnoldi2.setSymmetric(false);
    MyBlockArnoldi2.setDebugLevel(0);
    
    MyBlockArnoldi2.iterate(5);
    
    // Solve the problem to the specified tolerances or length
    MyBlockArnoldi2.solve();
    
    // Obtain results directly
    
    double * evali = MyBlockArnoldi2.getiEvals();
    
    MaxImag = abs(evali[0] / 2);

    delete [] evali;
    
    residuals  = MyBlockArnoldi2.getResiduals();
    if( output > 5 && MyPID == 0 ) {
      cout << "ML_Anasazi : Ritz Residual for A^T - A = " << residuals[0] << endl;
    }
    
    if( PrintCurrentStatus && MyPID == 0 ) MyBlockArnoldi2.currentStatus();

    if( output > 5 && MyPID == 0 ) {
      cout << "ML_Anasazi : Time = " << Time.ElapsedTime() << " (s)" << endl;
    }
    
  } else {

    MaxImag = 0;
    
  }
  
  return 0;
  
}

} // namespace ML_Anasazi

// ================================================ ====== ==== ==== == =

extern "C" {
  
int ML_Anasazi_Get_FiledOfValuesBox_Interface(ML_Operator * Amat,
					      struct ML_Field_Of_Values * fov )
{

  Epetra_CrsMatrix * CrsTemp;
  int MaxNumNonzeros;
  double CPUTime;
  
  ML_Operator2EpetraCrsMatrix(Amat,CrsTemp,MaxNumNonzeros,
			      true,CPUTime);
  
  double MaxReal,MaxImag;
  Teuchos::ParameterList * EigenList = (Teuchos::ParameterList *) fov->EigenList;
  
  ML_Anasazi::GetFieldOfValuesBox(CrsTemp,MaxReal,MaxImag,*EigenList);

  double eta = MaxImag/MaxReal;

  fov->eta = eta;
  fov->real_max = MaxReal;
  fov->imag_max = MaxImag;

  delete CrsTemp;
  
  return 0;
 
}
}

#endif /* ifdef HAVE_ML_TEUCHOS */

extern "C" {  
// ================================================ ====== ==== ==== == =

int ML_Anasazi_Get_SpectralNorm_Anasazi(ML_Operator * Amat,
					int MaxIters, double Tolerance,
					int IsProblemSymmetric,
					int UseDiagonalScaling,
					double * LambdaMax )
{

  Epetra_CrsMatrix * RowMatrix;
  int MaxNumNonzeros;
  double CPUTime;
  
  ML_Operator2EpetraCrsMatrix(Amat,RowMatrix,MaxNumNonzeros,
			      true,CPUTime);
  
  bool verbose( RowMatrix->Comm().MyPID() == 0 && (5 < ML_Get_PrintLevel() ) );
  
  int length = MaxIters;
  int restarts = 1;

  if( verbose ) {
    cout << "ML_Anasazi : Estimate Lambda Max, ";
    if( IsProblemSymmetric == ML_TRUE ) cout << "problem is symmetric, ";
    if( UseDiagonalScaling == ML_TRUE ) cout << "diagonal scaling.";
    cout << endl;
  }

  /* ********************************************************************** */
  /* First compute A + A^T to get the real bound                            */
  /* ********************************************************************** */
  
  Anasazi::PetraVec<double> Vectors(RowMatrix->RowMatrixRowMap(), 1);
  Vectors.MvRandom();
  
  int step = restarts*length*1;
  
  // Call the ctor that calls the petra ctor for a matrix

  bool flag = false;
  if( UseDiagonalScaling == ML_TRUE ) flag = true;
  
  MLMat<double> MLMatrix(*RowMatrix,A_MATRIX,flag);
  
  Anasazi::Eigenproblem<double> MyProblem(&MLMatrix, &Vectors);

  // Initialize the Block Arnoldi solver
  Anasazi::BlockArnoldi<double> MyBlockArnoldi(MyProblem, Tolerance, 1, length, 1, 
						"LM", step, restarts);

  // Inform the solver that the problem is symmetric
  if( IsProblemSymmetric == ML_TRUE ) MyBlockArnoldi.setSymmetric(true);
  else                                MyBlockArnoldi.setSymmetric(false);
  
  MyBlockArnoldi.setDebugLevel(0);

  //  MyBlockArnoldi.iterate(5);
  
  // Solve the problem to the specified tolerances or length
  MyBlockArnoldi.solve();
  
  // Obtain results direc3tly

  double * evalr = MyBlockArnoldi.getEvals();
  double * evali = MyBlockArnoldi.getiEvals(); 

  *LambdaMax = sqrt(pow(evalr[0],2) + pow(evali[0],2));
  
  if( verbose ) {
    double * residuals  = MyBlockArnoldi.getResiduals();
    cout << "ML_Anasazi : Ritz Residual = " << residuals[0] << endl;
    delete [] residuals;
   
  }
  
  delete RowMatrix;

  delete [] evalr;
  delete [] evali;
  
  return 0;
 
}
}
#else

#include <stdio.h>

extern "C" {  
  // ================================================ ====== ==== ==== == =
  
  int ML_Anasazi_Get_SpectralNorm_Anasazi(ML_Operator * Amat,
					int MaxIters, double Tolerance,
					  int IsProblemSymmetric,
					  int UseDiagonalScaling,
					double * LambdaMax )
  {
    
    puts("You must configure with options --with-ml_epetra and ");
    puts("--with-ml_anasazi to estimate lambda max with Anasazi.");

    exit( EXIT_FAILURE );

    return( -1 );
  
  }
  
}

#endif
