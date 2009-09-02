/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_common.h"
#include "ml_include.h"
#include "ml_smoother.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_ANASAxI) && defined(HAVE_ML_TEUCHOS)

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

#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziOutputManager.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziOperator.hpp"

#include "ml_epetra_utils.h"
#include "Epetra_CrsMatrix.h"
#include "ml_anasazi.h"
#include "Teuchos_ParameterList.hpp"

void Orthogonalize(Anasazi::EpetraMultiVec & vec)
{

  int NumVectors = vec.NumVectors();

  // replace each vector v_i by v_i - rho_ij v_j,
  // where rho_ij = <v_i, v_j>, for all j <i
  for( int i=1 ; i<NumVectors ; ++i ) {
    double norm;
    // scale previous (possibly modified at previous step)
    vec(i-1)->Norm2(&norm);
    if( norm < 1e-10 ) vec(i-1)->Scale(0.0);
    else               vec(i-1)->Scale(1.0/norm);
    // scale current (still unscaled)
    vec(i)->Norm2(&norm);
    if( norm < 1e-10 ) vec(i)->Scale(0.0);
    else               vec(i)->Scale(1.0/norm);
    // orthogonalize, only `i' is modified
    for( int j=0 ; j<i ; ++j ) {
      double rho;
      vec(i)->Dot(*(vec(j)),&rho);
      vec(i)->Update(-rho,*(vec(j)),1.0);
      vec(i)->Norm2(&norm);
      if( norm < 1e-10 ) vec(i)->Scale(0.0);
      else               vec(i)->Scale(1.0/norm);      
    }
  }

  return;
}

// ================================================ ====== ==== ==== == =

using namespace Anasazi;

// ================================================ ====== ==== ==== == =

enum MLMatOp { A_MATRIX, I_MINUS_A_MATRIX, A_PLUS_AT_MATRIX,
	       A_MINUS_AT_MATRIX, ITERATION_MATRIX,
               PREC_MATRIX, SMOOTHED_MATRIX };


// ================================================ ====== ==== ==== == =
//! ML/Anasazi matrix for matrix-vector product.
/*! Anasazi requires the input matrix to be derived from class
 * Anasazi::Operator<TYPE> (here TYPE is double).
 *
 * This class wraps a given Epetra_RowMatrix into an Anasazi matrix. As
 * Epetra_RowMatrix can perform a multiplication by A, as well as by A^T, the
 * class can perform the following operations:
 * - A * x
 * - (A + A^T) * x
 * - (A - A^T) * x
 * - (I - ML^{-1} A ) * x
 * - (ML^{-1} A) * x
 * 
 * Optionally, the matrix can be scaled by the diagonal.
 *
 */

class MLMat : public virtual Operator<double> {

public:
  //! Default constructor
  /*! Constructs an ML/Anasazi matrix.
    \param Matrix (In) : square matrix
    \param MatOp (In) : MatOp is an enum MLMatOp, that can assume the
      following values:
      - A_MATRIX (compute the eigenvalues of A)
      - I_MINUS_A_MATRIX (compute the eigenvalues of I-A
      - I_PLUS_A_MATRIX (compute the eigenvalue of I+A)
      - A_MINUS_A_MATRIX (compute the eigenvalue of A-A^T)
      - A_PLUS_A_MATRIX (compute the eigenvalue of A+A^T)
      - ITERATION_MATRIX (compute the eigenvalue of I-ML^{-1}A, 
        where ML is an already build ML preconditioner
      - PREC_MATRIX (compute the eigenvalues of ML^{-1}A).	
      - SMOOTHED_MATRIX (compute the eigenvalues of S^{-1}A).	
    \param UseDiagScaling (In) : if \c true, the matrix is scaled by the
      diagonal
    \param ml (In) : pointer to an already built ML hierarchy.
   */
  MLMat(const Epetra_RowMatrix & Matrix, const MLMatOp MatOp, 
	const bool UseDiagScaling,
	ML* ml = 0, ML_Smoother* smoother = 0);
  ~MLMat();

  //! Applies the flavor of ML/Anasazi matrix to Anasazi::Vector x, stores the result in y.
  ReturnType Apply ( const MultiVec<double>& x, 
		     MultiVec<double>& y ) const;

  //! Sets matrix operation.
  void SetMatOp(const MLMatOp MatOp) 
  {
    MatOp_ = MatOp;
  }

  //! Sets the diagonal scaling.
  // NOTE:  This function should not be used if UseDiagScaling_ is const.
/*  void SetDiagScaling(const bool UseDiagScaling) 
  {
    UseDiagScaling_ = UseDiagScaling;
  }
*/  
private:
  const Epetra_RowMatrix& Mat_;
  MLMatOp MatOp_;
  const bool UseDiagScaling_;
  double Scale_;
  Epetra_MultiVector* tmp_;
  Epetra_Vector* Diagonal_;

  Epetra_Vector* Mask_;
  int NumMyRows_;
  ML* ml_;
  ML_Smoother* smoother_;
};


MLMat::MLMat(const Epetra_RowMatrix & Matrix,
		   const MLMatOp MatOp, const bool UseDiagScaling,
		   ML* ml, ML_Smoother* smoother) 
  : Mat_(Matrix),
    MatOp_(MatOp),
    UseDiagScaling_(UseDiagScaling),
    Scale_(1.0),
    tmp_(0),
    ml_(ml),
    smoother_(smoother)
{
  NumMyRows_ = Matrix.RowMatrixRowMap().NumMyElements();
  
  if( UseDiagScaling_ ) {
    
    Diagonal_ = new Epetra_Vector(Matrix.RowMatrixRowMap());
    Matrix.ExtractDiagonalCopy(*Diagonal_);
    
    int NumMyElements = Matrix.RowMatrixRowMap().NumMyElements();
    for( int i=0 ; i<NumMyElements ; ++i ) {
      if( (*Diagonal_)[i] != 0.0 ) (*Diagonal_)[i] = 1.0/(*Diagonal_)[i];
      else  (*Diagonal_)[i] = 0.0;
    }

  }

  // define the Dirichlet boundary nodes
  Mask_ = new Epetra_Vector(Matrix.RowMatrixRowMap());

  for( int i=0 ; i<Matrix.NumMyRows() ; ++i ) {
    int Nnz;
    Matrix.NumMyRowEntries(i,Nnz);
    if( Nnz <= 1 ) (*Mask_)[i] = 0.0;
    else           (*Mask_)[i] = 1.0;
  }
  
}

// ================================================ ====== ==== ==== == =

MLMat::~MLMat() 
{
  if( UseDiagScaling_ ) delete Diagonal_;
  if( tmp_ )  delete tmp_;
  if( Mask_ ) delete Mask_;
}

// ================================================ ====== ==== ==== == =

ReturnType MLMat::Apply ( const MultiVec<double>& x, 
				MultiVec<double>& y ) const 
{

  int info = 0;
  MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
  Epetra_MultiVector * vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  Epetra_MultiVector * vec_y = dynamic_cast<Epetra_MultiVector* >(&y);

  if( vec_x==0 || vec_y==0 ) return Failed;

  int NumVectors = vec_x->NumVectors();
  double ** tmp_view, ** x_view, ** y_view, * mask_view;
  
  if( tmp_ == 0 ) {
    const_cast<MLMat *>(this)->tmp_ = new Epetra_MultiVector(*vec_x);
  }
  
  // extract views for x and tmp_
  vec_x->ExtractView(&x_view);
  vec_y->ExtractView(&y_view);
  tmp_->ExtractView(&tmp_view);
  Mask_->ExtractView(&mask_view);

  if( MatOp_ == PREC_MATRIX ) {

    assert (ml_ != 0);
    
    // Here I apply ML^{-1} A
    
    const Epetra_RowMatrix & RowMatrix = dynamic_cast<const Epetra_RowMatrix &>(Mat_);

    // 0- zero out tmp_
    tmp_->PutScalar(0.0);

    // 1- apply the linear system matrix to vec_x
    info = RowMatrix.Multiply(false,*vec_x,*tmp_);
    if( info ) return Failed;

    // 1.1- damp out Dirichlet nodes (FIXME: am I useful?)
    for( int j=0 ; j<NumVectors ; ++j ) {
      for( int i=0 ; i<NumMyRows_ ; ++i ) {
	tmp_view[j][i] *= mask_view[i];
      }    
    }

    // 2- apply the multilevel hierarchy
    vec_y->PutScalar(0.0);
    for (int i = 0 ; i < NumVectors ; ++i) 
      ML_Solve_MGV(ml_,tmp_view[i], y_view[i]);
    
    // 4- return and skip the crap below
    return Ok;
    
  }
  else if( MatOp_ == SMOOTHED_MATRIX ) {

    assert (smoother_ != 0);
    
    // Here I apply S^{-1} A
    const Epetra_RowMatrix& RowMatrix = dynamic_cast<const Epetra_RowMatrix&>(Mat_);

    // 1- apply the linear system matrix to vec_x
    info = RowMatrix.Multiply(false,*vec_x,*tmp_);
    if (info) return Failed;

    // 2- apply the smoother
    vec_y->PutScalar(0.0);
    for (int i = 0 ; i < NumVectors ; ++i) 
      ML_Smoother_Apply(smoother_,NumMyRows_,y_view[i],
                        NumMyRows_,tmp_view[i],ML_NONZERO);

    // 4- return and skip the crap below
    return Ok;
    
  }

  if( MatOp_ == ITERATION_MATRIX ) {

    assert( ml_ != 0 );
    
    // Here I apply I - ML^{-1} A
    
    const Epetra_RowMatrix & RowMatrix = dynamic_cast<const Epetra_RowMatrix &>(Mat_);

    // 0- zero out tmp_
    tmp_->PutScalar(0.0);

    // 1- apply the linear system matrix to vec_x
    info = RowMatrix.Multiply(false,*vec_x,*tmp_);
    if( info ) return Failed;

    // 1.1- damp out Dirichlet nodes (FIXME: am I useful?)
    for( int j=0 ; j<NumVectors ; ++j ) {
      for( int i=0 ; i<NumMyRows_ ; ++i ) {
	tmp_view[j][i] *= mask_view[i];
      }    
    }

    // 2- apply the multilevel hierarchy
    vec_y->PutScalar(0.0);
    for( int i=0 ; i<NumVectors ; ++i ) ML_Solve_MGV(ml_,tmp_view[i], y_view[i]);
    
    // 3- add the contribution of the identity matrix
    vec_y->Update (1.0,*vec_x,-1.0);

    // 4- return and skip the crap below
    return Ok;
    
  }

  for( int j=0 ; j<NumVectors ; ++j ) {
    for( int i=0 ; i<NumMyRows_ ; ++i ) {
      tmp_view[j][i] = x_view[j][i] * mask_view[i];
    }    
  }
  
  info = const_cast<Epetra_RowMatrix &>(Mat_).Apply( *tmp_, *vec_y );
  vec_y->Scale(Scale_);
  if( info ) return Failed;
  
  if( MatOp_ == A_PLUS_AT_MATRIX || MatOp_ == A_MINUS_AT_MATRIX ) {
    const Epetra_RowMatrix & RowMatrix = dynamic_cast<const Epetra_RowMatrix &>(Mat_);

    info = RowMatrix.Multiply(true,*vec_x,*tmp_);
    if( info ) return Failed;

    for( int j=0 ; j<NumVectors ; ++j ) {
      for( int i=0 ; i<NumMyRows_ ; ++i ) {
	tmp_view[j][i] *= mask_view[i];
      }    
    }
    
    tmp_->Scale(Scale_);
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

namespace ML_Anasazi {

// ================================================ ====== ==== ==== == =

int Interface(const Epetra_RowMatrix * RowMatrix, Epetra_MultiVector & EigenVectors,
	      double RealEigenvalues[], double ImagEigenvalues[],
	      Teuchos::ParameterList & List,
	      double RealEigenvectors[],  double ImagEigenvectors[],
	      int * NumRealEigenvectors, int *  NumImagEigenvectors,
	      ML * ml) 
{

  int MyPID = RowMatrix->Comm().MyPID();
  
  /* ********************************************************************** */
  /* Retrive parameters' choices                                            */
  /* ********************************************************************** */

  bool UseDiagScaling = true;
  UseDiagScaling = List.get("eigen-analysis: use diagonal scaling", UseDiagScaling);
  bool isSymmetric = List.get("eigen-analysis: symmetric problem", false);

  MLMatOp MatOp = A_MATRIX;
  string MatOpStr = "A";
  MatOpStr = List.get("eigen-analysis: matrix operation", MatOpStr);

  if( MatOpStr == "A" ) MatOp = A_MATRIX;
  else if( MatOpStr == "I-A" ) MatOp = I_MINUS_A_MATRIX;
  else if( MatOpStr == "A+A^T" ) MatOp = A_PLUS_AT_MATRIX;
  else if( MatOpStr == "A-A^T" ) MatOp = A_MINUS_AT_MATRIX;
  else if( MatOpStr == "ML^{-1}A" ) MatOp = PREC_MATRIX;
  else if( MatOpStr == "S^{-1}A" ) MatOp = SMOOTHED_MATRIX;
  else if( MatOpStr == "I-ML^{-1}A" ) {
    MatOp = ITERATION_MATRIX;
    UseDiagScaling = false; // doesn't work with iteration matrix
    // note that `which' is usually "LM"
    // also isSymmetric is probably false, set it to `false' for safety
    isSymmetric = false;
  }
  
  int length = List.get("eigen-analysis: length", 20);
  int BlockSize = List.get("eigen-analysis: block-size", 1);
  double tol = List.get("eigen-analysis: tolerance", 1.0e-5);
  string which = List.get("eigen-analysis: action", "LM");
  int restarts = List.get("eigen-analysis: restart", 100);

  int output = List.get("eigen-analysis: output", 6);
  
  if( output > 5 && MyPID == 0 ) {
    if( MatOp == A_MATRIX ) cout << "ML_Anasazi : Computing eigenvalues of A" << endl;
    if( MatOp == I_MINUS_A_MATRIX ) cout << "ML_Anasazi : Computing eigenvalues of I - A" << endl;
    if( MatOp == A_PLUS_AT_MATRIX ) cout << "ML_Anasazi : Computing eigenvalues of A + A^T" << endl;
    if( MatOp == A_MINUS_AT_MATRIX ) cout << "ML_Anasazi : Computing eigenvalues of A - A^T" << endl;
    if( MatOp == ITERATION_MATRIX ) cout << "ML_Anasazi : Computing eigenvalues of I-BA" << endl;
    if( UseDiagScaling ) cout << "ML_Anasazi : where A is scaled by D^{-1}" << endl;
    if( isSymmetric ) cout << "ML_Anasazi : Problem is symmetric" << endl;
    cout << "ML_Anasazi : Tolerance = " << tol << endl;
    cout << "ML_Anasazi : Required Action = " << which << endl;
	
  }

  bool PrintCurrentStatus =  List.get("eigen-analysis: print current status", false);

  /* ********************************************************************** */
  /* view mode for EpetraMultiVectors                                             */
  /* ********************************************************************** */
  
  int NumBlocks = EigenVectors.NumVectors();

  std::vector<int> indices( NumBlocks );
  for( int i=0 ; i<NumBlocks ; ++i ) indices[i] = i;
  
  Teuchos::RefCountPtr<Anasazi::EpetraMultiVec> Vectors =
    Teuchos::rcp( new Anasazi::EpetraMultiVec(View,EigenVectors,indices) );
  Vectors->MvRandom();
  
  /* ********************************************************************** */
  /* Perform required action                                                */
  /* ********************************************************************** */

  // Call the ctor that calls the petra ctor for a matrix
  Teuchos::RefCountPtr<MLMat> Amat =
    Teuchos::rcp( new MLMat(*RowMatrix,MatOp,UseDiagScaling,ml) );
  
  typedef Anasazi::MultiVec<double> MV;
  typedef Anasazi::Operator<double> OP;

  Teuchos::ParameterList AnasaziPL;
  AnasaziPL.set("Block Size", BlockSize);
  AnasaziPL.set("Max Blocks", length);
  AnasaziPL.set("Max Restarts", restarts);
  AnasaziPL.set("Tol", tol);

  // Create default output manager 
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM = 
    Teuchos::rcp( new Anasazi::OutputManager<double>() );

  // Create a sorting manager to handle the sorting of eigenvalues in the solver
  Teuchos::RefCountPtr<Anasazi::BasicSort<double, MV, OP> > MySort = 
    Teuchos::rcp( new Anasazi::BasicSort<double, MV, OP>(which) );
  
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(Amat, Vectors) );

  // Inform the solver that the problem is symmetric
  if( isSymmetric ) MyProblem->SetSymmetric(true);
  else              MyProblem->SetSymmetric(false);

  // Set the number of eigenvalues required
  MyProblem->SetNEV( NumBlocks );
  MyProblem->SetProblem();

  // Initialize the Block Arnoldi solver
  Anasazi::BlockKrylovSchur<double, MV, OP> MyBlockKrylovSchur1(MyProblem, MySort, MyOM, AnasaziPL);

  // Solve the problem to the specified tolerances or length
  MyBlockKrylovSchur1.solve();
  
  // Obtain results directly for eigenvalues
  // The vector is 2*nev in length if problem is non-symmetric
  Teuchos::RefCountPtr<std::vector<double> > evals = MyProblem->GetEvals();

  for( int i=0 ; i<NumBlocks ; ++i ) {
    RealEigenvalues[i] = (*evals)[i];
    if( isSymmetric )
      ImagEigenvalues[i] = 0.0;
    else
      ImagEigenvalues[i] = (*evals)[NumBlocks+i];    
  }
  
  // populate real and imaginary components of the eigenvectors
  // if the uses has passed RealEigenvalues or ImagEigenvalues not null
  if( RealEigenvectors && NumRealEigenvectors ) {

    *NumRealEigenvectors = 0;    

    Anasazi::EpetraMultiVec* evecs = dynamic_cast<Anasazi::EpetraMultiVec *>(MyProblem->GetEvecs().get()); 
    Orthogonalize(*evecs);
    
    int NumRows = EigenVectors.Map().NumMyPoints();
    for( int i=0 ; i<NumBlocks ; ++i ) {
      double norm;
      (*evecs)(i)->Norm2(&norm);
      if( fabs(norm)>1e-8 )  {
	for( int j=0 ; j<NumRows ; ++j ) {
	  RealEigenvectors[(*NumRealEigenvectors)*NumRows+j] = (*evecs)[i][j];
	}
	++(*NumRealEigenvectors);
      }    
    }
  }
  
  if( ImagEigenvectors && NumImagEigenvectors ) {

    *NumImagEigenvectors=0;

    if ( isSymmetric ) {
      Anasazi::EpetraMultiVec* evecs = dynamic_cast<Anasazi::EpetraMultiVec *>(MyProblem->GetEvecs().get()); 
      Orthogonalize(*evecs);
    
      int NumRows = EigenVectors.Map().NumMyPoints();    
      for( int i=0 ; i<NumBlocks ; ++i ) {
        double norm;
        (*evecs)(NumBlocks+i)->Norm2(&norm);
        if( fabs(norm)>1e-8 ) {
	  for( int j=0 ; j<EigenVectors.Map().NumMyPoints() ; ++j ) {
	    ImagEigenvectors[(*NumImagEigenvectors)*NumRows+j] = (*evecs)[i][j];
	  }
	  ++(*NumImagEigenvectors);
        }
      }
    }
  }
  
  // Output results to screen
  if( PrintCurrentStatus && MyPID == 0 ) MyBlockKrylovSchur1.currentStatus();

  /* ********************************************************************** */
  /* Scale vectors (real and imag part) so that max is 1 for each of them   */
  /* ********************************************************************** */

  /* ERASEME
  if( List.get("eigen-analysis: normalize eigenvectors",false) ) {
    
    double * MaxVals = new double[NumBlocks];
    double * MinVals = new double[NumBlocks];
    double * Vals = new double[NumBlocks];
    EigenVectors.MaxValue(MaxVals);
    EigenVectors.MinValue(MinVals);

    for( int i=0 ; i<NumBlocks ; ++i ) {
      Vals[i] = EPETRA_MAX(fabs(MaxVals[i]),fabs(MinVals[i]));
      if( Vals[i] == fabs(MaxVals[i]) && MaxVals[i]<0 ) Vals[i] *= -1.0;
      if( Vals[i] == fabs(MinVals[i]) && MinVals[i]<0 ) Vals[i] *= -1.0;
      Vals[i] = 1.0/Vals[i];
      EigenVectors(i)->Scale(Vals[i]);

    }

    delete [] MaxVals;
    delete [] MinVals;
    delete [] Vals;
  }
  */
  
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

  int output = List_.get("ML output", -47);  
  if (output == -47) output = List_.get("output", 5);  

  
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
  
  Teuchos::RefCountPtr<Anasazi::EpetraMultiVec> Vectors =
    Teuchos::rcp( new Anasazi::EpetraMultiVec(RowMatrix->RowMatrixRowMap(), 1) );
  Vectors->MvRandom();
  
  typedef Anasazi::MultiVec<double> MV;
  typedef Anasazi::Operator<double> OP;

  Teuchos::ParameterList AnasaziPL;
  AnasaziPL.set("Block Size", 1);
  AnasaziPL.set("Max Blocks", length);
  AnasaziPL.set("Max Restarts", restarts);
  AnasaziPL.set("Tol", tol);

  // Call the ctor that calls the petra ctor for a matrix
  Teuchos::RefCountPtr<MLMat> Amat = 
    Teuchos::rcp( new MLMat(*RowMatrix,A_PLUS_AT_MATRIX,UseDiagScaling) );
  
  // Create default output manager 
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM = 
    Teuchos::rcp( new Anasazi::OutputManager<double>() );

  // Create a sorting manager to handle the sorting of eigenvalues in the solver
  Teuchos::RefCountPtr<Anasazi::BasicSort<double, MV, OP> > MySort = 
    Teuchos::rcp( new Anasazi::BasicSort<double, MV, OP>("LM") );
  
  // Create the eigenproblem
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(Amat, Vectors) );
  MyProblem->SetSymmetric(true);
  MyProblem->SetNEV( 1 );
  MyProblem->SetProblem();

  // Initialize the Block Arnoldi solver
  Anasazi::BlockKrylovSchur<double, MV, OP> MyBlockKrylovSchur1(MyProblem, MySort, MyOM, AnasaziPL);

  // Solve the problem to the specified tolerances or length
  MyBlockKrylovSchur1.solve();
  
  // Obtain results directly
  Teuchos::RefCountPtr<std::vector<double> > evalr = MyProblem->GetEvals();

  MaxReal = (*evalr)[0] / 2;

  Teuchos::RefCountPtr<const std::vector<double> > residuals = MyBlockKrylovSchur1.GetRitzResiduals();

  if( output > 5 && MyPID == 0 ) {
    cout << "ML_Anasazi : Ritz Residual for A^T + A = " << (*residuals)[0] << endl;
  }

  if( PrintCurrentStatus && MyPID == 0 ) MyBlockKrylovSchur1.currentStatus();

  /* ********************************************************************** */
  /* First compute A - A^T to get the real bound                            */
  /* ********************************************************************** */
  
  if( isSymmetric == false ) {
    
    if( output > 5 && MyPID == 0 ) {
      cout << "ML_Anasazi : Computing eigenvalues of A - A^T" << endl;
    }
    
    Vectors->MvRandom();
  
    // Call the ctor that calls the petra ctor for a matrix
    Teuchos::RefCountPtr<MLMat> Amat2 = 
      Teuchos::rcp( new MLMat(*RowMatrix,A_MINUS_AT_MATRIX,UseDiagScaling) );
    
    Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem2 =
      Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(Amat2, Vectors) );
    MyProblem2->SetSymmetric(false);
    MyProblem2->SetNEV( 1 );
    MyProblem2->SetProblem();
    
    // Initialize the Block Arnoldi solver
    Anasazi::BlockKrylovSchur<double, MV, OP> MyBlockKrylovSchur2(MyProblem2, MySort, MyOM, AnasaziPL);
    
    // Solve the problem to the specified tolerances or length
    MyBlockKrylovSchur2.solve();
    
    // Obtain results directly
    // The vector is 2*nev in length if problem is non-symmetric
    Teuchos::RefCountPtr<std::vector<double> > evals = MyProblem2->GetEvals();
       
    MaxImag = fabs((*evals)[1] / 2);

    Teuchos::RefCountPtr<const std::vector<double> > residuals = MyBlockKrylovSchur2.GetRitzResiduals();
    if( output > 5 && MyPID == 0 ) {
      cout << "ML_Anasazi : Ritz Residual for A^T - A = " << (*residuals)[0] << endl;
    }
    
    if( PrintCurrentStatus && MyPID == 0 ) MyBlockKrylovSchur2.currentStatus();

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

  
#ifndef ML_CPP
#ifdef __cplusplus
extern "C"
{  
#endif
#endif

int ML_Anasazi_Get_FieldOfValuesBox_Interface(ML_Operator * Amat,
					      struct ML_Field_Of_Values * fov )
{

  Epetra_CrsMatrix * CrsTemp;
  int MaxNumNonzeros;
  double CPUTime;
  
  // FIXME: use ML_Epetra::RowMatrix!
  ML_Operator2EpetraCrsMatrix(Amat,CrsTemp,MaxNumNonzeros,
			      true,CPUTime);
  
  double MaxReal,MaxImag;
  Teuchos::ParameterList * EigenList = (Teuchos::ParameterList *) fov->EigenList;
  
  ML_Anasazi::GetFieldOfValuesBox(CrsTemp,MaxReal,MaxImag,*EigenList);

  if (MaxReal < 0.0) {
    if (CrsTemp->Comm().MyPID() == 0)
      cout << "Warning (GetFieldOfValuesBox) : MaxReal was negative!" << endl;
    MaxReal = -MaxReal;
  }
  if (MaxImag < 0.0) {
    if (CrsTemp->Comm().MyPID() == 0)
      cout << "Warning (GetFieldOfValuesBox) : MaxImag was negative!" << endl;
    MaxImag = -MaxImag;
  }

  double eta = MaxImag/MaxReal;

  fov->eta = eta;
  fov->real_max = MaxReal;
  fov->imag_max = MaxImag;

  delete CrsTemp;
  
  return 0;
 
}

int ML_Anasazi_Get_FieldOfValuesBoxNonScaled_Interface(ML_Operator * Amat,
						       struct ML_Field_Of_Values * fov )
{

  Epetra_CrsMatrix * CrsTemp;
  int MaxNumNonzeros;
  double CPUTime;
  
  // FIXME: use ML_Epetra::RowMatrix!
  ML_Operator2EpetraCrsMatrix(Amat,CrsTemp,MaxNumNonzeros,
			      true,CPUTime);
  
  double MaxReal,MaxImag;
  Teuchos::ParameterList * EigenList = (Teuchos::ParameterList *) fov->EigenList;

  bool UseDiagScaling = EigenList->get("field-of-values: use diagonal scaling", false);
  EigenList->set("field-of-values: use diagonal scaling", false);
  
  ML_Anasazi::GetFieldOfValuesBox(CrsTemp,MaxReal,MaxImag,*EigenList);

  EigenList->set("field-of-values: use diagonal scaling", UseDiagScaling);
  
  if (MaxReal < 0.0) {
    if (CrsTemp->Comm().MyPID() == 0)
      cout << "Warning (GetFieldOfValuesBox) : MaxReal was negative!" << endl;
    MaxReal = -MaxReal;
  }
  if (MaxImag < 0.0) {
    if (CrsTemp->Comm().MyPID() == 0)
      cout << "Warning (GetFieldOfValuesBox) : MaxReal was negative!" << endl;
    MaxImag = -MaxImag;
  }

  double eta = MaxImag/MaxReal;

  fov->eta = eta;
  fov->real_max = MaxReal;
  fov->imag_max = MaxImag;

  delete CrsTemp;
  
  return 0;
 
}

// ================================================ ====== ==== ==== == =

int ML_Anasazi_Get_SpectralNorm_Anasazi(ML_Operator* Amat,
                                        ML_Smoother* smoother,
					int MaxIters, double Tolerance,
                                        int IsProblemSymmetric,
                                        int UseDiagonalScaling,
					double* LambdaMax )
{

  Epetra_CrsMatrix * RowMatrix;
  int MaxNumNonzeros;
  double CPUTime;
  
  // FIXME: replace with ML_Epetra::RowMatrix
  ML_Operator2EpetraCrsMatrix(Amat,RowMatrix,MaxNumNonzeros,
			      true,CPUTime);
  
  bool verbose =  RowMatrix->Comm().MyPID() == 0 && (5 < ML_Get_PrintLevel());
  
  int length = MaxIters;
  int restarts = 1;

  if( verbose ) {
    cout << "ML_Anasazi : Estimate Lambda Max, ";
    if( IsProblemSymmetric == ML_TRUE ) cout << "problem is symmetric, ";
    if( UseDiagonalScaling == ML_TRUE ) cout << "diagonal scaling, ";
    cout << " max iters = " << MaxIters << endl;
  }

  /* ********************************************************************** */
  /* First compute A + A^T to get the real bound                            */
  /* ********************************************************************** */
  
  Teuchos::RefCountPtr<Anasazi::EpetraMultiVec> Vectors =
    Teuchos::rcp( new Anasazi::EpetraMultiVec(RowMatrix->RowMatrixRowMap(), 1) );
  Vectors->MvRandom();
  
  typedef Anasazi::MultiVec<double> MV;
  typedef Anasazi::Operator<double> OP;

  Teuchos::ParameterList AnasaziPL;
  AnasaziPL.set("Block Size", 1);
  AnasaziPL.set("Max Blocks", length);
  AnasaziPL.set("Max Restarts", restarts);
  AnasaziPL.set("Tol", Tolerance);

  // Create default output manager 
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM = 
    Teuchos::rcp( new Anasazi::OutputManager<double>() );

  // Create a sorting manager to handle the sorting of eigenvalues in the solver
  Teuchos::RefCountPtr<Anasazi::BasicSort<double, MV, OP> > MySort = 
    Teuchos::rcp( new Anasazi::BasicSort<double, MV, OP>("LM") );
  
  // Call the ctor that calls the petra ctor for a matrix

  bool flag = false;
  if( UseDiagonalScaling == ML_TRUE ) flag = true;
  
  // smoother is generally 0 (NULL)
  Teuchos::RefCountPtr<MLMat> MLMatrix =
    Teuchos::rcp( new MLMat(*RowMatrix,SMOOTHED_MATRIX,flag,0,smoother) );
  
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(MLMatrix, Vectors) );
  MyProblem->SetNEV( 1 );

  // Inform the problem that it is symmetric
  if( IsProblemSymmetric == ML_TRUE ) MyProblem->SetSymmetric(true);
  else                                MyProblem->SetSymmetric(false);

  MyProblem->SetProblem();
  
  // Initialize the Block Arnoldi solver
  Anasazi::BlockKrylovSchur<double, MV, OP> MyBlockKrylovSchur(MyProblem, MySort, MyOM, AnasaziPL);

  // Solve the problem to the specified tolerances or length
  MyBlockKrylovSchur.solve();
  
  // Obtain results direc3tly

  Teuchos::RefCountPtr<std::vector<double> > evals = MyProblem->GetEvals();

  if ( IsProblemSymmetric == ML_TRUE )
    *LambdaMax = sqrt((*evals)[0]*(*evals)[0]);
  else
    *LambdaMax = sqrt((*evals)[0]*(*evals)[0] + (*evals)[1]*(*evals)[1]);
  
  if (verbose) {
    Teuchos::RefCountPtr<const std::vector<double> > residuals = MyBlockKrylovSchur.GetRitzResiduals();
    cout << "ML_Anasazi : Ritz Residual = " << (*residuals)[0] << endl;
  }
  
  delete RowMatrix;

  return 0;
 
}
#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif
  
#endif
