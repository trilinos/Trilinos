/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/************************************************************************/
/*          Utilities for Trilinos/ML users                             */
/*----------------------------------------------------------------------*/
/* Author :  Marzio Sala (SNL)                                          */
/************************************************************************/

#include "ml_common.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

//#include <cstring>
#include "ml_amesos_wrap.h"
#include "ml_ifpack_wrap.h"
#include "ml_agg_METIS.h"
#include "ml_epetra_utils.h"

#include "ml_epetra.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_agg_ParMETIS.h"

#include "ml_anasazi.h"

using namespace Teuchos;

// ================================================ ====== ==== ==== == =

int ML_Epetra::MultiLevelPreconditioner::SetNullSpaceMaxwell()
{

  // FIXME...
  return(0);
}

// ================================================ ====== ==== ==== == =

int ML_Epetra::MultiLevelPreconditioner::SetNullSpace() 
{

  char parameter[80];
  
  const Epetra_VbrMatrix * VbrMatrix = dynamic_cast<const Epetra_VbrMatrix *>(RowMatrix_);
  if( VbrMatrix == 0 ) {
    sprintf(parameter,"%sPDE equations", Prefix_);
    NumPDEEqns_ = List_.get(parameter, 1);
  }
  else {
    int NumBlockRows = VbrMatrix->RowMap().NumGlobalElements();
    int NumRows = VbrMatrix->RowMap().NumGlobalPoints();
    if( NumRows % NumBlockRows ) {
      cerr << "# rows must be a multiple of # block rows ("
	   << NumRows << "," << NumBlockRows << ")" << endl;
      exit( EXIT_FAILURE );
    }
    NumPDEEqns_ = NumRows/NumBlockRows;
  }

  if( verbose_ ) cout << "Number of PDE equations = " << NumPDEEqns_ << endl;
  
  int NullSpaceDim = NumPDEEqns_;
  double * NullSpacePtr = NULL;

  sprintf(parameter,"%snull space: type", Prefix_);
  string option = List_.get(parameter, "default vectors");

  // to save time, the 1-level case will always use "default vectors"
  if( NumLevels_ == 1 ) option = "default vectors";
  
  // Null space can be obtained in 3 ways:
  // 1. default vectors, one constant vector for each physical unknown
  // 2. precomputed, the user is furnishing a pointer to a double vector,
  //    containing all the required components
  // 3. by computing the eigenvalues of a suitable matrix (for instance, the
  //    lowest of A, or the largest of I-A). Default space can be added
  //    if required.
  
  // FIXME: null space scaling now works only with default vectors...

  if (option != "default vectors" && Scaling_) {
    cerr << ErrorMsg_ << "Scaling the null spaceworks only" << endl
         << ErrorMsg_ << "for `default vectors' (at present)..." << endl;
  }

  if (option == "default vectors") {

    if (Scaling_) {

      if (verbose_)
	cout << PrintMsg_ << "Scaling default null space..." << endl;

      NullSpaceToFree_ = new double[NumPDEEqns_*NumMyRows()];
      // fill it with normal 0's and 1's for standard vectors
      for( int i=0 ; i<NumPDEEqns_ ; ++i )
	for( int j=0 ; j<NumMyRows() ; ++j )
	  if (j%NumPDEEqns_ == i) 
	    NullSpaceToFree_[j + i * NumMyRows()] = 1.0 * (*InvScaling_)[j];
	  else                     
	    NullSpaceToFree_[j + i * NumMyRows()] = 0.0;

      ML_Aggregate_Set_NullSpace(agg_,NumPDEEqns_,NumPDEEqns_,
				 NullSpaceToFree_,
				 RowMatrix_->NumMyRows());

    }
    else {
      // sanity check for default null-space
      if( NullSpacePtr == NULL ) NullSpaceDim = NumPDEEqns_;
      ML_Aggregate_Set_NullSpace(agg_,NumPDEEqns_,NumPDEEqns_,NULL,
				 RowMatrix_->NumMyRows());
    }
    
  } else if( option == "pre-computed" ) {

    sprintf(parameter,"%snull space: dimension", Prefix_);    
    NullSpaceDim = List_.get(parameter, NumPDEEqns_);
    sprintf(parameter,"%snull space: vectors", Prefix_);
    NullSpacePtr = List_.get(parameter, NullSpacePtr);

    if( NullSpacePtr == 0 ) {
      if( Comm().MyPID() == 0 ) cerr << ErrorMsg_ << "Null space vectors is NULL!" << endl;
      exit( EXIT_FAILURE );
    }
    
    ML_Aggregate_Set_NullSpace(agg_,NumPDEEqns_,NullSpaceDim,NullSpacePtr,
			       RowMatrix_->NumMyRows());
  
  } else if( option == "enriched" ) {

    sprintf(parameter,"%snull space: vectors to compute", Prefix_);    
    NullSpaceDim = List_.get(parameter, 1);

    // by default, 0 means to compute one eigenvector per equation.
    // This number can be doubled if the imaginary part is added.
    // Although, eigenvectors with 0-norm are discarded.
    if( NullSpaceDim == 0 ) NullSpaceDim = NumPDEEqns_;
    
#ifdef HAVE_ML_ANASAZI

    Epetra_Time Time(Comm());
    
    sprintf(parameter,"%snull space: add default vectors", Prefix_);
    bool UseDefaultVectors = List_.get(parameter, true);

    // NOTE: NullSpaceDim always refers to the number of eigenvectors,
    //       if this flag is true we will keep also the imaginary part
    sprintf(parameter,"%snull space: add imaginary components", Prefix_);
    bool UseImaginaryComponents = List_.get(parameter, true);
    
    if( verbose_ ) {
      cout << PrintMsg_ << "Enriching null space with " << NullSpaceDim << " vector(s)";
      if( UseImaginaryComponents ) cout << PrintMsg_ << ", both real and imaginary components";
      if( UseDefaultVectors ) cout << PrintMsg_ << endl << "plus " << NumPDEEqns_ << " constant vector(s)" << endl;
      else cout << endl;
    }

    // allocate space for the entire null space, that contains:
    // 1- the default one, a constant for each unknown (if required)
    // 2- the real part of the Anasazi computations
    // 3- the imaginary part of the Anasazi computions (if required)
    int TotalNullSpaceDim = NullSpaceDim;
    if( UseImaginaryComponents ) TotalNullSpaceDim *= 2;
    if( UseDefaultVectors ) TotalNullSpaceDim += NumPDEEqns_;

    // create a double vector hosting null space
    if( NullSpacePtr ) {
      cerr << ErrorMsg_ << "NullSpacePtr is not NULL. Is null space already defined?" << endl
	   << ErrorMsg_ << "Now I delete the old null space, and proceed with finger crossed..." << endl;
      delete [] NullSpacePtr;
    }

    NullSpacePtr = new double[TotalNullSpaceDim*NumMyRows()];
    
    if( NullSpacePtr == 0 ) {
      cerr << ErrorMsg_ << "Not enough space to allocate " << TotalNullSpaceDim*NumMyRows()*8
	   << " bytes" << endl
	   << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
    }

    // here NullSpaceDim is the number of eigenvalues/vectors that Anasazi has to compute

    // will contain the eigenVALUES
    double * RealEigenvalues = new double[NullSpaceDim];
    double * ImagEigenvalues = new double[NullSpaceDim];

    /*
    int offset;
    if( UseDefaultVectors ) offset = NumPDEEqns_;
    else                    offset = 0;
    */

    // create List for Anasazi (kept separate from List_, I don't want to pollute it)
    ParameterList AnasaziList;
    
    {
      
      sprintf(parameter,"%snull space: matrix operation", Prefix_);
      string opt = List_.get(parameter, "I-A");
      if( opt == "I-A" ) {
	AnasaziList.set("eigen-analysis: matrix operation", opt);
	AnasaziList.set("eigen-analysis: action", "LM");
	AnasaziList.set("eigen-analysis: use diagonal scaling",false);
      }	else if( opt == "I-D^{-1}A" ) {
	AnasaziList.set("eigen-analysis: matrix operation", "I-A");
	AnasaziList.set("eigen-analysis: action", "LM");
	AnasaziList.set("eigen-analysis: use diagonal scaling",true);
      } else if( opt == "A" ) {
	AnasaziList.set("eigen-analysis: matrix operation", opt);
	AnasaziList.set("eigen-analysis: action", "SM");	
	AnasaziList.set("eigen-analysis: use diagonal scaling",false);
      } else if( opt == "D^{-1}A" ) {
	AnasaziList.set("eigen-analysis: matrix operation", "A");
	AnasaziList.set("eigen-analysis: action", "SM");	
	AnasaziList.set("eigen-analysis: use diagonal scaling",true);
      } else {
	cerr << ErrorMsg_ << "value for `null space: matrix operation' not recognized" << endl
	     << ErrorMsg_ << "(" << opt << "). It should be: " << endl
	     << ErrorMsg_ << "<I-A> / <A> / <D^{-1}A> / <I-D^{-1}A>" << endl;
	exit( EXIT_FAILURE );
      }

      AnasaziList.set("eigen-analysis: length", List_.get("eigen-analysis: length", 20));
      AnasaziList.set("eigen-analysis: tolerance", List_.get("eigen-analysis: tolerance", 1.0e-1));
      AnasaziList.set("eigen-analysis: restart", List_.get("eigen-analysis: restart", 2));
    }
    
    // this is the starting value -- random
    Epetra_MultiVector EigenVectors(OperatorDomainMap(),NullSpaceDim);
    EigenVectors.Random();
    
    // call Anasazi. Real and imaginary part of the selected eigenvalues
    // will be copied into RealEigenvalues

    double * RealEigenvectors = 0, * ImagEigenvectors = 0;
    if( UseDefaultVectors ) RealEigenvectors = NullSpacePtr+NumPDEEqns_*NumMyRows();
    else                    RealEigenvectors = NullSpacePtr;
    if( UseImaginaryComponents ) ImagEigenvectors = RealEigenvectors+NullSpaceDim*NumMyRows();
    else                         ImagEigenvectors = 0;
    
    ML_Anasazi::Interface(RowMatrix_,EigenVectors,RealEigenvalues,
			  ImagEigenvalues, AnasaziList, RealEigenvectors,
			  ImagEigenvectors);

    // fill it with normal 0's and 1's for standard vectors
    if( UseDefaultVectors )
      for( int i=0 ; i<NumPDEEqns_ ; ++i )
	for( int j=0 ; j<NumMyRows() ; ++j )
	  if( j%NumPDEEqns_ == i ) NullSpacePtr[j+i*NumMyRows()] = 1.0;
	  else                     NullSpacePtr[j+i*NumMyRows()] = 0.0;

    NullSpaceToFree_ = NullSpacePtr; // this null space will be freed later

    int Discarded = 0;
    if( verbose_ ) {
      cout << PrintMsg_ << "\tComputed eigenvalues:" << endl;
      for( int i=0 ; i<NullSpaceDim ; ++i ) {
	cout << PrintMsg_ << "\t" << std::setw(10) << RealEigenvalues[i]
	     << " + " << std::setw(10) << ImagEigenvalues[i] << " i" << endl;
	if( RealEigenvalues[i] == 0 ) ++Discarded;
	if( ImagEigenvalues[i] == 0 ) ++Discarded;
      }
      cout << endl;
    }

    if( Discarded && verbose_ ) {
      cout << PrintMsg_ << "Discarded " << Discarded << " eigenvectors" << endl;
    }
    
      
    if( 0 ) { // debugging only, print all computed eigenvectors
      for( int i=0 ; i<NumMyRows() ; ++i ) {
	cout << i << ": ";
	for( int j=0 ; j<TotalNullSpaceDim-Discarded ; ++j ) {
	  cout << NullSpacePtr[i+j*NumMyRows()] << " ";
	}
	cout << endl;
      }
    }

    
    ML_Aggregate_Set_NullSpace(agg_,NumPDEEqns_,TotalNullSpaceDim-Discarded,
			       NullSpacePtr,
			       NumMyRows());

    delete [] RealEigenvalues;
    delete [] ImagEigenvalues;
    
    if( verbose_ ) cout << PrintMsg_ << "Total time for eigen-analysis = " << Time.ElapsedTime() << " (s)\n";
    
#else
     cout << "ML_Anasazi ERROR: you must compile with --with-ml_anasazi "  << endl
       << "ML_Anasazi ERROR: for eigen-analysis." << endl;
     exit( EXIT_FAILURE );
#endif

  } else {

    cerr << ErrorMsg_ << "Option `null space: type' not recognized ("
	 << option << ")" << endl
	 << ErrorMsg_ << "It should be:" << endl
	 << ErrorMsg_ << "<default vectors> / <pre-computed> / <enriched>" << endl;
    exit( EXIT_FAILURE );
  }
  
  // May need to scale the null space ??

  sprintf(parameter,"%snull space: scaling", Prefix_);
  double * NullSpaceScaling = List_.get(parameter, (double *)0);

  if( NullSpaceScaling != 0 ) {
    if( verbose_ ) cout << PrintMsg_ << "Scaling Null Space..." << endl;
    ML_Aggregate_Scale_NullSpace(agg_,NullSpaceScaling,RowMatrix_->RowMatrixRowMap().NumMyElements());
  } 

  return(0);
}

// ================================================ ====== ==== ==== == =

int ML_Epetra::MultiLevelPreconditioner::SetEigenList() 
{

  char parameter[80];
  
  // eigen-analysis:
  sprintf(parameter,"%seigen-analysis: use symmetric algorithm", Prefix_);
  bool IsSymmetric = List_.get(parameter,false);
    
  if( IsSymmetric ) EigenList_.set("eigen-analysis: symmetric problem",true);
  else              EigenList_.set("eigen-analysis: symmetric problem",false);

  sprintf(parameter,"%seigen-analysis: tolerance", Prefix_);
  EigenList_.set("eigen-analysis: tolerance", List_.get(parameter, 1e-2));

  sprintf(parameter,"%seigen-analysis: use diagonal scaling", Prefix_);    
  EigenList_.set("eigen-analysis: use diagonal scaling", List_.get(parameter,true));
    
  sprintf(parameter,"%seigen-analysis: restart", Prefix_);
  int itemp = List_.get(parameter, 100);
  EigenList_.set("eigen-analysis: restart", itemp);

  sprintf(parameter,"%seigen-analysis: length", Prefix_);
  itemp =  List_.get(parameter, 20);
  EigenList_.set("eigen-analysis: length", itemp);

  // field of values:

  sprintf(parameter,"%sfield-of-values: tolerance", Prefix_);
  EigenList_.set("field-of-values: tolerance", List_.get(parameter, 1e-2));

  sprintf(parameter,"%sfield-of-values: use diagonal scaling", Prefix_);    
  EigenList_.set("field-of-values: use diagonal scaling", List_.get(parameter,true));
    
  sprintf(parameter,"%sfield-of-values: restart", Prefix_);
  itemp = List_.get(parameter, 100);
  EigenList_.set("field-of-values: restart", itemp);

  sprintf(parameter,"%sfield-of-values: length", Prefix_);
  itemp =  List_.get(parameter, 20);
  EigenList_.set("field-of-values: ", itemp);

  sprintf(parameter,"%sfield-of-values: print current status", Prefix_);
  bool btemp =  List_.get(parameter, false);
  EigenList_.set("field-of-values: print current status", btemp);

  // general output
  
  sprintf(parameter,"%soutput", Prefix_);
  itemp =  List_.get(parameter, 10);
  EigenList_.set("output",itemp);
    
  return(0);
}

int ML_Operator_GetDiagonal(ML_Operator* Amat,
			    double* diagonal)
{

  int allocated = 100;
  int i;
  int j;
  int ierr;
  int NumNonzeros;
  int* colInd;
  double* colVal;

  colInd = (int*)   ML_allocate(sizeof(int)*allocated);
  colVal = (double*)ML_allocate(sizeof(double)*allocated);

  for (i = 0 ; i < Amat->invec_leng ; ++i) {

    ierr = ML_Operator_Getrow(Amat,1,&i,allocated,colInd,colVal,&NumNonzeros);

    if (ierr == 0) {
      do {
	ML_free(colInd);
	ML_free(colVal);
	allocated *= 2;
	colInd = (int*)   ML_allocate(sizeof(int)*allocated);
	colVal = (double*)ML_allocate(sizeof(double)*allocated);

	ierr = ML_Operator_Getrow(Amat,1,&i,allocated,colInd,colVal,
				  &NumNonzeros);
      } while (ierr == 0);
    }
    for (j = 0 ; j < NumNonzeros ; ++j)
      if (colInd[j] == i)
	diagonal[i] = colVal[j];

  }
  
  ML_free(colInd);
  ML_free(colVal);

  return(0);
}

// ================================================ ====== ==== ==== == =

int ML_Epetra::MultiLevelPreconditioner::SetScaling() 
{
  
  int ierr;
  string ScalingType;
 
  ScalingType = List_.get("scaling: type", "none");

  if (ScalingType == "none") 
    return(0);

  Scaling_ = new Epetra_Vector(RowMatrix_->RowMatrixRowMap());
  InvScaling_ = new Epetra_Vector(RowMatrix_->RowMatrixRowMap());
  assert(InvScaling_ != 0);

  if (ScalingType == "col sum") {

    ierr = RowMatrix_->InvColSums(*Scaling_);
    if (ierr) {
      cerr << endl;
      cerr << ErrorMsg_ << "Method InvColSums() returns "
	<< ierr << "." << endl
	<< ErrorMsg_ << "Is this method implemented in your"
	<< " Epetra_RowMatrix-derived class?" << endl;
      cerr << ErrorMsg_ << "Sorry, I must skip the scaling..." << endl;
      cerr << endl;
      ML_CHK_ERR(-1);
    }

    ierr = InvScaling_->Reciprocal(*Scaling_);
    assert (ierr == 0);

  }
  else if (ScalingType == "diagonal" ) {

    InvScaling_->PutScalar(1.0);

    ierr = ML_Operator_GetDiagonal(&(ml_->Amat[LevelID_[0]]),
				   InvScaling_->Values());

    // FIXME: only for non-VBR matrices??
    /* ierr = RowMatrix_->ExtractDiagonalCopy(*InvScaling_); */

    if (ierr) {
      cerr << endl;
      cerr << ErrorMsg_ << "Method ExtractDiagonalCopy() returns "
	<< ierr << "." << endl
	<< ErrorMsg_ << "Is this method implemented in your"
	<< " Epetra_RowMatrix-derived class?" << endl;
      cerr << ErrorMsg_ << "Sorry, I must skip the scaling..." << endl;
      cerr << endl;
      return(-2);
    }

    ierr = Scaling_->Reciprocal(*InvScaling_);
    assert (ierr == 0);

  }
  else {
    cerr << ErrorMsg_ << "Parameter `scaling type' as an incorrect" << endl
         << ErrorMsg_ << "value (" << ScalingType << "). It can be:" << endl
	 << ErrorMsg_ << "<none> / <col sum>" << endl;
    return(-2);
  }

  if (verbose_)
    cout << "Scaling type = " << ScalingType << endl;

  // I scale the matrix right now, it will be scaled back
  // as I return from ComputePreconditioner().
  // If Scaling_ is different from 0, in SetNullSpace() I
  // will take care of creating (or scaling) the null space
  Epetra_RowMatrix* RM = const_cast<Epetra_RowMatrix*>(RowMatrix_);
  ierr = RM->RightScale(*Scaling_);
  assert(ierr == 0);

  return(0);
}
 
#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
