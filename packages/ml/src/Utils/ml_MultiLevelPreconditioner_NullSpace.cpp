/*
 * Methods to set and/or compute the null space.
 * 
 * \author Marzio Sala, SNL 9214
 * 
 * \date Last updated on 19-Jan-05
 */

#include "ml_common.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "ml_epetra.h"
#include "ml_epetra_utils.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_anasazi.h"

using namespace Teuchos;

// ================================================ ====== ==== ==== == =

int ML_Epetra::MultiLevelPreconditioner::SetNullSpace() 
{
  int NullSpaceDim = NumPDEEqns_;
  double* NullSpacePtr = 0;

  string option = List_.get("null space: type", "default vectors");

  // to save time, the 1-level case will always use "default vectors"
  if (NumLevels_ == 1) option = "default vectors";
  
  // Null space can be obtained in 3 ways:
  // 1. default vectors, one constant vector for each physical unknown
  // 2. precomputed, the user is furnishing a pointer to a double vector,
  //    containing all the required components
  // 3. by computing the eigenvalues of a suitable matrix (for instance, the
  //    lowest of A, or the largest of I-A). Default space can be added
  //    if required.
  
  if (option == "default vectors") 
  {
    // sanity check for default null-space
    if( NullSpacePtr == NULL ) NullSpaceDim = NumPDEEqns_;
    ML_Aggregate_Set_NullSpace(agg_,NumPDEEqns_,NumPDEEqns_,NULL,
                               RowMatrix_->NumMyRows());

  } 
  else if (option == "pre-computed") 
  {
    NullSpaceDim = List_.get("null space: dimension", NumPDEEqns_);
    NullSpacePtr = List_.get("null space: vectors", NullSpacePtr);

    if (verbose_) {
      cout << PrintMsg_ << "Using pre-computed null space of dimension "
           << NullSpaceDim << endl;
    }
    if (NullSpacePtr == 0) {
      if (Comm().MyPID() == 0) 
        cerr << ErrorMsg_ << "Null space vectors is NULL!" << endl;
      ML_EXIT(EXIT_FAILURE);
    }
    
    ML_Aggregate_Set_NullSpace(agg_,NumPDEEqns_,NullSpaceDim,NullSpacePtr,
			       RowMatrix_->NumMyRows());
  
  } 
#if FIXME
  else if (option == "enriched") 
  {
    NullSpaceDim = List_.get("null space: vectors to compute", 1);

    // by default, 0 means to compute one eigenvector per equation.
    // This number can be doubled if the imaginary part is added.
    // Although, eigenvectors with 0-norm are discarded.
    if( NullSpaceDim == 0 ) NullSpaceDim = NumPDEEqns_;
    
#ifdef HAVE_ML_ANASAZI

    Epetra_Time Time(Comm());
    
    bool UseDefaultVectors = List_.get("null space: add default vectors", true);

    // NOTE: NullSpaceDim always refers to the number of eigenvectors,
    //       if this flag is true we will keep also the imaginary part
    bool UseImaginaryComponents = List_.get("null space: add imaginary components", true);
    
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
      
      string opt = List_.get("null space: matrix operation", "I-A");
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

  } 
#endif
  else 
  {
    cerr << ErrorMsg_ << "Option `null space: type' not recognized ("
	 << option << ")" << endl
	 << ErrorMsg_ << "It should be:" << endl
	 << ErrorMsg_ << "<default vectors> / <pre-computed> / <enriched>" << endl;
    exit(EXIT_FAILURE);
  }
  
  // May need to scale the null space ??

  double * NullSpaceScaling = List_.get("null space: scaling", (double *)0);

  if (NullSpaceScaling != 0) 
  {
    if (verbose_) 
      cout << PrintMsg_ << "Scaling Null Space..." << endl;
    ML_Aggregate_Scale_NullSpace(agg_,NullSpaceScaling,
                                 RowMatrix_->RowMatrixRowMap().NumMyElements());
  } 

  return(0);
}

#endif /*ifdef HAVE_ML_EPETRA && HAVE_ML_TEUCHOS*/
