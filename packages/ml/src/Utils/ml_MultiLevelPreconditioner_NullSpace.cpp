/*
 * Methods to set and/or compute the null space.
 *
 * \author Marzio Sala, SNL 9214
 *
 * \date Last updated on 19-Jan-05
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

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
#include "ml_rbm.h"

using namespace Teuchos;

// ================================================ ====== ==== ==== == =

int ML_Epetra::MultiLevelPreconditioner::SetNullSpace()
{
  int NullSpaceDim = NumPDEEqns_;
  double* NullSpacePtr = 0;

  std::string option = List_.get("null space: type", "default vectors");

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
    ML_Aggregate_Set_NullSpace(agg_,NumPDEEqns_,NullSpaceDim,NULL,
        RowMatrix_->NumMyRows());
    if (verbose_)
      std::cout << PrintMsg_ << "Null space type      = default (constants)" << std::endl;
  } //"default vectors"

  else if (option == "elasticity from coordinates")
  {
    //TODO could add an option to provide center of rotations
    double *in_x_coord = List_.get("x-coordinates", (double *)0);
    double *in_y_coord = List_.get("y-coordinates", (double *)0);
    double *in_z_coord = List_.get("z-coordinates", (double *)0);

    //NullSpaceDim = List_.get("null space: dimension", NumPDEEqns_);
    if      (in_y_coord == 0) NullSpaceDim = 1;
    else if (in_z_coord == 0) NullSpaceDim = 3;
    else                      NullSpaceDim = 6;

    if (in_x_coord == 0) {
      if (Comm().MyPID() == 0)
        std::cerr << ErrorMsg_ << "You asked for the near null-space from coordinates," << std::endl
          << ErrorMsg_ << "but x-coordinate vector is NULL!" << std::endl;
      ML_EXIT(EXIT_FAILURE);
    }

    NullSpacePtr = new double[NullSpaceDim*NumMyRows()];
    ML_Coord2RBM(NumMyRows()/NumPDEEqns_, in_x_coord, in_y_coord, in_z_coord, NullSpacePtr,
        NumPDEEqns_, 0);
    ML_Aggregate_Set_NullSpace(agg_,NumPDEEqns_,NullSpaceDim,NullSpacePtr, RowMatrix_->NumMyRows());
    if (verbose_) {
      std::cout << PrintMsg_ << "Null space type      = elasticity from coordinates" << std::endl;
      std::cout << PrintMsg_ << "  (This option ignores any user-specified nullspace dimension.)" << std::endl;
    }

  } //"elasticity from coordinates"

  else if (option == "pre-computed")
  {
    NullSpaceDim = List_.get("null space: dimension", NumPDEEqns_);
    NullSpacePtr = List_.get("null space: vectors", NullSpacePtr);

    if (NullSpacePtr == 0) {
      if (Comm().MyPID() == 0)
        std::cerr << ErrorMsg_ << "Null space vectors is NULL!" << std::endl;
      ML_EXIT(EXIT_FAILURE);
    }

    ML_Aggregate_Set_NullSpace(agg_,NumPDEEqns_,NullSpaceDim,NullSpacePtr,
        RowMatrix_->NumMyRows());
    if (verbose_)
      std::cout << PrintMsg_ << "Null space type      = user-supplied" << std::endl;

  } //"pre-computed"
#if FIXME
  else if (option == "enriched")
  {
    NullSpaceDim = List_.get("null space: vectors to compute", 1);

    // by default, 0 means to compute one eigenvector per equation.
    // This number can be doubled if the imaginary part is added.
    // Although, eigenvectors with 0-norm are discarded.
    if( NullSpaceDim == 0 ) NullSpaceDim = NumPDEEqns_;

#ifdef HAVE_ML_ANASAxI

    Epetra_Time Time(Comm());

    bool UseDefaultVectors = List_.get("null space: add default vectors", true);

    // NOTE: NullSpaceDim always refers to the number of eigenvectors,
    //       if this flag is true we will keep also the imaginary part
    bool UseImaginaryComponents = List_.get("null space: add imaginary components", true);

    if( verbose_ ) {
      std::cout << PrintMsg_ << "Enriching null space with " << NullSpaceDim << " vector(s)";
      if( UseImaginaryComponents ) std::cout << PrintMsg_ << ", both real and imaginary components";
      if( UseDefaultVectors ) std::cout << PrintMsg_ << std::endl << "plus " << NumPDEEqns_ << " constant vector(s)" << std::endl;
      else std::cout << std::endl;
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
      std::cerr << ErrorMsg_ << "NullSpacePtr is not NULL. Is null space already defined?" << std::endl
        << ErrorMsg_ << "Now I delete the old null space, and proceed with finger crossed..." << std::endl;
      delete [] NullSpacePtr;
    }

    NullSpacePtr = new double[TotalNullSpaceDim*NumMyRows()];

    if( NullSpacePtr == 0 ) {
      std::cerr << ErrorMsg_ << "Not enough space to allocate " << TotalNullSpaceDim*NumMyRows()*8
        << " bytes" << std::endl
        << "(file " << __FILE__ << ", line " << __LINE__ << ")" << std::endl;
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

      std::string opt = List_.get("null space: matrix operation", "I-A");
      if( opt == "I-A" ) {
        AnasaziList.set("eigen-analysis: matrix operation", opt);
        AnasaziList.set("eigen-analysis: action", "LM");
        AnasaziList.set("eigen-analysis: use diagonal scaling",false);
      } else if( opt == "I-D^{-1}A" ) {
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
        std::cerr << ErrorMsg_ << "value for `null space: matrix operation' not recognized" << std::endl
          << ErrorMsg_ << "(" << opt << "). It should be: " << std::endl
          << ErrorMsg_ << "<I-A> / <A> / <D^{-1}A> / <I-D^{-1}A>" << std::endl;
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
      std::cout << PrintMsg_ << "\tComputed eigenvalues:" << std::endl;
      for( int i=0 ; i<NullSpaceDim ; ++i ) {
        std::cout << PrintMsg_ << "\t" << std::setw(10) << RealEigenvalues[i]
          << " + " << std::setw(10) << ImagEigenvalues[i] << " i" << std::endl;
        if( RealEigenvalues[i] == 0 ) ++Discarded;
        if( ImagEigenvalues[i] == 0 ) ++Discarded;
      }
      std::cout << std::endl;
    }

    if( Discarded && verbose_ ) {
      std::cout << PrintMsg_ << "Discarded " << Discarded << " eigenvectors" << std::endl;
    }


    if( 0 ) { // debugging only, print all computed eigenvectors
      for( int i=0 ; i<NumMyRows() ; ++i ) {
        std::cout << i << ": ";
        for( int j=0 ; j<TotalNullSpaceDim-Discarded ; ++j ) {
          std::cout << NullSpacePtr[i+j*NumMyRows()] << " ";
        }
        std::cout << std::endl;
      }
    }


    ML_Aggregate_Set_NullSpace(agg_,NumPDEEqns_,TotalNullSpaceDim-Discarded,
        NullSpacePtr,
        NumMyRows());

    delete [] RealEigenvalues;
    delete [] ImagEigenvalues;

    if( verbose_ ) std::cout << PrintMsg_ << "Total time for eigen-analysis = " << Time.ElapsedTime() << " (s)\n";

#else
    std::cout << "ML_Anasazi ERROR: you must compile with --with-ml_anasazi "  << std::endl
      << "ML_Anasazi ERROR: for eigen-analysis." << std::endl;
    exit( EXIT_FAILURE );
#endif

  }
#endif
  else
  {
    std::cerr << ErrorMsg_ << "Option `null space: type' not recognized ("
      << option << ")" << std::endl
      << ErrorMsg_ << "It should be:" << std::endl
      << ErrorMsg_ << "<default vectors> / <pre-computed> / <elasticity from coordinates> / <enriched>" << std::endl;
    exit(EXIT_FAILURE);
  }

  // May need to scale the null space ??

  double * NullSpaceScaling = List_.get("null space: scaling", (double *)0);

  if (NullSpaceScaling != 0)
  {
    if (verbose_)
      std::cout << PrintMsg_ << "Scaling Null Space..." << std::endl;
    ML_Aggregate_Scale_NullSpace(agg_,NullSpaceScaling,
        RowMatrix_->RowMatrixRowMap().NumMyElements());
  }

  if (verbose_) std::cout << PrintMsg_ << "Null space dimension = " << NullSpaceDim << std::endl;

  return(0);
}

#endif /*ifdef HAVE_ML_EPETRA && HAVE_ML_TEUCHOS*/
