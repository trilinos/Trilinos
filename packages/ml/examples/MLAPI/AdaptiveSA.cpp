
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_config.h"
#include "ml_common.h"
#ifdef HAVE_ML_MLAPI
#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_MultiLevelAdaptiveSA.h"
#include "MLAPI_DistributedMatrix.h"
#include "MLAPI_Krylov.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"
#include "Epetra_Vector.h"

using namespace Teuchos;
using namespace MLAPI;

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif
  Teuchos::CommandLineProcessor clp;
  std::string matrixFile;        clp.setOption("matrix", &matrixFile, "Matrix to solve");
  std::string rhsFile;           clp.setOption("rhs",    &rhsFile,    "RHS to solve");
  int AdditionalCandidates=1;    clp.setOption("numvecs",&AdditionalCandidates,    "Number of additional adaptive candidate vectors");
  int NumPDEs=1;                 clp.setOption("numpdes",&NumPDEs,    "Number of PDE equations");
  clp.recogniseAllOptions(true);
  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }

  try {

    // Initialize the workspace and set the output level
    Init();

    RCP<Operator> A;
    RCP<Epetra_Vector> rhsE;
    if(!matrixFile.empty()) {
      Epetra_CrsMatrix * temp;
      EpetraExt::MatrixMarketFileToCrsMatrix(matrixFile.c_str(),GetEpetra_Comm(),temp);
      temp->FillComplete();
      Space DomainSpace(temp->DomainMap());
      Space RangeSpace(temp->RangeMap());
      A = rcp(new Operator(DomainSpace,RangeSpace,temp));
      if(!rhsFile.empty()) {
        Epetra_Vector * tempV;
        EpetraExt::MatrixMarketFileToVector(rhsFile.c_str(),temp->RowMap(),tempV);
        rhsE = rcp(tempV);
      }
    }
    else { 
      int NX = 1000;
      // define the space for fine level vectors and operators.
      Space FineSpace(2*NX);
      
      DistributedMatrix * Atemp = new DistributedMatrix(FineSpace, FineSpace);
      
      // assemble the matrix on processor 0 only
      if (GetMyPID() == 0)
        {
          for (int i = 0 ; i < NX ; ++i)
            {
              (*Atemp)(2*i, 2*i)     = 2.0;
              (*Atemp)(2*i+1, 2*i+1) = 2.0;
              if (i)
                {
                  (*Atemp)(2*i, 2*(i - 1))     = - 1.0;
                  (*Atemp)(2*i+1, 2*(i - 1)+1) = - 1.0;
                }
              if (i != NX - 1)
                {
                  (*Atemp)(2*i, 2*(i + 1))     = - 1.0;
                  (*Atemp)(2*i+1, 2*(i + 1)+1) = - 1.0;
                }
            }
        }
      Atemp->FillComplete();
      A = rcp(Atemp);
    }// end else
    

    int MaxLevels = 10;

    Teuchos::ParameterList List;
    List.set("additional candidates", AdditionalCandidates);
    List.set("PDE equations",NumPDEs);
    List.set("use default null space", true);
    List.set("krylov: type", "cg");
    SetPrintLevel(10);
    RCP<MultiLevelAdaptiveSA> Prec;
    Prec=rcp(new MultiLevelAdaptiveSA(*A, List, NumPDEs, MaxLevels));

    // =============================================================== //
    // setup the hierarchy:                                            //
    // - `UseDefaultNullSpace' toggles the use of default candidates.  //
    // - AdditionalCandidates = 2' means to compute two additionals.   //
    // - the final null space dimension is 3.                          //
    // =============================================================== //

    bool UseDefaultNullSpace = true;
    Prec->AdaptCompute(UseDefaultNullSpace, AdditionalCandidates);

    MultiVector LHS(A->GetDomainSpace());
    MultiVector RHS(A->GetRangeSpace());


    if(!rhsE.is_null()) {
      for(int i=0; i<rhsE->MyLength(); i++)
        RHS(i,0) = (*rhsE)[i];
      LHS = 0.0;
    }
    else {
      LHS.Random();
      RHS = 0.0;
    }

    Krylov(*A, LHS, RHS, *Prec, List);

    Finalize();

  }
  catch (const int e) {
    std::cerr << "Caught integer exception, code = " << e << std::endl;
  }
  catch (...) {
    std::cerr << "Caught exception..." << std::endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(0);

}

#else

#include "ml_include.h"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("This MLAPI example requires the following configuration options:");
  puts("\t--enable-epetra");
  puts("\t--enable-teuchos");
  puts("\t--enable-ifpack");
  puts("\t--enable-amesos");
  puts("Please check your configure line.");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);
}
#endif // #if defined(HAVE_ML_MLAPI)
