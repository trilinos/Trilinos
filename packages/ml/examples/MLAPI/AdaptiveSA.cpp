
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

#include "EpetraExt_CrsMatrixIn.h"

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

  try {

    // Initialize the workspace and set the output level
    Init();

    RCP<Operator> A;
    if(argc == 2){
      Epetra_CrsMatrix * temp;
      EpetraExt::MatrixMarketFileToCrsMatrix(argv[1],GetEpetra_Comm(),temp);
      temp->FillComplete();
      Space DomainSpace(temp->DomainMap());
      Space RangeSpace(temp->RangeMap());
      A = rcp(new Operator(DomainSpace,RangeSpace,temp));
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
    

    int NumPDEEqns = 2;
    int MaxLevels = 10;

    Teuchos::ParameterList List;
    List.set("additional candidates", 2);
    List.set("use default null space", true);
    List.set("krylov: type", "cg");
    RCP<MultiLevelAdaptiveSA> Prec;
    Prec=rcp(new MultiLevelAdaptiveSA(*A, List, NumPDEEqns, MaxLevels));

    // =============================================================== //
    // setup the hierarchy:                                            //
    // - `UseDefaultNullSpace' toggles the use of default candidates.  //
    // - AdditionalCandidates = 2' means to compute two additionals.   //
    // - the final null space dimension is 3.                          //
    // =============================================================== //

    bool UseDefaultNullSpace = true;
    int AdditionalCandidates = 1;
    Prec->AdaptCompute(UseDefaultNullSpace, AdditionalCandidates);

    MultiVector LHS(A->GetDomainSpace());
    MultiVector RHS(A->GetRangeSpace());

    LHS.Random();
    RHS = 0.0;

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
