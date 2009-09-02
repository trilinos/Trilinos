#ifndef TRIKOTA_DRIVER
#define TRIKOTA_DRIVER

// Dakota includes

#include "system_defs.h"
#include "ParallelLibrary.H"
#include "ProblemDescDB.H"
#include "DakotaStrategy.H"
#include "DakotaModel.H"
#include "DakotaInterface.H"
#include "DirectApplicInterface.H"


//Trilinos includes
#include "Teuchos_RCP.hpp"

using namespace Dakota;

namespace TriKota {

  //! Class which wraps library-mode Dakota calls into a few simple steps
class Driver {
public:

  //! Constructor, with all dakota filenames having default values

  Driver(const char* dakota_in="dakota.in",  
    const char* dakota_out="dakota.out",
    const char* dakota_err="dakota.err",
    const char* dakota_restart_out="dakota_restart.out");

  ~Driver() {};

  /*! \brief Accessor to get an MPI_Comm from Dakota. This allows Dakota to
     choose the parallelism, and the application to be constructed as a 
     second step using this communicator. If the application is built 
     on MPI_COMM_WORLD, then this call can be used to verify that Dakota
     is running in that mode as well.*/
  MPI_Comm getAnalysisComm(); 

  /*! \brief Accessor to get problem description from Dakota. 
    This hook is used within TriKota::DirectApplicInterface to
    (re)set the initial parameters in Dakota using those selected
    in the model evaluator.
  */
  ProblemDescDB& getProblemDescDB();

  /*! \brief Main call to execute the dakota analysis. 
     The argument may be of type TriKota::DirectApplicInterface, 
     whoch wraps an EpetraExt::ModelEvaluator.
  */
  void run(Dakota::DirectApplicInterface* appInterface);

  //! Accessor for final parameters after an optimization run.
  const Dakota::Variables getFinalSolution() const;

private:

  ParallelLibrary parallel_lib;
  ProblemDescDB problem_db;
  Strategy selected_strategy;

}; // end of class Driver

} // namespace TriKota

#endif //TRIKOTA_DRIVER

