#ifndef TRIKOTA_DIRECTAPPLICINTERFACE
#define TRIKOTA_DIRECTAPPLICINTERFACE

#include "DirectApplicInterface.H"
#include "CommandLineHandler.H"
#include "DakotaStrategy.H"
#include "DakotaModel.H"
#include "ParallelLibrary.H"
#include "ProblemDescDB.H"

#include "EpetraExt_ModelEvaluator.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

/** \brief Interface class inherited from Dakota base class*/
namespace TriKota {

class DirectApplicInterface : public Dakota::DirectApplicInterface
{
public:

  // Constructor and destructor

   DirectApplicInterface(Dakota::ProblemDescDB& problem_db_,
                         const Teuchos::RCP<EpetraExt::ModelEvaluator> App_);

  ~DirectApplicInterface() {};

protected:

  // Virtual function redefinitions

  //int derived_map_if(const Dakota::String& if_name);
  int derived_map_ac(const Dakota::String& ac_name);
  int derived_map_of(const Dakota::String& of_name);

private:

  // Data
    Teuchos::RCP<EpetraExt::ModelEvaluator> App;
    Teuchos::RCP<Epetra_Vector> model_p;
    Teuchos::RCP<Epetra_Vector> model_g;
    Teuchos::RCP<Epetra_MultiVector> model_dgdp;
    int numParameters;
    int numResponses;
};

} // namespace TriKota

#endif //TRIKOTA_DIRECTAPPLICINTERFACE
