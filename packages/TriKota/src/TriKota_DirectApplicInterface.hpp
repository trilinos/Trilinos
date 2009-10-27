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

//!  TriKota namespace
namespace TriKota {

/*! \brief Adapter class that transates from a Trilinos interface toa Dakota interface. 
  An object of this class IS a Dakota::DirectApplicInterface
  and wraps an EpetraExt::ModelEvaluator. It can then be passed in
  as the argument to the TriKota::Driver::run method.
*/
class DirectApplicInterface : public Dakota::DirectApplicInterface
{
public:

  //! Constructor that takes the Model Evaluator to wrap

   DirectApplicInterface(Dakota::ProblemDescDB& problem_db_,
                         const Teuchos::RCP<EpetraExt::ModelEvaluator> App_);

  ~DirectApplicInterface() {};

protected:

  //! Virtual function redefinition from Dakota::DirectApplicInterface
  int derived_map_ac(const Dakota::String& ac_name);

  //! Virtual function redefinition from Dakota::DirectApplicInterface
  int derived_map_of(const Dakota::String& of_name);

  //int derived_map_if(const Dakota::String& if_name);

private:

  // Data
    Teuchos::RCP<EpetraExt::ModelEvaluator> App;
    Teuchos::RCP<Epetra_Vector> model_p;
    Teuchos::RCP<Epetra_Vector> model_g;
    Teuchos::RCP<Epetra_MultiVector> model_dgdp;
    int numParameters;
    int numResponses;
    bool supportsSensitivities;
    EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation orientation;

};

} // namespace TriKota

#endif //TRIKOTA_DIRECTAPPLICINTERFACE
