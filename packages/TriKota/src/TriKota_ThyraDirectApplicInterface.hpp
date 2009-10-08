#ifndef TRIKOTA_DIRECTAPPLICINTERFACE
#define TRIKOTA_DIRECTAPPLICINTERFACE

#include "DirectApplicInterface.H"
#include "CommandLineHandler.H"
#include "DakotaStrategy.H"
#include "DakotaModel.H"
#include "ParallelLibrary.H"
#include "ProblemDescDB.H"

#include "Thyra_ModelEvaluatorDefaultBase.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

//!  TriKota namespace
namespace TriKota {

/*! \brief Adapter class that transates from a Trilinos interface toa Dakota interface. 
  An object of this class IS a Dakota::DirectApplicInterface
  and wraps an EpetraExt::ModelEvaluator. It can then be passed in
  as the argument to the TriKota::Driver::run method.
*/
class ThyraDirectApplicInterface : public Dakota::DirectApplicInterface
{
public:

  //! Constructor that takes the Model Evaluator to wrap

   ThyraDirectApplicInterface(Dakota::ProblemDescDB& problem_db_,
                         const Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > App_);
                         //const Teuchos::RCP<Thyra::ModelEvaluator<double> > App_);

  ~ThyraDirectApplicInterface() {};

protected:

  //! Virtual function redefinition from Dakota::DirectApplicInterface
  int derived_map_ac(const Dakota::String& ac_name);

  //! Virtual function redefinition from Dakota::DirectApplicInterface
  int derived_map_of(const Dakota::String& of_name);

  //int derived_map_if(const Dakota::String& if_name);

private:

  // Data
    Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > App;
    Teuchos::RCP<Thyra::VectorBase<double> > model_p;
    Teuchos::RCP<Thyra::VectorBase<double> > model_g;
    Teuchos::RCP<Thyra::MultiVectorBase<double> > model_dgdp;
    int numParameters;
    int numResponses;
};

} // namespace TriKota

#endif //TRIKOTA_DIRECTAPPLICINTERFACE
