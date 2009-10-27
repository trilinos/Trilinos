
#include "Diagonal_ThyraROME_def.hpp"

#include "TriKota_Driver.hpp"
#include "TriKota_ThyraDirectApplicInterface.hpp"
#include "Thyra_VectorStdOps.hpp"


#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"


namespace {


double sqr(const double v)
{
  return v*v;
}


} // namespace



int main(int argc, char* argv[])
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::OSTab;

  bool success = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  const RCP<FancyOStream>
    out = VerboseObjectBase::getDefaultOStream();

  try {

    *out << "\nStarting TriKota Example!" << std::endl;
    
    // Construct driver with default file names
    TriKota::Driver dakota("dakota_conmin.in");
    
    // Construct a ModelEvaluator for your application with the
    // MPI_Comm chosen by Dakota. This example ModelEvaluator 
    // only takes an MPI_Comm to construct.

    const int num_p = 16;
    const RCP<TriKota::DiagonalROME<double> > thyraApp =
      TriKota::createModel<double>(num_p,5.0);

    // Construct a concrete Dakota interface with an EpetraExt::ModelEvaluator
    Teuchos::RCP<TriKota::ThyraDirectApplicInterface> trikota_interface =
      Teuchos::rcp(new TriKota::ThyraDirectApplicInterface(dakota.getProblemDescDB(), thyraApp), false);
    
    // Run the requested Dakota strategy using this interface
    dakota.run(trikota_interface.get());

    // Get the final solution and check it!

    if (mpiSession.getRank() == 0)  {
      Dakota::Variables finalVariables = dakota.getFinalSolution();
      Dakota::RealVector finalValues = finalVariables.continuous_variables();
    
      std::cout << "\nfinalValues =\n" << finalValues;

      const double errorTol = 1e-6;
      double finalError = 0.0;
      for (int i=0; i< num_p; i++) finalError += sqr(finalValues[0] - 2.0);
      finalError = std::sqrt(finalError);

      std::cout << "\nfinalError = "<<finalError<<"\n";
      
      if (finalError > errorTol) {
        std::cout << "\nError:  finalError > errorTol = "<<errorTol<<"\n";
        success = false;
      }
      if (success) std::cout << "\nEnd Result: TEST PASSED\n";
    }

    *out << std::flush;

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  if(success)
    *out << "\nEnd Result: TEST PASSED\n";
  else
    *out << "\nEnd Result: TEST FAILED\n";
    
  return ( success ? 0 : 1 );


}
