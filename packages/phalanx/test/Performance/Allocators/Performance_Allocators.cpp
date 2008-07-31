#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "CellData.hpp"
#include "Traits.hpp"
#include "Traits_ContiguousAllocator.hpp"
#include "FactoryTraits.hpp"

using namespace std;
using namespace Teuchos;
using namespace PHX;
  
// **********************************************************************
// FieldManager Builder Function
// **********************************************************************
template<typename Traits>
RCP< FieldManager<Traits> > buildFieldManager()
{
  RCP<DataLayout> scalar_qp = 
    rcp(new Generic<typename Traits::MY_SCALAR>("Q1_QP", 4));
  RCP<DataLayout> vector_qp = 
    rcp(new Generic<typename Traits::MY_VECTOR>("Q1_QP", 4));
  RCP<DataLayout> scalar_node = 
    rcp(new Generic<typename Traits::MY_SCALAR>("Q1_NODE", 4));

  // Parser will build parameter list that determines the field
  // evaluators to build
  map<string, RCP<ParameterList> > evaluators_to_build;
      
  { // Temperature
    RCP<ParameterList> p = rcp(new ParameterList);
    int type = MyFactoryTraits<Traits>::id_constant;
    p->set<int>("Type", type);
    p->set<string>("Name", "Temperature");
    p->set<double>("Value", 2.0);
    p->set< RCP<DataLayout> >("Data Layout", scalar_node);
    evaluators_to_build["DOF_Temperature"] = p;
  }
  { // Density
    RCP<ParameterList> p = rcp(new ParameterList);
    int type = MyFactoryTraits<Traits>::id_density;
    p->set<int>("Type", type);
    p->set< RCP<DataLayout> >("Data Layout", scalar_qp);
    evaluators_to_build["Density"] = p;
  }
  
  { // Constant Diffusion Coefficient
    RCP<ParameterList> p = rcp(new ParameterList);
    int type = MyFactoryTraits<Traits>::id_constant;
    p->set<int>("Type", type);
    p->set<string>("Name", "Diffusion Coefficient");
    p->set<double>("Value", 2.0);
    p->set< RCP<DataLayout> >("Data Layout", scalar_qp);
    evaluators_to_build["Diffusion Coefficient"] = p;
  }
  
  { // Nonlinear Source
    RCP<ParameterList> p = rcp(new ParameterList);
    int type = MyFactoryTraits<Traits>::id_nonlinearsource;
    p->set<int>("Type", type);
    p->set< RCP<DataLayout> >("Data Layout", scalar_qp);
    evaluators_to_build["Nonlinear Source"] = p;
  }

  { // Fourier Energy Flux
    RCP<ParameterList> p = rcp(new ParameterList);
    int type = MyFactoryTraits<Traits>::id_fourier;
    p->set<int>("Type", type);
    p->set< RCP<DataLayout> >("Scalar Data Layout", scalar_qp);
    p->set< RCP<DataLayout> >("Vector Data Layout", vector_qp);
    evaluators_to_build["Energy Flux"] = p;
  }
  
  { // FE Interpolation
    RCP<ParameterList> p = rcp(new ParameterList);
    
    int type = MyFactoryTraits<MyTraits>::id_feinterpolation;
    p->set<int>("Type", type);
    
    p->set<string>("Node Variable Name", "Temperature");
    p->set<string>("QP Variable Name", "Temperature");
    p->set<string>("Gradient QP Variable Name", "Temperature Gradient");
    
    p->set< RCP<DataLayout> >("Node Data Layout", scalar_node);
    p->set< RCP<DataLayout> >("QP Data Layout", scalar_qp);
    p->set< RCP<DataLayout> >("Gradient QP Data Layout", vector_qp);
    
    evaluators_to_build["FE Interpolation"] = p;
  }
  
  // Build Field Evaluators
  EvaluatorFactory<Traits,MyFactoryTraits<Traits> > factory;
  RCP< vector< RCP<Evaluator_TemplateManager<Traits> > > > 
    providers;
  providers = factory.buildEvaluators(evaluators_to_build);
  
  
  // Request quantities to assemble PDE operators
  RCP< FieldManager<Traits> > vm = rcp(new FieldManager<Traits>);
  FieldTag energy_flux("Energy_Flux", vector_qp);
  vm->requireFieldForAllTypes(energy_flux);
  FieldTag source("Nonlinear Source", scalar_qp);
  vm->requireFieldForAllTypes(source);
  
  // Register all Evaluators
  registerEvaluators(providers, *vm);
  
  return vm;

}

// **********************************************************************
// Main
// **********************************************************************

int main(int argc, char *argv[]) 
{
  try {
    
    RCP<Time> time_total = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*time_total);
    
    RCP< FieldManager<MyTraits> > fmn = buildFieldManager<MyTraits>();
    RCP< FieldManager<MyCTraits> > fmc = buildFieldManager<MyCTraits>();
    
    const std::size_t num_eval_loops = 50;
    const std::size_t num_cells = 100;
    std::vector<CellData> cells(num_cells);

    RCP<Time> time_fmn_prs = 
      TimeMonitor::getNewTimer("NEW: Post Registration Setup Time");
    {
      TimeMonitor t(*time_fmn_prs);
      fmn->postRegistrationSetup(num_cells);
    }
    
    RCP<Time> time_fmc_prs = 
      TimeMonitor::getNewTimer("CONTIGUOUS: Post Registration Setup Time");
    {
      TimeMonitor t(*time_fmc_prs);
      fmc->postRegistrationSetup(num_cells);
    }

    RCP<Time> time_fmn = 
      TimeMonitor::getNewTimer("NEW: Evaluation Time");
    {
      TimeMonitor t(*time_fmn);
      for (std::size_t i=0; i < num_eval_loops; ++i)
	fmn->evaluateFields<MyTraits::Residual>(cells);
    }

    RCP<Time> time_fmc = 
      TimeMonitor::getNewTimer("CONTIGUOUS: Evaluation Time");
    {
      TimeMonitor t(*time_fmc);
      for (std::size_t i=0; i < num_eval_loops; ++i)
	fmc->evaluateFields<MyCTraits::Residual>(cells);
    }

    // *********************************************************************
    // Finished all testing
    // *********************************************************************
    std::cout << "\nTest passed!\n" << std::endl; 
    // *********************************************************************
    // *********************************************************************

  }
  catch (const std::exception& e) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Exception Caught!" << endl;
    std::cout << "Error message is below\n " << e.what() << endl;
    std::cout << "************************************************" << endl;
  }
  catch (...) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Unknown Exception Caught!" << endl;
    std::cout << "************************************************" << endl;
  }

  TimeMonitor::summarize();
    
  return 0;
}
