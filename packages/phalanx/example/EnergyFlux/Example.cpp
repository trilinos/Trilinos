#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "CellData.hpp"
#include "Traits.hpp"
#include "FactoryTraits.hpp"

int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    // *********************************************************************
    // Start of FieldManager testing
    // *********************************************************************
    {
      cout << "\nStarting FieldManager Testing !" << endl;

      RCP<DataLayout> scalar_qp = 
	rcp(new Generic<MyTraits::MY_SCALAR>("Q1_QP", 4));
      RCP<DataLayout> vector_qp = 
	rcp(new Generic<MyTraits::MY_VECTOR>("Q1_QP", 4));
      RCP<DataLayout> scalar_node = 
	rcp(new Generic<MyTraits::MY_SCALAR>("Q1_NODE", 4));

      // Parser will build parameter list that determines the field
      // evaluators to build
      map<string, RCP<ParameterList> > evaluators_to_build;
      
      { // Temperature
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_constant;
	p->set<int>("Type", type);
	p->set<string>("Name", "Temperature");
	p->set<double>("Value", 2.0);
	p->set< RCP<DataLayout> >("Data Layout", scalar_node);
	evaluators_to_build["DOF_Temperature"] = p;
      }
      { // Density
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_density;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Data Layout", scalar_qp);
	evaluators_to_build["Density"] = p;
      }

      { // Constant Diffusion Coefficient
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_constant;
	p->set<int>("Type", type);
	p->set<string>("Name", "Diffusion Coefficient");
	p->set<double>("Value", 2.0);
	p->set< RCP<DataLayout> >("Data Layout", scalar_qp);
	evaluators_to_build["Diffusion Coefficient"] = p;
      }
      
      { // Nonlinear Source
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_nonlinearsource;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Data Layout", scalar_qp);
	evaluators_to_build["Nonlinear Source"] = p;
      }

      { // Fourier Energy Flux
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_fourier;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Scalar Data Layout", scalar_qp);
	p->set< RCP<DataLayout> >("Vector Data Layout", vector_qp);
	evaluators_to_build["Energy Flux"] = p;
      }

      { // FE Interpolation
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_feinterpolation;
	p->set<int>("Type", type);
	RCP< vector<FieldTag> > dof_node = rcp(new vector<FieldTag>);
	dof_node->push_back(FieldTag("Temperature", scalar_node));
	RCP< vector<FieldTag> > dof_qp = rcp(new vector<FieldTag>);
	dof_qp->push_back(FieldTag("Temperature", scalar_qp));
	RCP< vector<FieldTag> > grad_dof_qp = rcp(new vector<FieldTag>);
	grad_dof_qp->push_back(FieldTag("Temperature Gradient", vector_qp));

	p->set< RCP< vector<FieldTag> > >("Scalar Node", dof_node);
	p->set< RCP< vector<FieldTag> > >("Scalar QP", dof_qp); 
	p->set< RCP< vector<FieldTag> > >("Grad Scalar Node", dof_node);
	p->set< RCP< vector<FieldTag> > >("Grad Vector QP", grad_dof_qp);
	evaluators_to_build["FE Interpolation"] = p;
      }

      // Build Field Evaluators
      FieldEvaluatorFactory<MyTraits,MyFactoryTraits<MyTraits> > factory;
      RCP< vector< RCP<FieldEvaluator_TemplateManager<MyTraits> > > > 
	providers;
      providers = factory.buildFieldEvaluators(evaluators_to_build);
 
          
      // Request quantities to assemble PDE operators
      FieldManager<MyTraits> vm;
      FieldTag energy_flux("Energy_Flux", vector_qp);
      vm.requireFieldForAllTypes(energy_flux);
      FieldTag source("Nonlinear Source", scalar_qp);
      vm.requireFieldForAllTypes(source);
      
      // Register all FieldEvaluators
      
      // Loop over each provider template manager
      vector< RCP<FieldEvaluator_TemplateManager<MyTraits> > >::iterator 
	tm = providers->begin();
      for (; tm != providers->end(); ++tm) {
	
	// Loop over Scalar Types
	PHX::FieldManager<MyTraits>::iterator vmit = vm.begin();
	FieldEvaluator_TemplateManager<MyTraits>::iterator vpit = 
	  (*tm)->begin();
	for (; vpit != (*tm)->end(); ++vpit) {
	  RCP<PHX::FieldEvaluator<MyTraits> > vp =
	    rcp_dynamic_cast<PHX::FieldEvaluator<MyTraits> >(vpit.rcp());
	  vm.registerEvaluatorForScalarType(vmit, vp);
	  ++vmit;
	}
	
      }

      const std::size_t num_cells = 10;
      const std::size_t num_eval_loops = 1;

      RCP<Time> registration_time = TimeMonitor::getNewTimer("Post Registration Setup Time");
      {
	TimeMonitor t(*registration_time);
	vm.postRegistrationSetup(num_cells);
      }

      cout << vm << endl;
      
      std::vector<CellData> cells(num_cells);

      RCP<Time> eval_time = TimeMonitor::getNewTimer("Evaluation Time");

      vm.preEvaluate<double>(NULL);
      {
	TimeMonitor t(*eval_time);
	for (std::size_t i=0; i < num_eval_loops; ++i)
	  vm.evaluateFields<double>(cells);
      }
      vm.postEvaluate<double>(NULL);

      // Test data retrieval
      cout << "Testing data members" << endl;
      FieldTag d_var("Density", scalar_qp);
      Field<double> den(d_var); 
      vm.setFieldData(den);
      cout << "size of density = " << den.size() << ", should be " 
	   << num_cells * d_var.dataLayout()->size() << "." << endl;
      TEST_FOR_EXCEPTION(den.size() != static_cast<Teuchos::ArrayRCP<double>::Ordinal>(num_cells * d_var.dataLayout()->size()),
			 std::runtime_error, 
			 "Returned arrays are not sized correctly!");
      
      
      //Field<double> temp("Temperature", scalar_node);
      //vm.setFieldData(temp);
      //for (int i=0; i < temp.size(); ++i)
      //cout << "temperature_node[" << i << "] = " << temp[i] << endl;
    }
    
    // *********************************************************************
    // Finished all testing
    // *********************************************************************
    std::cout << "\nRun has completed successfully!\n" << std::endl; 
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
