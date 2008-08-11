// @HEADER
// ************************************************************************
// 
//            Phalanx: A Partial Differential Equation Assembly 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

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

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template<typename ScalarT>
void compareScalarFields(PHX::Field<ScalarT>& field1, 
			 PHX::Field<ScalarT>& field2,
			 double tol)
{
  std::cout << "Comparing scalar fields\n" << field1.fieldTag() << "\n" 
	    << field2.fieldTag() << std::endl; 

  TEST_FOR_EXCEPTION(field1.size() != field2.size(), std::logic_error,
		     "Fields for comparison do not have the same size!");

  double error = 0.0;
  for (int i=0; i < field1.size(); ++i)
    error += fabs(field1[i]-field2[i]);
  
  TEST_FOR_EXCEPTION(error > tol, std::runtime_error,
		     "Fields are not equal in comparison!");

  std::cout << "Passed: " << error << " < " << tol << std::endl; 
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template<typename ScalarT>
void compareVectorFields(PHX::Field< MyVector<ScalarT> >& field1, 
			 PHX::Field< MyVector<ScalarT> >& field2,
			 double tol)
{
  std::cout << "Comparing vector fields\n" << field1.fieldTag() << "\n" 
	    << field2.fieldTag() << std::endl; 

  TEST_FOR_EXCEPTION(field1.size() != field2.size(), std::logic_error,
		     "Fields for comparison do not have the same size!");

  double error = 0.0;
  for (int i=0; i < field1.size(); ++i)
    for (int j=0; j < 3; ++j)
      error += fabs(field1[i][j]-field2[i][j]);
  
  TEST_FOR_EXCEPTION(error > tol, std::runtime_error,
		     "Fields are not equal in comparison!");

  std::cout << "Passed: " << error << " < " << tol << std::endl; 
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    {
      cout << "\nStarting EnergyFlux Example!" << endl;

      RCP<DataLayout> qp = rcp(new Generic("Q1_QP", 4));
      RCP<DataLayout> node = rcp(new Generic("Q1_NODE", 4));

      // Parser will build parameter list that determines the field
      // evaluators to build
      map<string, RCP<ParameterList> > evaluators_to_build;
      
      { // Temperature
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_constant;
	p->set<int>("Type", type);
	p->set<string>("Name", "Temperature");
	p->set<double>("Value", 2.0);
	p->set< RCP<DataLayout> >("Data Layout", node);
	evaluators_to_build["DOF_Temperature"] = p;
      }
      { // Density
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_density;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Data Layout", qp);
	evaluators_to_build["Density"] = p;
      }

      { // Constant Diffusion Coefficient
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_constant;
	p->set<int>("Type", type);
	p->set<string>("Name", "Diffusion Coefficient");
	p->set<double>("Value", 2.0);
	p->set< RCP<DataLayout> >("Data Layout", qp);
	evaluators_to_build["Diffusion Coefficient"] = p;
      }
      
      { // Nonlinear Source
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_nonlinearsource;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Data Layout", qp);
	evaluators_to_build["Nonlinear Source"] = p;
      }

      { // Fourier Energy Flux
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_fourier;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Data Layout", qp);
	evaluators_to_build["Energy Flux"] = p;
      }

      { // FE Interpolation
	RCP<ParameterList> p = rcp(new ParameterList);

	int type = MyFactoryTraits<MyTraits>::id_feinterpolation;
	p->set<int>("Type", type);

	p->set<string>("Node Variable Name", "Temperature");
	p->set<string>("QP Variable Name", "Temperature");
	p->set<string>("Gradient QP Variable Name", "Temperature Gradient");

	p->set< RCP<DataLayout> >("Node Data Layout", node);
	p->set< RCP<DataLayout> >("QP Data Layout", qp);

	evaluators_to_build["FE Interpolation"] = p;
      }

      // Build Field Evaluators
      EvaluatorFactory<MyTraits,MyFactoryTraits<MyTraits> > factory;
      RCP< vector< RCP<Evaluator_TemplateManager<MyTraits> > > > 
	evaluators;
      evaluators = factory.buildEvaluators(evaluators_to_build);
 
      // Create a FieldManager
      FieldManager<MyTraits> vm;

      // Request quantities to assemble RESIDUAL PDE operators
      {
	typedef MyTraits::Residual::ScalarT ResScalarT;
	Tag< MyVector<ResScalarT> > energy_flux("Energy_Flux", qp);
	vm.requireField<MyTraits::Residual>(energy_flux);
	Tag<ResScalarT> source("Nonlinear Source", qp);
	vm.requireField<MyTraits::Residual>(source);
      }

      // Request quantities to assemble JACOBIAN PDE operators
      {
	typedef MyTraits::Jacobian::ScalarT JacScalarT;
	Tag< MyVector<JacScalarT> > energy_flux("Energy_Flux", qp);
	vm.requireField<MyTraits::Jacobian>(energy_flux);
	Tag<JacScalarT> source("Nonlinear Source", qp);
	vm.requireField<MyTraits::Jacobian>(source);
      }

      // Register all Evaluators 
      registerEvaluators(evaluators, vm);

      const std::size_t num_cells = 10;
      const std::size_t num_eval_loops = 1;

      RCP<Time> registration_time = 
	TimeMonitor::getNewTimer("Post Registration Setup Time");
      {
	TimeMonitor t(*registration_time);
	vm.postRegistrationSetup(num_cells);
      }

      cout << vm << endl;
      
      std::vector<CellData> cells(num_cells);

      RCP<Time> eval_time = TimeMonitor::getNewTimer("Evaluation Time");

      vm.preEvaluate<MyTraits::Residual>(NULL);
      {
	TimeMonitor t(*eval_time);
	for (std::size_t i=0; i < num_eval_loops; ++i)
	  vm.evaluateFields<MyTraits::Residual>(cells);
      }
      vm.postEvaluate<MyTraits::Residual>(NULL);

      // Test data retrieval
      cout << "Testing data members" << endl;
      Tag<double> d_var("Density", qp);
      Field<double> den(d_var); 
      vm.getFieldData<double,MyTraits::Residual>(den);
      cout << "size of density = " << den.size() << ", should be " 
	   << num_cells * d_var.dataLayout().size() << "." << endl;
      TEST_FOR_EXCEPTION(den.size() != static_cast<Teuchos::ArrayRCP<double>::Ordinal>(num_cells * d_var.dataLayout().size()),
			 std::runtime_error, 
			 "Returned arrays are not sized correctly!");
      
      
      cout << endl;

      // Compare temperature fields, should be 2.0
      Field<double> temp("Temperature", node);
      vm.getFieldData<double,MyTraits::Residual>(temp);
      
      Field<double> temp_base("Temperature Baseline", node);
      ArrayRCP<double> temp_base_data = 
	arcp<double>(num_cells * node->size());
      temp_base.setFieldData(temp_base_data);
      for (int i=0; i<temp_base.size(); ++i)
	temp_base[i] = 2.0;
      
      compareScalarFields(temp, temp_base, 1.0e-12);

      cout << endl;

      // Compare temperature gradient fields, should be 2.0
      Field< MyVector<double> > tg("Temperature Gradient", qp);
      vm.getFieldData<MyVector<double>,MyTraits::Residual>(tg);

      Field< MyVector<double> > 
	tg_base("Temperature Gradient Baseline", qp);
      ArrayRCP< MyVector<double> > tg_base_data = 
	arcp< MyVector<double> >(num_cells * qp->size());
      tg_base.setFieldData(tg_base_data);
      for (int i=0; i<tg_base.size(); ++i)
	tg_base[i] = 2.0;
      
      compareVectorFields(tg, tg_base, 1.0e-12);

      cout << endl;

      // Compare energy flux fields, should be -16.0
      Field< MyVector<double> > ef("Energy_Flux", qp);
      vm.getFieldData<MyVector<double>,MyTraits::Residual>(ef);

      Field< MyVector<double> > ef_base("Energy_Flux Baseline", qp);
      ArrayRCP< MyVector<double> > ef_base_data = 
	arcp< MyVector<double> >(num_cells * qp->size());
      ef_base.setFieldData(ef_base_data);
      for (int i=0; i<ef_base.size(); ++i)
	ef_base[i] = -16.0;
      
      compareVectorFields(ef, ef_base, 1.0e-12);

      cout << endl;

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

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
