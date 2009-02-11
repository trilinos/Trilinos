// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
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

// User defined objects
#include "Cell.hpp"
#include "Edge.hpp"
#include "MultiWorkset.hpp"
#include "Traits.hpp"
#include "MultiWorkset_FactoryTraits.hpp"

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
int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    bool print_debug_info = false;

    {
      cout << "\nStarting MultiWorkset Example!" << endl;

      RCP<DataLayout> qp = rcp(new FlatLayout("QP", 4, "CELL"));
      RCP<DataLayout> edge = rcp(new FlatLayout("Centroid", 1, "EDGE"));

      // Parser will build parameter list that determines the field
      // evaluators to build
      map<string, RCP<ParameterList> > evaluators_to_build;
      
      // *******************
      // Evaluators
      // *******************

      { // Constant Temperature for finite element
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_constant;
	p->set<int>("Type", type);
	p->set<string>("Name", "Temperature");
	p->set<double>("Value", 2.0);
	p->set< RCP<DataLayout> >("Data Layout", qp);
	evaluators_to_build["Finite Element Temperature"] = p;
      }
      
      { // Constant Temperature for finite volume
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_constant;
	p->set<int>("Type", type);
	p->set<string>("Name", "Temperature");
	p->set<double>("Value", 4.0);
	p->set< RCP<DataLayout> >("Data Layout", edge);
	evaluators_to_build["Finite Volume Temperature"] = p;
      }
      
      { // Density for finite element
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_density;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Data Layout", qp);
	evaluators_to_build["Finite Element Density"] = p;
      }

      { // Density for finite volume
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_density;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Data Layout", edge);
	evaluators_to_build["Finite Volume Density"] = p;
      }

      // Build Field Evaluators for each evaluation type
      EvaluatorFactory<MyTraits,MyFactoryTraits<MyTraits> > factory;
      RCP< vector< RCP<Evaluator_TemplateManager<MyTraits> > > > 
	evaluators;
      evaluators = factory.buildEvaluators(evaluators_to_build);
 
      // Create a FieldManager
      FieldManager<MyTraits> fm;

      // Register all Evaluators 
      registerEvaluators(evaluators, fm);

      // Request quantities to assemble RESIDUAL PDE operators
      {
	typedef MyTraits::Residual::ScalarT ResScalarT;
	
	Tag<ResScalarT> density_fe("Density", qp);
	fm.requireField<MyTraits::Residual>(density_fe);
	
	Tag<ResScalarT> density_fv("Density", edge);
	fm.requireField<MyTraits::Residual>(density_fv);
      }

      // Request quantities to assemble JACOBIAN PDE operators
      {
	typedef MyTraits::Jacobian::ScalarT JacScalarT;
	
	Tag<JacScalarT> j_density_fe("Density", qp);
	fm.requireField<MyTraits::Jacobian>(j_density_fe);

	Tag<JacScalarT> j_density_fv("Density", edge);
	fm.requireField<MyTraits::Jacobian>(j_density_fv);
      }

      // ***********
      // Workset construction
      // ***********

      const std::size_t num_cells_in_workset = 20;
      const std::size_t num_edges_in_workset = 25;
      std::map<std::string,std::size_t> workset_types;
      workset_types["CELL"] = num_cells_in_workset;
      workset_types["EDGE"] = num_edges_in_workset;

      RCP<Time> registration_time = 
	TimeMonitor::getNewTimer("Post Registration Setup Time");
      {
	TimeMonitor t(*registration_time);
	fm.postRegistrationSetup(workset_types);
      }

      if (print_debug_info)
	cout << fm << endl;
      
      // Create Workset information: Cells and EvalData objects
      std::vector<MyCell> cells(num_cells_in_workset);
      std::vector<MyEdge> edges(num_edges_in_workset);
      MyMultiWorkset workset;
      workset.cells = cells;
      workset.edges = edges;

      // ***********
      // Evaluation
      // ***********

      RCP<Time> eval_time = TimeMonitor::getNewTimer("Evaluation Time");
      {
	TimeMonitor t(*eval_time);
	fm.evaluateFields<MyTraits::Residual>(workset);
	fm.evaluateFields<MyTraits::Jacobian>(workset);
      }

      // ************************************************************
      // * Tests to make sure fields are correct
      // ************************************************************

      
      {// Compare temperature fe fields, should be 2.0
	cout << endl;

	Field<double> temp("Temperature", qp);
	fm.getFieldData<double,MyTraits::Residual>(temp);
	
	Field<double> temp_base("Temperature Baseline", qp);
	ArrayRCP<double> temp_base_data = 
	  arcp<double>(num_cells_in_workset * qp->size());
	temp_base.setFieldData(temp_base_data);
	for (int i=0; i<temp_base.size(); ++i)
	  temp_base[i] = 2.0;
	compareScalarFields(temp, temp_base, 1.0e-12);
	
	cout << endl;
      }

      {// Compare temperature fv fields, should be 4.0
	Field<double> temp("Temperature", edge);
	fm.getFieldData<double,MyTraits::Residual>(temp);
	
	
	Field<double> temp_base("Temperature Baseline", edge);
	ArrayRCP<double> temp_base_data = 
	  arcp<double>(num_edges_in_workset * edge->size());
	temp_base.setFieldData(temp_base_data);
	for (int i=0; i<temp_base.size(); ++i)
	  temp_base[i] = 4.0;
	compareScalarFields(temp, temp_base, 1.0e-12);
	
	cout << endl;
      }

      
      {  // Compare density fe fields, should be 2.0^2
	Field<double> d("Density", qp);
	fm.getFieldData<double,MyTraits::Residual>(d);
	
	Field<double> d_base("Density Baseline", qp);
	ArrayRCP<double> d_base_data = 
	  arcp<double>(num_cells_in_workset * qp->size());
	d_base.setFieldData(d_base_data);
	for (int i=0; i<d_base.size(); ++i)
	  d_base[i] = 4.0;
	compareScalarFields(d, d_base, 1.0e-12);
	
	cout << endl;
      }      

      {  // Compare density fv fields, should be 4.0^2
	Field<double> d("Density", edge);
	fm.getFieldData<double,MyTraits::Residual>(d);
	
	Field<double> d_base("Density Baseline", edge);
	ArrayRCP<double> d_base_data = 
	  arcp<double>(num_edges_in_workset * edge->size());
	d_base.setFieldData(d_base_data);
	for (int i=0; i<d_base.size(); ++i)
	  d_base[i] = 16.0;
	compareScalarFields(d, d_base, 1.0e-12);
	
	cout << endl;
      }


      // ************************************************************
      // EXAMPLE IS FINISHED!!!!
      // Below are some unit tests for full code coverage.
      // Please ignore.
      // ************************************************************

      cout << "Beginning Unit tests!" << endl;

      // Test post registration setup command for individual evaluation types
      // need a new field manager
      {
	// Create a FieldManager
	FieldManager<MyTraits> fm2;
	
	// Reuse original Evaluators 
	registerEvaluators(evaluators, fm2);

	{
	  typedef MyTraits::Residual::ScalarT ResScalarT;
	  
	  Tag<ResScalarT> density_fe("Density", qp);
	  fm2.requireField<MyTraits::Residual>(density_fe);
	  
	  Tag<ResScalarT> density_fv("Density", edge);
	  fm2.requireField<MyTraits::Residual>(density_fv);
	}
	
	{
	  typedef MyTraits::Jacobian::ScalarT JacScalarT;
	  
	  Tag<JacScalarT> j_density_fe("Density", qp);
	  fm2.requireField<MyTraits::Jacobian>(j_density_fe);
	  
	  Tag<JacScalarT> j_density_fv("Density", edge);
	  fm2.requireField<MyTraits::Jacobian>(j_density_fv);
	}
	
	fm2.postRegistrationSetupForType<MyTraits::Residual>(workset_types);
	fm2.postRegistrationSetupForType<MyTraits::Jacobian>(workset_types);
	
	
	RCP<Time> eval_time = TimeMonitor::getNewTimer("Evaluation Time");
	{
	  TimeMonitor t(*eval_time);
	  fm.evaluateFields<MyTraits::Residual>(workset);
	  fm.evaluateFields<MyTraits::Jacobian>(workset);
	}

	{// Compare temperature fe fields, should be 2.0
	  cout << endl;
	  
	  Field<double> temp("Temperature", qp);
	  fm2.getFieldData<double,MyTraits::Residual>(temp);
	  
	  Field<double> temp_base("Temperature Baseline", qp);
	  ArrayRCP<double> temp_base_data = 
	    arcp<double>(num_cells_in_workset * qp->size());
	  temp_base.setFieldData(temp_base_data);
	  for (int i=0; i<temp_base.size(); ++i)
	    temp_base[i] = 2.0;
	  compareScalarFields(temp, temp_base, 1.0e-12);
	
	  cout << endl;
	}
	
	{// Compare temperature fv fields, should be 4.0
	  Field<double> temp("Temperature", edge);
	  fm2.getFieldData<double,MyTraits::Residual>(temp);
	  
	  
	  Field<double> temp_base("Temperature Baseline", edge);
	  ArrayRCP<double> temp_base_data = 
	    arcp<double>(num_edges_in_workset * edge->size());
	  temp_base.setFieldData(temp_base_data);
	  for (int i=0; i<temp_base.size(); ++i)
	    temp_base[i] = 4.0;
	  compareScalarFields(temp, temp_base, 1.0e-12);
	  
	  cout << endl;
	}
	
	
	{  // Compare density fe fields, should be 2.0^2
	  Field<double> d("Density", qp);
	  fm2.getFieldData<double,MyTraits::Residual>(d);
	  
	  Field<double> d_base("Density Baseline", qp);
	  ArrayRCP<double> d_base_data = 
	    arcp<double>(num_cells_in_workset * qp->size());
	  d_base.setFieldData(d_base_data);
	  for (int i=0; i<d_base.size(); ++i)
	    d_base[i] = 4.0;
	  compareScalarFields(d, d_base, 1.0e-12);
	  
	  cout << endl;
	}      
	
	{  // Compare density fv fields, should be 4.0^2
	  Field<double> d("Density", edge);
	  fm2.getFieldData<double,MyTraits::Residual>(d);
	  
	  Field<double> d_base("Density Baseline", edge);
	  ArrayRCP<double> d_base_data = 
	    arcp<double>(num_edges_in_workset * edge->size());
	  d_base.setFieldData(d_base_data);
	  for (int i=0; i<d_base.size(); ++i)
	    d_base[i] = 16.0;
	  compareScalarFields(d, d_base, 1.0e-12);
	  
	  cout << endl;
	}

      }

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
