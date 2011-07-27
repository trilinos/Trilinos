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
#include "Teuchos_GlobalMPISession.hpp"

// User defined objects
#include "Cell.hpp"
#include "Workset.hpp"
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
int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  GlobalMPISession mpi_session(&argc, &argv);

  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    bool print_debug_info = false;

    {
      cout << "\nStarting MultiDimensionalArray EnergyFlux Example!" << endl;

      // Assume we have 102 cells on processor and can fit 20 cells in cache
      const std::size_t num_local_cells = 102;
      const std::size_t workset_size = 20;

      RCP<DataLayout> qp_scalar = 
	rcp(new MDALayout<Cell,Node>(workset_size,4));
      RCP<DataLayout> node_scalar = 
	rcp(new MDALayout<Cell,QuadPoint>(workset_size,4));

      RCP<DataLayout> qp_vec = 
	rcp(new MDALayout<Cell,Node,Dim>(workset_size,4,3));
      RCP<DataLayout> node_vec = 
	rcp(new MDALayout<Cell,QuadPoint,Dim>(workset_size,4,3));

      RCP<DataLayout> qp_mat = 
	rcp(new MDALayout<Cell,Node,Dim,Dim>(workset_size,4,3,3));
      RCP<DataLayout> node_mat = 
	rcp(new MDALayout<Cell,QuadPoint,Dim,Dim>(workset_size,4,3,3));

      // Parser will build parameter list that determines the field
      // evaluators to build
      map<string, RCP<ParameterList> > evaluators_to_build;
      
      { // Temperature
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_constant;
	p->set<int>("Type", type);
	p->set<string>("Name", "Temperature");
	p->set<double>("Value", 2.0);
	p->set< RCP<DataLayout> >("Data Layout", node_scalar);
	evaluators_to_build["DOF_Temperature"] = p;
      }
      { // Density
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_density;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Data Layout", qp_scalar);
	evaluators_to_build["Density"] = p;
      }

      { // Constant Thermal Conductivity
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_constant;
	p->set<int>("Type", type);
	p->set<string>("Name", "Thermal Conductivity");
	p->set<double>("Value", 2.0);
	p->set< RCP<DataLayout> >("Data Layout", qp_scalar);
	evaluators_to_build["Thermal Conductivity"] = p;
      }
      
      { // Nonlinear Source
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_nonlinearsource;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Data Layout", qp_scalar);
	evaluators_to_build["Nonlinear Source"] = p;
      }

      { // Fourier Energy Flux
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_fourier;
	p->set<int>("Type", type);
	p->set< RCP<DataLayout> >("Scalar Data Layout", qp_scalar);
	p->set< RCP<DataLayout> >("Vector Data Layout", qp_vec);
	evaluators_to_build["Energy Flux"] = p;
      }

      { // FE Interpolation
	RCP<ParameterList> p = rcp(new ParameterList);

	int type = MyFactoryTraits<MyTraits>::id_feinterpolation;
	p->set<int>("Type", type);

	p->set<string>("Node Variable Name", "Temperature");
	p->set<string>("QP Variable Name", "Temperature");
	p->set<string>("Gradient QP Variable Name", "Temperature Gradient");

	p->set< RCP<DataLayout> >("Node Data Layout", node_scalar);
	p->set< RCP<DataLayout> >("QP Scalar Data Layout", qp_scalar);
	p->set< RCP<DataLayout> >("QP Vector Data Layout", qp_vec);

	evaluators_to_build["FE Interpolation"] = p;
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
      typedef MyTraits::Residual::ScalarT ResScalarT;
      
      Tag<ResScalarT> energy_flux_tag("Energy_Flux", qp_vec);
      fm.requireField<MyTraits::Residual>(energy_flux_tag);
      
      Tag<ResScalarT> source_tag("Nonlinear Source", qp_scalar);
      fm.requireField<MyTraits::Residual>(source_tag);

      // Request quantities to assemble JACOBIAN PDE operators
      typedef MyTraits::Jacobian::ScalarT JacScalarT;
      
      Tag<JacScalarT> j_energy_flux_tag("Energy_Flux", qp_vec);
      fm.requireField<MyTraits::Jacobian>(j_energy_flux_tag);
      
      Tag<JacScalarT> j_source_tag("Nonlinear Source", qp_scalar);
      fm.requireField<MyTraits::Jacobian>(j_source_tag);

      RCP<Time> registration_time = 
	TimeMonitor::getNewTimer("Post Registration Setup Time");
      {
	TimeMonitor t(*registration_time);
	fm.postRegistrationSetup(NULL);
      }

      cout << fm << endl;

      fm.writeGraphvizFile<MyTraits::Residual>("graph_residual.dot");
      fm.writeGraphvizFile<MyTraits::Jacobian>("graph_jacobian.dot");
      fm.writeGraphvizFile("all_graph", ".dot");      

      // Create Workset information: Cells and EvalData objects
      std::vector<MyCell> cells(num_local_cells);
      for (std::size_t i = 0; i < cells.size(); ++i)
	cells[i].setLocalIndex(i);
      std::vector<MyWorkset> worksets;

      std::vector<MyCell>::iterator cell_it = cells.begin();
      std::size_t count = 0;
      MyWorkset w;
      w.local_offset = cell_it->localIndex();
      w.begin = cell_it;
      for (; cell_it != cells.end(); ++cell_it) {
	++count;
	std::vector<MyCell>::iterator next = cell_it;
	++next;
	
	if ( count == workset_size || next == cells.end()) {
	  w.end = next;
	  w.num_cells = count;
	  worksets.push_back(w);
	  count = 0;

	  if (next != cells.end()) {
	    w.local_offset = next->localIndex();
	    w.begin = next;
	  }
	}
      }
      
      if (print_debug_info) {
	for (std::size_t i = 0; i < worksets.size(); ++i) {
	  cout << "worksets[" << i << "]" << endl;
	  cout << "  num_cells =" << worksets[i].num_cells << endl;
	  cout << "  local_offset =" << worksets[i].local_offset << endl;
	  std::vector<MyCell>::iterator it = worksets[i].begin;
	  for (; it != worksets[i].end; ++it)
	    cout << "  cell_local_index =" << it->localIndex() << endl;
	}
	cout << endl;
      }

      // Create local vectors to hold flux and source at quad points
      std::vector<double> 
	local_energy_flux_at_qp(num_local_cells * qp_vec->size());
      std::vector<double> 
	local_source_at_qp(num_local_cells * qp_scalar->size());

      // Fields we require
      MDField<double,Cell,QuadPoint,Dim> energy_flux(energy_flux_tag);
      MDField<double,Cell,QuadPoint> source(source_tag);
      fm.getFieldData<double,MyTraits::Residual,Cell,QuadPoint,Dim>(energy_flux);
      fm.getFieldData<double,MyTraits::Residual,Cell,QuadPoint>(source);

      RCP<Time> eval_time = TimeMonitor::getNewTimer("Evaluation Time");

      fm.preEvaluate<MyTraits::Residual>(NULL);
      {
	TimeMonitor t(*eval_time);

	for (std::size_t i = 0; i < worksets.size(); ++i) {
	  fm.evaluateFields<MyTraits::Residual>(worksets[i]);
	  
	  // Use values: in this example, move values into local arrays
	  for (int j = 0; j < energy_flux.size(); ++j) {
	    std::size_t index = worksets[i].local_offset + j;
	    local_energy_flux_at_qp[index] =  energy_flux[j];
	  }
	  for (int j = 0; j < source.size(); ++j) {
	    std::size_t index = worksets[i].local_offset + j;
	    local_source_at_qp[index] =  source[j];
	  }

	}

      }
      fm.postEvaluate<MyTraits::Residual>(NULL);

      // Test data retrieval
      cout << "Testing data members" << endl;
      Tag<double> d_var("Density", qp_scalar);
      Field<double> den(d_var); 
      fm.getFieldData<double,MyTraits::Residual>(den);
      cout << "size of density = " << den.size() << ", should be " 
	   << d_var.dataLayout().size() << "." << endl;
      TEST_FOR_EXCEPTION(den.size() != d_var.dataLayout().size(),
			 std::runtime_error, 
			 "Returned arrays are not sized correctly!");

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
