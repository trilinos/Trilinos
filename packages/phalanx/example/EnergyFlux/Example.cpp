// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Cell)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Cell)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Node)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Node)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(QP)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(QP)

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
      cout << "\nStarting EnergyFlux Example!" << endl;

      // Assume we have 102 cells on processor and can fit 20 cells in cache
      const std::size_t num_local_cells = 102;
      const std::size_t workset_size = 20;

      RCP<DataLayout> qp = rcp(new MDALayout<Cell,QP>(workset_size, 4));
      RCP<DataLayout> node = rcp(new MDALayout<Cell,Node>(workset_size, 4));

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

      { // Constant Heat Capacity
	RCP<ParameterList> p = rcp(new ParameterList);
	int type = MyFactoryTraits<MyTraits>::id_constant;
	p->set<int>("Type", type);
	p->set<string>("Name", "Thermal Conductivity");
	p->set<double>("Value", 2.0);
	p->set< RCP<DataLayout> >("Data Layout", qp);
	evaluators_to_build["Thermal Conductivity"] = p;
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

      Tag< MyVector<ResScalarT> > energy_flux_tag("Energy_Flux", qp);
      fm.requireField<MyTraits::Residual>(energy_flux_tag);

      Tag<ResScalarT> source_tag("Nonlinear Source", qp);
      fm.requireField<MyTraits::Residual>(source_tag);

      // Request quantities to assemble JACOBIAN PDE operators
      typedef MyTraits::Jacobian::ScalarT JacScalarT;

      Tag< MyVector<JacScalarT> > j_energy_flux_tag("Energy_Flux", qp);
      fm.requireField<MyTraits::Jacobian>(j_energy_flux_tag);

      Tag<JacScalarT> j_source_tag("Nonlinear Source", qp);
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
      std::vector< MyVector<double> > 
	local_energy_flux_at_qp(num_local_cells * qp->size());
      std::vector<double> local_source_at_qp(num_local_cells * qp->size());

      // Fields we require
      Field< MyVector<double> > energy_flux(energy_flux_tag);
      Field<double> source(source_tag);
      fm.getFieldData< MyVector<double>,MyTraits::Residual >(energy_flux);
      fm.getFieldData<double,MyTraits::Residual>(source);

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
      Tag<double> d_var("Density", qp);
      Field<double> den(d_var); 
      fm.getFieldData<double,MyTraits::Residual>(den);
      cout << "size of density = " << den.size() << ", should be " 
	   << d_var.dataLayout().size() << "." << endl;
      TEST_FOR_EXCEPTION(den.size() != d_var.dataLayout().size(),
			 std::runtime_error, 
			 "Returned arrays are not sized correctly!");
      
      cout << endl;

      // ************************************************************
      // * Tests to make sure solution is correct
      // ************************************************************

      // Compare temperature fields, should be 2.0
      Field<double> temp("Temperature", node);
      fm.getFieldData<double,MyTraits::Residual>(temp);
      
      Field<double> temp_base("Temperature Baseline", node);
      ArrayRCP<double> temp_base_data = 
	arcp<double>(node->size());
      temp_base.setFieldData(temp_base_data);
      for (int i=0; i<temp_base.size(); ++i)
	temp_base[i] = 2.0;
      
      compareScalarFields(temp, temp_base, 1.0e-12);

      cout << endl;

      // Compare temperature gradient fields, should be 2.0
      Field< MyVector<double> > tg("Temperature Gradient", qp);
      fm.getFieldData<MyVector<double>,MyTraits::Residual>(tg);

      Field< MyVector<double> > 
	tg_base("Temperature Gradient Baseline", qp);
      ArrayRCP< MyVector<double> > tg_base_data = 
	arcp< MyVector<double> >(qp->size());
      tg_base.setFieldData(tg_base_data);
      for (int i=0; i<tg_base.size(); ++i)
	tg_base[i] = 2.0;
      
      compareVectorFields(tg, tg_base, 1.0e-12);

      cout << endl;

      // Compare energy flux fields, should be -16.0
      Field< MyVector<double> > ef("Energy_Flux", qp);
      fm.getFieldData<MyVector<double>,MyTraits::Residual>(ef);

      Field< MyVector<double> > ef_base("Energy_Flux Baseline", qp);
      ArrayRCP< MyVector<double> > ef_base_data = 
	arcp< MyVector<double> >(qp->size());
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
