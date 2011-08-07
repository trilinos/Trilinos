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

#include "CellData.hpp"
#include "Traits.hpp"
#include "Traits_ContiguousAllocator.hpp"
#include "FactoryTraits.hpp"

using namespace std;
using namespace Teuchos;
using namespace PHX;
  
SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Cell)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Cell)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Node)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Node)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(QP)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(QP)

// **********************************************************************
// FieldManager Builder Function
// **********************************************************************
template<typename Traits>
RCP< FieldManager<Traits> > buildFieldManager(std::size_t num_cells)
{
  RCP<DataLayout> qp = rcp(new MDALayout<Cell,QP>(num_cells,4));
  RCP<DataLayout> node = rcp(new MDALayout<Cell,Node>(num_cells,4));

  // Parser will build parameter list that determines the field
  // evaluators to build
  map<string, RCP<ParameterList> > evaluators_to_build;
      
  { // Temperature
    RCP<ParameterList> p = rcp(new ParameterList);
    int type = MyFactoryTraits<Traits>::id_constant;
    p->set<int>("Type", type);
    p->set<string>("Name", "Temperature");
    p->set<double>("Value", 2.0);
    p->set< RCP<DataLayout> >("Data Layout", node);
    evaluators_to_build["DOF_Temperature"] = p;
  }
  { // Density
    RCP<ParameterList> p = rcp(new ParameterList);
    int type = MyFactoryTraits<Traits>::id_density;
    p->set<int>("Type", type);
    p->set< RCP<DataLayout> >("Data Layout", qp);
    evaluators_to_build["Density"] = p;
  }
  
  { // Constant Diffusion Coefficient
    RCP<ParameterList> p = rcp(new ParameterList);
    int type = MyFactoryTraits<Traits>::id_constant;
    p->set<int>("Type", type);
    p->set<string>("Name", "Diffusion Coefficient");
    p->set<double>("Value", 2.0);
    p->set< RCP<DataLayout> >("Data Layout", qp);
    evaluators_to_build["Diffusion Coefficient"] = p;
  }
  
  { // Nonlinear Source
    RCP<ParameterList> p = rcp(new ParameterList);
    int type = MyFactoryTraits<Traits>::id_nonlinearsource;
    p->set<int>("Type", type);
    p->set< RCP<DataLayout> >("Data Layout", qp);
    evaluators_to_build["Nonlinear Source"] = p;
  }

  { // Fourier Energy Flux
    RCP<ParameterList> p = rcp(new ParameterList);
    int type = MyFactoryTraits<Traits>::id_fourier;
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
  EvaluatorFactory<Traits,MyFactoryTraits<Traits> > factory;
  RCP< vector< RCP<Evaluator_TemplateManager<Traits> > > > 
    providers;
  providers = factory.buildEvaluators(evaluators_to_build);
  
  
  // Request quantities to assemble PDE operators
  RCP< FieldManager<Traits> > vm = rcp(new FieldManager<Traits>);
  Tag<MyVector<double> > energy_flux("Energy_Flux", qp);
  vm->template requireField<typename Traits::Residual>(energy_flux);
  Tag<double> source("Nonlinear Source", qp);
  vm->template requireField<typename Traits::Residual>(source);
  
  // Register all Evaluators
  registerEvaluators(providers, *vm);
  
  return vm;

}

// **********************************************************************
// Main
// **********************************************************************

int main(int argc, char *argv[]) 
{
  GlobalMPISession mpi_session(&argc, &argv);

  try {
    
    RCP<Time> time_total = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*time_total);
    
    // WARNING: For timings, we should be flushing the cache in
    // between each evaluation loop, to eliminate cache reuse.  Not a
    // big deal as we really just wanted to make sure alignment is
    // correct.

    // 1 * 50 * 100 for checked in test
    // 4 * 1 * 5000000 for timing results
    const std::size_t num_samples = 1;
    const std::size_t num_eval_loops = 50;
    const std::size_t num_cells = 100;
    const std::size_t total_work = num_samples * num_eval_loops * num_cells;
    std::vector<CellData> cells(num_cells);

    TEST_FOR_EXCEPTION(total_work != 5000,
		       std::logic_error,
		       "Total work is not consistent!");

    RCP< FieldManager<MyTraits> > fmn = 
      buildFieldManager<MyTraits>(num_cells);

    RCP< FieldManager<MyCTraits> > fmc = 
      buildFieldManager<MyCTraits>(num_cells);
    
    RCP<Time> time_fmn_prs = 
      TimeMonitor::getNewTimer("NEW: Post Registration Setup Time");
    {
      TimeMonitor t(*time_fmn_prs);
      fmn->postRegistrationSetup(NULL);
    }
    
    RCP<Time> time_fmc_prs = 
      TimeMonitor::getNewTimer("CONTIGUOUS: Post Registration Setup Time");
    {
      TimeMonitor t(*time_fmc_prs);
      fmc->postRegistrationSetup(NULL);
    }

    RCP<Time> time_fmn = 
      TimeMonitor::getNewTimer("NEW: Evaluation Time");
    {
      TimeMonitor t(*time_fmn);
      for (std::size_t i=0; i < num_samples; ++i)
	for (std::size_t j=0; j < num_eval_loops; ++j)
	  fmn->evaluateFields<MyTraits::Residual>(cells);
    }

    RCP<Time> time_fmc = 
      TimeMonitor::getNewTimer("CONTIGUOUS: Evaluation Time");
    {
      TimeMonitor t(*time_fmc);
      for (std::size_t i=0; i < num_samples; ++i)
	for (std::size_t j=0; j < num_eval_loops; ++j)
	  fmc->evaluateFields<MyCTraits::Residual>(cells);
    }

    // *********************************************************************
    // Finished all testing
    // *********************************************************************
    std::cout << "\nTest passed!\n" << std::endl; 
    // *********************************************************************
    // *********************************************************************

    std::cout << num_samples << " X " 
	      << num_eval_loops << " X " 
	      << num_cells << " = " 
	      << total_work << std::endl;

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
