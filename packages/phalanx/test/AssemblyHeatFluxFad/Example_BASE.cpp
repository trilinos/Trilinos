// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Phalanx_Evaluator_Factory.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Kokkos_Core.hpp"

// User defined objects
#include "Cell.hpp"
#include "Workset.hpp"
#include "FactoryTraits.hpp"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  GlobalMPISession mpi_session(&argc, &argv);

  try {

    Kokkos::initialize(argc,argv);
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    bool print_debug_info = false;

    {
      cout << "\nStarting MultiDimensionalArray EnergyFlux Example!" << endl;

      // Assume we have 102 cells on processor and can fit 20 cells in cache
      typedef PHX::Device::size_type size_type;
      const size_type num_local_cells = 102;
      const size_type workset_size = 20;
      // const size_type num_local_cells = 80000;
      // const size_type workset_size = 100;

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
      RCP< vector< RCP<Evaluator_TemplateManager<MyTraits> > > >  evaluators;
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

      // For Kokkos extended types (Sacado FAD) set derivtive array
      // size
      std::vector<PHX::index_size_type> derivative_dimensions;
      derivative_dimensions.push_back(8);
      fm.setKokkosExtendedDataTypeDimensions<MyTraits::Jacobian>(derivative_dimensions);

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
      fm.getFieldData<MyTraits::Residual,double,Cell,QuadPoint,Dim>(energy_flux);
      fm.getFieldData<MyTraits::Residual,double,Cell,QuadPoint>(source);

      // Get host arrays to copy from device back to host
      auto host_energy_flux = Kokkos::create_mirror_view(energy_flux.get_static_view());
      auto host_source = Kokkos::create_mirror_view(source.get_static_view());

      RCP<Time> eval_time = TimeMonitor::getNewTimer("Evaluation Time Residual");

      fm.preEvaluate<MyTraits::Residual>(NULL);
      {
	TimeMonitor t(*eval_time);

	for (std::size_t i = 0; i < worksets.size(); ++i) {
#ifdef PHX_ENABLE_KOKKOS_AMT
	  fm.evaluateFieldsTaskParallel<MyTraits::Residual>(worksets[i].num_cells,worksets[i]);
#else
	  fm.evaluateFields<MyTraits::Residual>(worksets[i]);
#endif	  
	  // Use values: in this example, move values into local host arrays
	  Kokkos::deep_copy(host_energy_flux,energy_flux.get_static_view());
	  Kokkos::deep_copy(host_source,source.get_static_view());

	  for (size_type cell = 0; cell < energy_flux.dimension(0); ++cell) {
	    for (size_type ip = 0; ip < energy_flux.dimension(1); ++ip) {
	      for (size_type dim = 0; dim < energy_flux.dimension(2); ++dim) {
		std::size_t index = cell * energy_flux.dimension(0) +
		  ip * energy_flux.dimension(1) + dim;
		local_energy_flux_at_qp[index] =  host_energy_flux(cell,ip,dim);
	      }
	    }
	  }
	  for (size_type cell = 0; cell < energy_flux.dimension(0); ++cell) {
	    for (size_type ip = 0; ip < energy_flux.dimension(1); ++ip) {
	      std::size_t index = cell * energy_flux.dimension(0) + ip;
	      local_source_at_qp[index] =  host_source(cell,ip);
	    }
	  }

	}

      }
      fm.postEvaluate<MyTraits::Residual>(NULL);

      // Test data retrieval
      cout << "Testing data members" << endl;
      Tag<double> d_var("Density", qp_scalar);
      MDField<double> den(d_var); 
      fm.getFieldData<MyTraits::Residual>(den);
      cout << "size of density = " << den.size() << ", should be " 
	   << d_var.dataLayout().size() << "." << endl;
      TEUCHOS_TEST_FOR_EXCEPTION(den.size() != d_var.dataLayout().size(),
			 std::runtime_error, 
			 "Returned arrays are not sized correctly!");

      cout << endl;

     //Jacobian:
     std::vector<JacScalarT>
        j_local_energy_flux_at_qp(num_local_cells * qp_vec->size());
      std::vector<JacScalarT>
        j_local_source_at_qp(num_local_cells * qp_scalar->size());
      // Fields we require
      MDField<JacScalarT,Cell,QuadPoint,Dim> j_energy_flux(j_energy_flux_tag);
      MDField<JacScalarT,Cell,QuadPoint> j_source(j_source_tag);
      fm.getFieldData<MyTraits::Jacobian,JacScalarT,Cell,QuadPoint,Dim>(j_energy_flux);
      fm.getFieldData<MyTraits::Jacobian,JacScalarT,Cell,QuadPoint>(j_source);
      RCP<Time> eval_time2 = TimeMonitor::getNewTimer("Evaluation Time Jacobian");

      auto host_j_energy_flux = Kokkos::create_mirror_view(j_energy_flux.get_static_view());
      auto host_j_source = Kokkos::create_mirror_view(j_source.get_static_view());

      fm.preEvaluate<MyTraits::Jacobian>(NULL);
      {
        TimeMonitor t(*eval_time2);

        for (std::size_t i = 0; i < worksets.size(); ++i) {
#ifdef PHX_ENABLE_KOKKOS_AMT
	  fm.evaluateFieldsTaskParallel<MyTraits::Jacobian>(worksets[i].num_cells,worksets[i]);
#else
          fm.evaluateFields<MyTraits::Jacobian>(worksets[i]);
#endif
	  Kokkos::deep_copy(host_j_energy_flux,j_energy_flux.get_static_view());
	  Kokkos::deep_copy(host_j_source,j_source.get_static_view());

           for (size_type cell = 0; cell < j_energy_flux.dimension(0); ++cell) {
            for (size_type ip = 0; ip < j_energy_flux.dimension(1); ++ip) {
              for (size_type dim = 0; dim < j_energy_flux.dimension(2); ++dim) {
                std::size_t index = cell * j_energy_flux.dimension(0) +
                  ip * j_energy_flux.dimension(1) + dim;
                j_local_energy_flux_at_qp[index] =  host_j_energy_flux(cell,ip,dim);
              }
            }
          }
          for (size_type cell = 0; cell < j_energy_flux.dimension(0); ++cell) {
            for (size_type ip = 0; ip < j_energy_flux.dimension(1); ++ip) {
              std::size_t index = cell * j_energy_flux.dimension(0) + ip;
              j_local_source_at_qp[index] =  host_j_source(cell,ip);
            }
          }


        }
       }
      fm.postEvaluate<MyTraits::Jacobian>(NULL);

      double scalability = 0.0;
      double parallelizability = 1.0;
      fm.analyzeGraph<MyTraits::Residual>(scalability,parallelizability);
      std::cout << "Residual Scalability = " << scalability << std::endl;
      std::cout << "Residual Parallelizability = " 
		<< parallelizability << std::endl;
      fm.analyzeGraph<MyTraits::Jacobian>(scalability,parallelizability);
      std::cout << "Residual Scalability = " << scalability << std::endl;
      std::cout << "Residual Parallelizability = " 
		<< parallelizability << std::endl;
    }

    Kokkos::finalize();
    
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
