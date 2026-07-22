// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Phalanx_DataLayout_DynamicLayout.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_Field.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Mesh.hpp"
#include "WorksetBuilder.hpp"
#include "LinearObjectFactory.hpp"
#include "CommandLineParser.hpp"
#include "PrintValues.hpp"
#include "MyTraits.hpp"
#include "Constant.hpp"
#include "GatherSolution.hpp"
#include "ProjectValueToQP.hpp"
#include "ProjectGradientToQP.hpp"
#include "IntegrateDiffusionTerm.hpp"
#include "IntegrateSourceTerm.hpp"
#include "ZeroContributedField.hpp"
#include "ScatterResidual.hpp"

#include <sstream>
#include <vector>

// User defined objects
int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  using PHX::MyTraits;
  using Residual = PHX::MyTraits::Residual;
  using Jacobian = PHX::MyTraits::Jacobian;
  
  Teuchos::GlobalMPISession mpi_session(&argc,&argv);
  int success = 0;

  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    Kokkos::initialize(argc,argv);
    PHX::exec_space().print_configuration(std::cout);
    
    // *********************************************************
    // * Build the Finite Element data structures
    // *********************************************************

    // Create the mesh
    phx_example::CommandLineParser p(argc,argv);
    const int nx = p.nx();
    const int ny = p.ny();
    const int nz = p.nz();
    const double lx = p.lx();
    const double ly = p.ly();
    const double lz = p.lz();
    const int num_equations = p.numEquations();
    const int workset_size = p.worksetSize();
    RCP<phx_example::Mesh> mesh = rcp(new phx_example::Mesh(nx, ny, nz, lx, ly, lz));
    std::vector<Workset> worksets;
    {
      WorksetBuilder builder;
      builder.buildWorksets(workset_size,*mesh,worksets);
    }

    phx_example::LinearObjectFactory lof(mesh->getNumNodes(),
                                         num_equations,
                                         mesh->getGlobalIndices());

    // Print statistics
    {
      std::cout << "Number of Elements = " << mesh->getNumElements() << std::endl;
      std::cout << "Number of Nodes = " << mesh->getNumNodes() << std::endl;
      std::cout << "Number of equations = " << num_equations << std::endl;
      std::cout << "Number of DOFs = " << lof.getNumDOFs() << std::endl;
      std::cout << "Matrix Size = " << lof.getMatrixSize() << std::endl;
      std::cout << "Workset Size = " << workset_size << std::endl;
      std::cout << "Number of Worksets = " << worksets.size() << std::endl;
    }
        
    // We will size the layouts later 
    RCP<PHX::DataLayout> qp_layout = rcp(new Layout("DL<CELL,QP>"));
    RCP<PHX::DataLayout> grad_qp_layout = rcp(new Layout("DL<CELL,QP,DIM>"));
    RCP<PHX::DataLayout> basis_layout = rcp(new Layout("DL<CELL,BASIS>"));
    RCP<PHX::DataLayout> scatter_layout = rcp(new Layout("DL<SCATTER>"));
    RCP<PHX::DataLayout> src_layout = rcp(new Layout("DL<CELL,QP>"));

    // For unit testing, force a user defined layout that is not optimal on every device type
    Teuchos::rcp_dynamic_cast<Layout>(grad_qp_layout,true)->setKokkosLayout(PHX::DataLayout::KokkosLayoutType::Left);
    Teuchos::rcp_dynamic_cast<Layout>(src_layout,true)->setKokkosLayout(PHX::DataLayout::KokkosLayoutType::Left);

    PHX::FieldManager<MyTraits> fm;
    {
      std::vector<PHX::index_size_type> derivative_dimensions;
      derivative_dimensions.push_back(8 * num_equations);
      fm.setKokkosExtendedDataTypeDimensions<MyTraits::Jacobian>(derivative_dimensions);
    }
    
    // Gather DOFs
    Kokkos::View<double*,PHX::Device> x = lof.createSolutionVector("x"); // solution
    for (int eq=0; eq < num_equations; ++eq) {
      std::stringstream s;
      s << "equation_" << eq;
      RCP<GatherSolution<Residual,MyTraits>> r =
        rcp(new GatherSolution<Residual,MyTraits>(s.str(),basis_layout,num_equations,eq,x));
      fm.registerEvaluator<Residual>(r);
      RCP<GatherSolution<Jacobian,MyTraits>> j =
        rcp(new GatherSolution<Jacobian,MyTraits>(s.str(),basis_layout,num_equations,eq,x));
      fm.registerEvaluator<Jacobian>(j);
    }

    // Project DOFs to qps
    for (int eq=0; eq < num_equations; ++eq) {
      std::stringstream s;
      s << "equation_" << eq;
      RCP<ProjectValueToQP<Residual,MyTraits>> r =
        rcp(new ProjectValueToQP<Residual,MyTraits>(s.str(),basis_layout,qp_layout));
      fm.registerEvaluator<Residual>(r);
      RCP<ProjectValueToQP<Jacobian,MyTraits>> j =
        rcp(new ProjectValueToQP<Jacobian,MyTraits>(s.str(),basis_layout,qp_layout));
      fm.registerEvaluator<Jacobian>(j);
    }

    // Project DOF Gradients to qps
    for (int eq=0; eq < num_equations; ++eq) {
      std::stringstream s;
      s << "equation_" << eq;
      RCP<ProjectGradientToQP<Residual,MyTraits>> r =
        rcp(new ProjectGradientToQP<Residual,MyTraits>(s.str(),basis_layout,grad_qp_layout));
      fm.registerEvaluator<Residual>(r);
      RCP<ProjectGradientToQP<Jacobian,MyTraits>> j =
        rcp(new ProjectGradientToQP<Jacobian,MyTraits>(s.str(),basis_layout,grad_qp_layout));
      fm.registerEvaluator<Jacobian>(j);
    }

    // source term
    {
      RCP<Constant<Residual,MyTraits>> r = rcp(new Constant<Residual,MyTraits>("s_dot",src_layout,5.0));
      fm.registerEvaluator<Residual>(r);
      RCP<Constant<Jacobian,MyTraits>> j = rcp(new Constant<Jacobian,MyTraits>("s_dot",src_layout,5.0));
      fm.registerEvaluator<Jacobian>(j);
    }
    
    // Zero contributed field
    for (int eq=0; eq < num_equations; ++eq) {
      std::stringstream s;
      s << "residual_" << eq;

      RCP<ZeroContributedField<Residual,MyTraits>> r =
        rcp(new ZeroContributedField<Residual,MyTraits>(s.str(),basis_layout));
      fm.registerEvaluator<Residual>(r);

      RCP<ZeroContributedField<Jacobian,MyTraits>> j =
        rcp(new ZeroContributedField<Jacobian,MyTraits>(s.str(),basis_layout));

      fm.registerEvaluator<Jacobian>(j);
    }
    
    // Integrate source term
    for (int eq=0; eq < num_equations; ++eq) {
      std::stringstream s;
      s << "s_dot";
      std::stringstream s_r;
      s_r << "residual_" << eq;
      
      RCP<IntegrateSourceTerm<Residual,MyTraits>> r =
        rcp(new IntegrateSourceTerm<Residual,MyTraits>(s.str(),src_layout,s_r.str(),basis_layout));

      fm.registerEvaluator<Residual>(r);
      
      RCP<IntegrateSourceTerm<Jacobian,MyTraits>> j =
        rcp(new IntegrateSourceTerm<Jacobian,MyTraits>(s.str(),src_layout,s_r.str(),basis_layout));

      fm.registerEvaluator<Jacobian>(j);
    }

    // Integrate diffusion term
    for (int eq=0; eq < num_equations; ++eq) {
      std::stringstream s;
      s << "equation_" << eq;
      std::stringstream s_r;
      s_r << "residual_" << eq;
      
      RCP<IntegrateDiffusionTerm<Residual,MyTraits>> r =
        rcp(new IntegrateDiffusionTerm<Residual,MyTraits>(s.str(),grad_qp_layout,s_r.str(),basis_layout));

      fm.registerEvaluator<Residual>(r);
      
      RCP<IntegrateDiffusionTerm<Jacobian,MyTraits>> j =
        rcp(new IntegrateDiffusionTerm<Jacobian,MyTraits>(s.str(),grad_qp_layout,s_r.str(),basis_layout));

      fm.registerEvaluator<Jacobian>(j);
    }

    // Scatter DOFs
    Kokkos::View<double*,PHX::Device> f = lof.createSolutionVector("global_residual"); // residual
    KokkosSparse::CrsMatrix<double,int,PHX::Device> J = lof.createJacobianMatrix("global_jacobian"); // Jacobian
    for (int eq=0; eq < num_equations; ++eq) {
      std::stringstream s;
      s << "residual_" << eq;
      
      RCP<PHX::FieldTag> scatter_tag_r = rcp(new Tag<MyTraits::Residual::ScalarT>(s.str(),scatter_layout));
      RCP<ScatterResidual<Residual,MyTraits>> r =
        rcp(new ScatterResidual<Residual,MyTraits>(scatter_tag_r,s.str(),basis_layout,eq,num_equations,f));
      fm.registerEvaluator<Residual>(r);
 
      RCP<PHX::FieldTag> scatter_tag_j = rcp(new Tag<MyTraits::Jacobian::ScalarT>(s.str(),scatter_layout));
      RCP<ScatterResidual<Jacobian,MyTraits>> j =
        rcp(new ScatterResidual<Jacobian,MyTraits>(scatter_tag_j,s.str(),basis_layout,eq,num_equations,f,J));
      fm.registerEvaluator<Jacobian>(j);

      // Require fields to be evalauted
      fm.requireField<Residual>(*scatter_tag_r);
      fm.requireField<Jacobian>(*scatter_tag_j);
    }

    // ********************************************* 
    // This block is not in the base FiniteElementAssembly example. It
    // is used to demonstrate and unit test a power user feature for
    // building the DAG separate from setting the field sizes and
    // allocating/binding the field memory.
    {
      const auto& tags = fm.getFieldTagsForSizing<Residual>();
      TEUCHOS_ASSERT(tags.size() == 0);
      fm.buildDagForType<Residual>();
      TEUCHOS_ASSERT(tags.size() == static_cast<std::size_t>(num_equations * 4 + 1)); 

      for (auto& t : tags) {
        if (t->dataLayout().identifier() == "DL<CELL,QP>")
          t->nonConstDataLayout().setExtents(workset_size,8);
        else if (t->dataLayout().identifier() == "DL<CELL,QP,DIM>")
          t->nonConstDataLayout().setExtents(workset_size,8,3);
        else if (t->dataLayout().identifier() == "DL<CELL,BASIS>")
          t->nonConstDataLayout().setExtents(workset_size,8);
        else if (t->dataLayout().identifier() == "DL<SCATTER>")
          t->nonConstDataLayout().setExtents(0);
        else
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                                     "ERROR: unknown DataLayout identifier:" <<
                                     t->dataLayout().identifier());
      }
    }
    // *********************************************

    fm.postRegistrationSetup(nullptr);
    fm.writeGraphvizFile("example_fem",".dot",true,true);

    // In debug mode, print the evaluator start stop messages
    fm.printEvaluatorStartStopMessage<Residual>(Teuchos::rcpFromRef(std::cout));

    // Set the team and vector size on the workset for Host DAG
    for (auto& w : worksets) {
      w.team_size_ = p.teamSize();
      w.vector_size_ = p.vectorSize();
    }

    // Kokkos::deep_copy(x,1.0);
    Kokkos::parallel_for(x.extent(0),KOKKOS_LAMBDA (const int& i) {x(i)=static_cast<double>(i);});
    Kokkos::deep_copy(f,0.0);
    RCP<Time> residual_eval_time = TimeMonitor::getNewTimer("Residual Evaluation Time");
    typename PHX::exec_space().fence();
    if (p.doResidual()) {
      TimeMonitor tm_r(*residual_eval_time);
      for (const auto& workset : worksets)
        fm.evaluateFields<Residual>(workset);
      typename PHX::exec_space().fence();
    }

    if (p.printResidual())
      phx_example::printResidual(f,"FEA_MFT: <Residual>",p.printToFile(),"FEA_MFT.Residual.txt");
    
    // Jacobian does both f and J
    Kokkos::deep_copy(f,0.0);
    Kokkos::deep_copy(J.values,0.0);
    RCP<Time> jacobian_eval_time = TimeMonitor::getNewTimer("Jacobian Evaluation Time");
    typename PHX::exec_space().fence();
    if (p.doJacobian()) {
      TimeMonitor tm_r(*jacobian_eval_time);
      for (const auto& workset : worksets)
        fm.evaluateFields<Jacobian>(workset);
      typename PHX::exec_space().fence();
    }

    if (p.printJacobian())
      phx_example::printResidualAndJacobian(f,J,"FEA_MFT: <Jacobian>",p.printToFile(),"FEA_MFT.Jacobian.txt");
    
    // Graph analysis
    if (p.doGraphAnalysis()) {
      double scalability = 0.0;
      double parallelizability = 1.0;
      fm.analyzeGraph<MyTraits::Residual>(scalability,parallelizability);
      std::cout << "Task Scalability       (Residual) = " << scalability << std::endl;
      std::cout << "Task Parallelizability (Residual) = " 
		<< parallelizability << std::endl;
      fm.analyzeGraph<MyTraits::Jacobian>(scalability,parallelizability);
      std::cout << "Task Scalability       (Jacobian) = " << scalability << std::endl;
      std::cout << "Task Parallelizability (Jacobian) = " 
		<< parallelizability << std::endl;
    }
    
  }
  catch (const std::exception& e) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Exception Caught!" << endl;
    std::cout << "Error message is below\n " << e.what() << endl;
    std::cout << "************************************************" << endl;
    success = -1;
  }
  catch (...) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Unknown Exception Caught!" << endl;
    std::cout << "************************************************" << endl;
    success = -1;
  }

  TimeMonitor::summarize();

  return success;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
