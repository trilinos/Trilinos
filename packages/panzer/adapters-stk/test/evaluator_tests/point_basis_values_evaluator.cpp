// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_PointValues_Evaluator.hpp"
#include "Panzer_BasisValues_Evaluator.hpp"
#include "Panzer_DOF.hpp"
#include "Panzer_DOF_PointValues.hpp"
#include "Panzer_Constant.hpp"
#include "Panzer_IntegrationValues2.hpp"
#include "Panzer_BasisValues2.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_WorksetFactory.hpp"

#include "RandomFieldEvaluator.hpp"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "user_app_EquationSetFactory.hpp"

#include <cstdio> // for get char
#include <vector>
#include <string>

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

  Teuchos::RCP<panzer::PureBasis> buildBasis(std::size_t worksetSize,const std::string & basisName);
  void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb, int integration_order);
  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY);
  Teuchos::RCP<panzer::IntegrationRule> buildIR(std::size_t worksetSize,int cubature_degree);

  TEUCHOS_UNIT_TEST(point_values_evaluator, eval)
  {

    const std::size_t workset_size = 4;
    const std::string fieldName_q1 = "U";
    const std::string fieldName_qedge1 = "V";

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(2,2);

    // build input physics block
    Teuchos::RCP<panzer::PureBasis> basis_q1 = buildBasis(workset_size,"Q1");
    panzer::CellData cell_data(workset_size, mesh->getCellTopology("eblock-0_0"));

    const int integration_order = 1;
    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList();
    testInitialization(ipb,integration_order);

    const int default_int_order = integration_order;
    std::string eBlockID = "eblock-0_0";
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();
    Teuchos::RCP<panzer::PhysicsBlock> physicsBlock =
      Teuchos::rcp(new PhysicsBlock(ipb,eBlockID,default_int_order,cell_data,eqset_factory,gd,false));

    Teuchos::RCP<std::vector<panzer::Workset> > work_sets = panzer_stk::buildWorksets(*mesh,physicsBlock->elementBlockID(),
                                                                                            physicsBlock->getWorksetNeeds());
    TEST_EQUALITY(work_sets->size(),1);

    int num_points = 3;
    RCP<const panzer::PointRule> point_rule = rcp(new panzer::PointRule("RandomPoints",num_points, cell_data));
    RCP<const panzer::PointRule> point_rule_basis = rcp(new panzer::PointRule("BasisPoints",basis_q1->cardinality(), cell_data));

    Teuchos::RCP<Kokkos::DynRankView<double,PHX::Device> > userArray
      = Teuchos::rcp(new Kokkos::DynRankView<double,PHX::Device>("userArray",num_points,2));
    Kokkos::DynRankView<double,PHX::Device>&  point_coordinates = *userArray;
    auto point_coordinates_h = Kokkos::create_mirror_view(point_coordinates);
    point_coordinates_h(0,0) =  0.0; point_coordinates_h(0,1) = 0.0; // mid point
    point_coordinates_h(1,0) =  0.5; point_coordinates_h(1,1) = 0.5; // mid point of upper left quadrant
    point_coordinates_h(2,0) = -0.5; point_coordinates_h(2,1) = 0.0; // mid point of line from center to left side
    Kokkos::deep_copy(point_coordinates, point_coordinates_h);

    // setup field manager, add evaluator under test
    /////////////////////////////////////////////////////////////

    PHX::FieldManager<panzer::Traits> fm;

    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator
          = Teuchos::rcp(new panzer::PointValues_Evaluator<panzer::Traits::Residual,panzer::Traits>(point_rule,*userArray));

       TEST_EQUALITY(evaluator->evaluatedFields().size(),6);

       fm.registerEvaluator<panzer::Traits::Residual>(evaluator);
       auto required_fields = evaluator->evaluatedFields();
       for (const auto& f : required_fields)
         fm.requireField<panzer::Traits::Residual>(*f);
    }

    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator
          = Teuchos::rcp(new panzer::PointValues_Evaluator<panzer::Traits::Residual,panzer::Traits>(point_rule_basis,basis_q1));

       TEST_EQUALITY(evaluator->evaluatedFields().size(),6);

       fm.registerEvaluator<panzer::Traits::Residual>(evaluator);
       auto required_fields = evaluator->evaluatedFields();
       for (const auto& f : required_fields)
         fm.requireField<panzer::Traits::Residual>(*f);
    }

    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator
          = Teuchos::rcp(new panzer::PointValues_Evaluator<panzer::Traits::Jacobian,panzer::Traits>(point_rule,*userArray));

       TEST_EQUALITY(evaluator->evaluatedFields().size(),6);

       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
       fm.requireField<panzer::Traits::Jacobian>(*evaluator->evaluatedFields()[0]);
    }

    panzer::Traits::SD sd;
    sd.worksets_ = work_sets;
    // run tests
    /////////////////////////////////////////////////////////////

    panzer::Workset & workset = (*work_sets)[0];
    workset.alpha = 0.0;
    workset.beta = 0.0;
    workset.time = 0.0;
    workset.evaluate_transient_terms = false;

    std::vector<PHX::index_size_type> derivative_dimensions;
    derivative_dimensions.push_back(8);
    fm.setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);

    fm.postRegistrationSetup(sd);
    fm.print(out);

    fm.evaluateFields<panzer::Traits::Residual>(workset);
    fm.evaluateFields<panzer::Traits::Jacobian>(workset);

    PHX::MDField<panzer::Traits::Residual::ScalarT>
       point_coords(point_rule->getName()+"_"+"point_coords",point_rule->dl_vector);
    fm.getFieldData<panzer::Traits::Residual>(point_coords);

    PHX::MDField<panzer::Traits::Residual::ScalarT>
       point_coords_basis(point_rule_basis->getName()+"_"+"point_coords",point_rule_basis->dl_vector);
    fm.getFieldData<panzer::Traits::Residual>(point_coords_basis);

    PHX::MDField<panzer::Traits::Residual::ScalarT,Cell,IP>
       point_coords_jac_det(point_rule_basis->getName()+"_"+"jac_det",point_rule_basis->dl_scalar);
    fm.getFieldData<panzer::Traits::Residual,panzer::Traits::Residual::ScalarT,Cell,IP>(point_coords_jac_det);

    PHX::MDField<panzer::Traits::Residual::ScalarT,Cell,IP,Dim,Dim>
       point_coords_jac(point_rule_basis->getName()+"_"+"jac",point_rule_basis->dl_tensor);
    fm.getFieldData<panzer::Traits::Residual,panzer::Traits::Residual::ScalarT,Cell,IP,Dim,Dim>(point_coords_jac);

    PHX::MDField<panzer::Traits::Residual::ScalarT,Cell,IP,Dim,Dim>
       point_coords_jac_inv(point_rule_basis->getName()+"_"+"jac_inv",point_rule_basis->dl_tensor);
    fm.getFieldData<panzer::Traits::Residual,panzer::Traits::Residual::ScalarT,Cell,IP,Dim,Dim>(point_coords_jac_inv);

    auto point_coords_h = Kokkos::create_mirror_view(point_coords.get_view());
    Kokkos::deep_copy(point_coords_h, point_coords.get_view());
    auto point_coords_basis_h = Kokkos::create_mirror_view(point_coords_basis.get_view());
    Kokkos::deep_copy(point_coords_basis_h, point_coords_basis.get_view());
    auto point_coords_jac_det_h = Kokkos::create_mirror_view(point_coords_jac_det.get_view());
    Kokkos::deep_copy(point_coords_jac_det_h, point_coords_jac_det.get_view());
    auto point_coords_jac_h = Kokkos::create_mirror_view(point_coords_jac.get_view());
    Kokkos::deep_copy(point_coords_jac_h, point_coords_jac.get_view());
    auto point_coords_jac_inv_h = Kokkos::create_mirror_view(point_coords_jac_inv.get_view());
    Kokkos::deep_copy(point_coords_jac_inv_h, point_coords_jac_inv.get_view());
    auto cell_node_coordinates_h = Kokkos::create_mirror_view(workset.cell_node_coordinates.get_view());
    Kokkos::deep_copy(cell_node_coordinates_h, workset.cell_node_coordinates.get_view());
    auto basis_coordinates_ref_h = Kokkos::create_mirror_view(workset.bases[1]->basis_coordinates_ref.get_view());
    Kokkos::deep_copy(basis_coordinates_ref_h, workset.bases[1]->basis_coordinates_ref.get_view());
    for(int c=0;c<basis_q1->numCells();c++) {
       double dx = 0.5;
       double dy = 0.5;
       for(int p=0;p<num_points;p++) {
          double x = dx*(point_coordinates_h(p,0)+1.0)/2.0 + cell_node_coordinates_h(c,0,0);
          double y = dy*(point_coordinates_h(p,1)+1.0)/2.0 + cell_node_coordinates_h(c,0,1);
          TEST_FLOATING_EQUALITY(point_coords_h(c,p,0),x,1e-10);
          TEST_FLOATING_EQUALITY(point_coords_h(c,p,1),y,1e-10);
       }

       for(int p=0;p<basis_q1->cardinality();p++) {
          double x = dx*(basis_coordinates_ref_h(p,0)+1.0)/2.0 + cell_node_coordinates_h(c,0,0);
          double y = dy*(basis_coordinates_ref_h(p,1)+1.0)/2.0 + cell_node_coordinates_h(c,0,1);
          TEST_FLOATING_EQUALITY(point_coords_basis_h(c,p,0),x,1e-10);
          TEST_FLOATING_EQUALITY(point_coords_basis_h(c,p,1),y,1e-10);
       }

       for(int p=0;p<num_points;p++)
          TEST_FLOATING_EQUALITY(point_coords_jac_det_h(c,p),dx*dy/4.0,1e-10);

       for(int p=0;p<num_points;p++) {
          TEST_FLOATING_EQUALITY(point_coords_jac_h(c,p,0,0),dx/2.0,1e-10);
          TEST_FLOATING_EQUALITY(point_coords_jac_h(c,p,0,1),0.0,1e-10);
          TEST_FLOATING_EQUALITY(point_coords_jac_h(c,p,1,0),0.0,1e-10);
          TEST_FLOATING_EQUALITY(point_coords_jac_h(c,p,1,1),dy/2.0,1e-10);
       }

       for(int p=0;p<num_points;p++) {
          TEST_FLOATING_EQUALITY(point_coords_jac_inv_h(c,p,0,0),1.0/(dx/2.0),1e-10);
          TEST_FLOATING_EQUALITY(point_coords_jac_inv_h(c,p,0,1),0.0,1e-10);
          TEST_FLOATING_EQUALITY(point_coords_jac_inv_h(c,p,1,0),0.0,1e-10);
          TEST_FLOATING_EQUALITY(point_coords_jac_inv_h(c,p,1,1),1.0/(dy/2.0),1e-10);
       }
    }
  }

  TEUCHOS_UNIT_TEST(basis_values_evaluator, eval)
  {

    const std::size_t workset_size = 4;
    const std::string fieldName_q1 = "U";
    const std::string fieldName_qedge1 = "V";

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(2,2);

    // build input physics block
    Teuchos::RCP<panzer::PureBasis> basis_q1 = buildBasis(workset_size,"Q1");
    panzer::CellData cell_data(basis_q1->numCells(), basis_q1->getCellTopology());

    int integration_order = 4;
    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList();
    testInitialization(ipb,integration_order);

    const int default_int_order = integration_order;
    std::string eBlockID = "eblock-0_0";
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    panzer::CellData cellData(workset_size,mesh->getCellTopology("eblock-0_0"));
    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();
    Teuchos::RCP<panzer::PhysicsBlock> physicsBlock =
      Teuchos::rcp(new PhysicsBlock(ipb,eBlockID,default_int_order,cellData,eqset_factory,gd,false));

    Teuchos::RCP<std::vector<panzer::Workset> > work_sets = panzer_stk::buildWorksets(*mesh,physicsBlock->elementBlockID(),
                                                                                            physicsBlock->getWorksetNeeds());
    panzer::Workset & workset = (*work_sets)[0];
    TEST_EQUALITY(work_sets->size(),1);

    Teuchos::RCP<panzer::IntegrationRule> point_rule = buildIR(workset_size,integration_order);
    panzer::IntegrationValues2<double> int_values("",true);
    int_values.setupArrays(point_rule);
    int_values.evaluateValues(workset.cell_node_coordinates);

    // Teuchos::RCP<Kokkos::DynRankView<double,PHX::Device> > userArray = Teuchos::rcpFromRef(int_values.cub_points);
    auto userArray = int_values.cub_points;

    Teuchos::RCP<panzer::BasisIRLayout> layout = Teuchos::rcp(new panzer::BasisIRLayout(basis_q1,*point_rule));

    // setup field manager, add evaluator under test
    /////////////////////////////////////////////////////////////

    PHX::FieldManager<panzer::Traits> fm;

    Teuchos::RCP<const std::vector<Teuchos::RCP<PHX::FieldTag > > > evalJacFields;
    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator
          = Teuchos::rcp(new panzer::PointValues_Evaluator<panzer::Traits::Jacobian,panzer::Traits>(point_rule,userArray));
       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
    }
    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator
          = Teuchos::rcp(new panzer::BasisValues_Evaluator<panzer::Traits::Jacobian,panzer::Traits>(point_rule,basis_q1));

       TEST_EQUALITY(evaluator->evaluatedFields().size(),4);
       evalJacFields = Teuchos::rcpFromRef(evaluator->evaluatedFields());

       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
       fm.requireField<panzer::Traits::Jacobian>(*evaluator->evaluatedFields()[0]);
    }

    std::vector<PHX::index_size_type> derivative_dimensions;
    derivative_dimensions.push_back(8);
    fm.setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);

    panzer::Traits::SD sd;
    fm.postRegistrationSetup(sd);
    fm.print(out);

    // run tests
    /////////////////////////////////////////////////////////////

    workset.alpha = 0.0;
    workset.beta = 0.0;
    workset.time = 0.0;
    workset.evaluate_transient_terms = false;

    fm.evaluateFields<panzer::Traits::Jacobian>(workset);

    PHX::MDField<panzer::Traits::Residual::ScalarT,panzer::Cell,panzer::BASIS,panzer::IP>
       basis(basis_q1->name()+"_"+point_rule->getName()+"_"+"basis",layout->basis);
    fm.getFieldData<panzer::Traits::Jacobian>(basis);
    out << basis << std::endl;

    WorksetDetailsAccessor wda;
    std::size_t basisIndex = panzer::getBasisIndex(layout->name(), workset, wda);
    Teuchos::RCP<panzer::BasisValues2<double> > bases = workset.bases[basisIndex];
    TEST_ASSERT(bases!=Teuchos::null);
    // TEST_EQUALITY(bases->basis.size(),basis.size());
    auto basis_h = Kokkos::create_mirror_view(basis.get_view());
    Kokkos::deep_copy(basis_h, basis.get_view());
    auto basis_scalar_h = Kokkos::create_mirror_view(bases->basis_scalar.get_static_view());
    Kokkos::deep_copy(basis_scalar_h, bases->basis_scalar.get_static_view());
    for(int i=0;i<4;i++) {
      for(int j=0;j<4;j++) {
        for(unsigned int k=0;k<bases->basis_scalar.extent(2);k++) {
          TEST_FLOATING_EQUALITY(basis_scalar_h(i,j,k),basis_h(i,j,k),1e-10);
        }
      }
    }
  }

  TEUCHOS_UNIT_TEST(basis_values_evaluator, eval_vector)
  {
    typedef Intrepid2::Basis<PHX::Device::execution_space,double,double> IntrepidBasis;

    const std::size_t workset_size = 4;
    const std::string fieldName = "U";

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(2,2);

    // build input physics block
    Teuchos::RCP<panzer::PureBasis> basis_edge = buildBasis(workset_size,"HCurl");
    panzer::CellData cell_data(basis_edge->numCells(), basis_edge->getCellTopology());

    int integration_order = 4;
    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList();
    testInitialization(ipb,integration_order);

    std::string eBlockID = "eblock-0_0";
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    panzer::CellData cellData(workset_size,mesh->getCellTopology("eblock-0_0"));
    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();
    Teuchos::RCP<panzer::PhysicsBlock> physicsBlock =
      Teuchos::rcp(new PhysicsBlock(ipb,eBlockID,integration_order,cellData,eqset_factory,gd,false));

    std::map<std::string,WorksetNeeds> needs;
    needs[physicsBlock->elementBlockID()] = physicsBlock->getWorksetNeeds();

    // build DOF Manager (with a single HDiv basis)
    /////////////////////////////////////////////////////////////

    // build the connection manager
    const RCP<panzer::ConnManager>
      conn_manager = rcp(new panzer_stk::STKConnManager(mesh));

    RCP<panzer::DOFManager> dof_manager
        = rcp(new panzer::DOFManager(conn_manager,MPI_COMM_WORLD));

    // build an intrepid basis and a related field pattern for seeding the DOFManager
    {
       RCP<IntrepidBasis> hgrad_intrepid_basis;
       hgrad_intrepid_basis
           = panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>("HGrad",1,
                                                                                      *mesh->getCellTopology(physicsBlock->elementBlockID()));
      RCP<panzer::Intrepid2FieldPattern> hgrad_field_pattern = rcp(new panzer::Intrepid2FieldPattern(hgrad_intrepid_basis));

      dof_manager->addField(physicsBlock->elementBlockID(), fieldName, hgrad_field_pattern);
    }
    dof_manager->buildGlobalUnknowns();

    /////////////////////////////////////////////////////////////

    RCP<panzer_stk::WorksetFactory> wkstFactory
       = rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = rcp(new panzer::WorksetContainer(wkstFactory,needs));
    wkstContainer->setGlobalIndexer(dof_manager);
    wkstContainer->setWorksetSize(workset_size);

    Teuchos::RCP<std::vector<panzer::Workset> > work_sets = wkstContainer->getWorksets(blockDescriptor(physicsBlock->elementBlockID()));
    panzer::Workset & workset = (*work_sets)[0];
    TEST_EQUALITY(work_sets->size(),1);

    Teuchos::RCP<panzer::IntegrationRule> point_rule = buildIR(workset_size,integration_order);
    panzer::IntegrationValues2<double> int_values("",true);
    int_values.setupArrays(point_rule);
    int_values.evaluateValues(workset.cell_node_coordinates);

    // Teuchos::RCP<Kokkos::DynRankView<double,PHX::Device> > userArray = Teuchos::rcpFromRef(int_values.cub_points);
    auto userArray = int_values.cub_points;

    Teuchos::RCP<panzer::BasisIRLayout> layout = Teuchos::rcp(new panzer::BasisIRLayout(basis_edge,*point_rule));

    // setup field manager, add evaluator under test
    /////////////////////////////////////////////////////////////

    PHX::FieldManager<panzer::Traits> fm;

    Teuchos::RCP<const std::vector<Teuchos::RCP<PHX::FieldTag > > > evalJacFields;
    {
       Teuchos::ParameterList input;
       input.set("Name",  "HCurl:1 Orientation");
       input.set("Value", 1.0);
       input.set("Data Layout", basis_edge->functional);
       Teuchos::RCP< PHX::Evaluator<panzer::Traits> > evaluator
          = rcp(new panzer::Constant<panzer::Traits::Jacobian,panzer::Traits>(input));
       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
    }
    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator
          = Teuchos::rcp(new panzer::PointValues_Evaluator<panzer::Traits::Jacobian,panzer::Traits>(point_rule,userArray));
       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
    }
    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator
          = Teuchos::rcp(new panzer::BasisValues_Evaluator<panzer::Traits::Jacobian,panzer::Traits>(point_rule,basis_edge));

       TEST_EQUALITY(evaluator->evaluatedFields().size(),4);
       evalJacFields = Teuchos::rcpFromRef(evaluator->evaluatedFields());

       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
       fm.requireField<panzer::Traits::Jacobian>(*evaluator->evaluatedFields()[0]);
       out << "REQUIRED FIELD = \"";
       evaluator->evaluatedFields()[0]->print(out);
       out << "\"" << std::endl;
    }

    std::vector<PHX::index_size_type> derivative_dimensions;
    derivative_dimensions.push_back(4);
    fm.setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);

    panzer::Traits::SD sd;
    sd.orientations_ = wkstContainer->getOrientations();
    fm.postRegistrationSetup(sd);
    fm.print(out);

    // run tests
    /////////////////////////////////////////////////////////////

    workset.alpha = 0.0;
    workset.beta = 0.0;
    workset.time = 0.0;
    workset.evaluate_transient_terms = false;

    fm.evaluateFields<panzer::Traits::Jacobian>(workset);

    // NOTE: basis values are always double, not Fad.
    PHX::MDField<panzer::Traits::Residual::ScalarT,Cell,BASIS,IP,Dim>
       basis(basis_edge->name()+"_"+point_rule->getName()+"_basis",layout->basis_grad);
    PHX::MDField<panzer::Traits::Residual::ScalarT,Cell,BASIS,IP>
       curl_basis(basis_edge->name()+"_"+point_rule->getName()+"_curl_basis",layout->basis);
    PHX::MDField<panzer::Traits::Residual::ScalarT,Cell,IP,Dim,Dim>
       jac_inv("CubaturePoints (Degree=4,volume)_jac_inv",point_rule->dl_tensor);

    fm.getFieldData<panzer::Traits::Jacobian,panzer::Traits::Residual::ScalarT,Cell,BASIS,IP,Dim>(basis);
    fm.getFieldData<panzer::Traits::Jacobian,panzer::Traits::Residual::ScalarT,Cell,BASIS,IP>(curl_basis);
    fm.getFieldData<panzer::Traits::Jacobian,panzer::Traits::Residual::ScalarT,Cell,IP,Dim,Dim>(jac_inv);

    WorksetDetailsAccessor wda;
    std::size_t basisIndex = panzer::getBasisIndex(layout->name(), workset, wda);
    Teuchos::RCP<panzer::BasisValues2<double> > bases = workset.bases[basisIndex];
    TEST_ASSERT(bases!=Teuchos::null);
    TEST_EQUALITY(bases->basis_vector.size(),basis.size());
    TEST_EQUALITY(bases->curl_basis_scalar.size(),curl_basis.size());

    auto basis_h = Kokkos::create_mirror_view(basis.get_view());
    Kokkos::deep_copy(basis_h, basis.get_view());
    auto basis_vector_h = Kokkos::create_mirror_view(bases->basis_vector.get_static_view());
    Kokkos::deep_copy(basis_vector_h,bases->basis_vector.get_static_view());
    for(int i=0;i<4;i++) {
      for(int j=0;j<4;j++) {
        for(unsigned int k=0;k<bases->curl_basis_scalar.extent(2);k++) {
          for(unsigned int d=0;d<bases->basis_vector.extent(3);d++) {
            TEST_FLOATING_EQUALITY(basis_vector_h(i,j,k,d),basis_h(i,j,k,d),1e-10);
          }
        }
      }
    }

    auto curl_basis_h = Kokkos::create_mirror_view(curl_basis.get_view());
    Kokkos::deep_copy(curl_basis_h, curl_basis.get_view());
    auto curl_basis_scalar_h = Kokkos::create_mirror_view(bases->curl_basis_scalar.get_static_view());
    Kokkos::deep_copy(curl_basis_scalar_h,bases->curl_basis_scalar.get_static_view());
    for(int i=0;i<4;i++) {
      for(int j=0;j<4;j++) {
        for(unsigned int k=0;k<bases->curl_basis_scalar.extent(2);k++) {
          TEST_FLOATING_EQUALITY(curl_basis_scalar_h(i,j,k),curl_basis_h(i,j,k),1e-10);
        }
      }
    }
  }

  TEUCHOS_UNIT_TEST(dof_point_values_evaluator, eval)
  {
    typedef Intrepid2::Basis<PHX::Device::execution_space,double,double> IntrepidBasis;

    const std::size_t workset_size = 4;
    const std::string fieldName_q1 = "U";
    const std::string fieldName_qedge1 = "V";

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(2,2);

    // build input physics block
    Teuchos::RCP<panzer::PureBasis> basis_q1 = buildBasis(workset_size,"Q1");
    panzer::CellData cell_data(basis_q1->numCells(), basis_q1->getCellTopology());

    int integration_order = 4;
    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList();
    testInitialization(ipb,integration_order);

    const int default_int_order = integration_order;
    std::string eBlockID = "eblock-0_0";
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    panzer::CellData cellData(workset_size,mesh->getCellTopology("eblock-0_0"));
    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();
    Teuchos::RCP<panzer::PhysicsBlock> physicsBlock =
      Teuchos::rcp(new PhysicsBlock(ipb,eBlockID,default_int_order,cellData,eqset_factory,gd,false));

    std::map<std::string,WorksetNeeds> needs;
    needs[physicsBlock->elementBlockID()] = physicsBlock->getWorksetNeeds();

    // build DOF Manager (with a single HDiv basis)
    /////////////////////////////////////////////////////////////

    // build the connection manager
    const RCP<panzer::ConnManager>
      conn_manager = rcp(new panzer_stk::STKConnManager(mesh));

    RCP<panzer::DOFManager> dof_manager
        = rcp(new panzer::DOFManager(conn_manager,MPI_COMM_WORLD));

    // build an intrepid basis and a related field pattern for seeding the DOFManager
    {
       RCP<IntrepidBasis> hgrad_intrepid_basis, hcurl_intrepid_basis;
       hgrad_intrepid_basis
           = panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>("HGrad",1,
                                                                                      *mesh->getCellTopology(physicsBlock->elementBlockID()));
       hcurl_intrepid_basis
           = panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>("HCurl",1,
                                                                                      *mesh->getCellTopology(physicsBlock->elementBlockID()));
      RCP<panzer::Intrepid2FieldPattern> hgrad_field_pattern = rcp(new panzer::Intrepid2FieldPattern(hgrad_intrepid_basis));
      RCP<panzer::Intrepid2FieldPattern> hcurl_field_pattern = rcp(new panzer::Intrepid2FieldPattern(hcurl_intrepid_basis));

      dof_manager->addField(physicsBlock->elementBlockID(), fieldName_q1, hgrad_field_pattern);
      dof_manager->addField(physicsBlock->elementBlockID(), fieldName_qedge1, hcurl_field_pattern);
    }
    dof_manager->buildGlobalUnknowns();

    /////////////////////////////////////////////////////////////

    RCP<panzer_stk::WorksetFactory> wkstFactory
       = rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = rcp(new panzer::WorksetContainer(wkstFactory,needs));
    wkstContainer->setGlobalIndexer(dof_manager);
    wkstContainer->setWorksetSize(workset_size);

    // Teuchos::RCP<std::vector<panzer::Workset> > work_sets = panzer_stk::buildWorksets(*mesh,physicsBlock->elementBlockID(),
    //                                                                                         physicsBlock->getWorksetNeeds());
    Teuchos::RCP<std::vector<panzer::Workset> > work_sets = wkstContainer->getWorksets(blockDescriptor(physicsBlock->elementBlockID()));
    panzer::Workset & workset = (*work_sets)[0];

    TEST_EQUALITY(work_sets->size(),1);

    Teuchos::RCP<panzer::IntegrationRule> point_rule = buildIR(workset_size,integration_order);
    panzer::IntegrationValues2<double> int_values("",true);
    int_values.setupArrays(point_rule);
    int_values.evaluateValues(workset.cell_node_coordinates);

    // Teuchos::RCP<Kokkos::DynRankView<double,PHX::Device> > userArray = Teuchos::rcpFromRef(int_values.cub_points);
    auto userArray = int_values.cub_points;

    Teuchos::RCP<panzer::BasisIRLayout> layout = Teuchos::rcp(new panzer::BasisIRLayout(basis_q1,*point_rule));

    // setup field manager, add evaluator under test
    /////////////////////////////////////////////////////////////

    PHX::FieldManager<panzer::Traits> fm;
    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator
          = Teuchos::rcp(new RandomFieldEvaluator<panzer::Traits::Jacobian,panzer::Traits>("TEMPERATURE",basis_q1->functional));
       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
    }
    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator
          = Teuchos::rcp(new panzer::PointValues_Evaluator<panzer::Traits::Jacobian,panzer::Traits>(point_rule,userArray));
       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
    }
    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator
          = Teuchos::rcp(new panzer::BasisValues_Evaluator<panzer::Traits::Jacobian,panzer::Traits>(point_rule,basis_q1));
       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
    }
    {
       Teuchos::ParameterList p;
       p.set("Name","TEMPERATURE");
       p.set("Basis",layout);
       p.set("IR",point_rule);
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator
          = Teuchos::rcp(new panzer::DOF<panzer::Traits::Jacobian,panzer::Traits>(p));
       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
       fm.requireField<panzer::Traits::Jacobian>(*evaluator->evaluatedFields()[0]); // require DOF
    }
    {
       Teuchos::ParameterList p;
       p.set("Name","TEMPERATURE");
       p.set("Basis",layout->getBasis().getConst());
       p.set("Point Rule",Teuchos::rcp_static_cast<panzer::PointRule>(point_rule).getConst());
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator
          = Teuchos::rcp(new panzer::DOF_PointValues<panzer::Traits::Jacobian,panzer::Traits>(p));

       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
       fm.requireField<panzer::Traits::Jacobian>(*evaluator->evaluatedFields()[0]); // require DOF_PointValues
    }

    std::vector<PHX::index_size_type> derivative_dimensions;
    derivative_dimensions.push_back(8);
    fm.setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);

    panzer::Traits::SD sd;
    sd.orientations_ = wkstContainer->getOrientations();
    sd.worksets_ = work_sets;
    fm.postRegistrationSetup(sd);
    fm.print(out);

    // run tests
    /////////////////////////////////////////////////////////////

    workset.alpha = 0.0;
    workset.beta = 0.0;
    workset.time = 0.0;
    workset.evaluate_transient_terms = false;

    panzer::Traits::PED ped;
    fm.preEvaluate<panzer::Traits::Jacobian>(ped);
    fm.evaluateFields<panzer::Traits::Jacobian>(workset);
    fm.postEvaluate<panzer::Traits::Jacobian>(NULL);

    PHX::MDField<panzer::Traits::Jacobian::ScalarT> ref_field("TEMPERATURE",point_rule->dl_scalar);
    PHX::MDField<panzer::Traits::Jacobian::ScalarT> fut_field("TEMPERATURE_"+point_rule->getName(),point_rule->dl_scalar); // "field under test"
    fm.getFieldData<panzer::Traits::Jacobian>(ref_field);
    fm.getFieldData<panzer::Traits::Jacobian>(fut_field);

    TEST_EQUALITY(ref_field.size(),fut_field.size());
    // for(int i=0;i<ref_field.size();i++) {
    //   TEST_FLOATING_EQUALITY(fut_field[i].val(),ref_field[i].val(),1e-10);
    // }
    auto fut_field_h = Kokkos::create_mirror_view(fut_field.get_view());
    Kokkos::deep_copy(fut_field_h, fut_field.get_view());
    auto ref_field_h = Kokkos::create_mirror_view(ref_field.get_view());
    Kokkos::deep_copy(ref_field_h, ref_field.get_view());
    for(int i=0;i<4;i++) {
      for(int j=0;j<4;j++) {
        TEST_FLOATING_EQUALITY(fut_field_h(i,j).val(),ref_field_h(i,j).val(),1e-10);
      }
    }
  }

  Teuchos::RCP<panzer::PureBasis> buildBasis(std::size_t worksetSize,const std::string & basisName)
  {
     Teuchos::RCP<shards::CellTopology> topo =
        Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

     panzer::CellData cellData(worksetSize,topo);
     return Teuchos::rcp(new panzer::PureBasis(basisName,1,cellData));
  }

  Teuchos::RCP<panzer::IntegrationRule> buildIR(std::size_t workset_size,int cubature_degree)
  {
     Teuchos::RCP<shards::CellTopology> topo =
        Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

     const panzer::CellData cell_data(workset_size, topo);

     return Teuchos::rcp(new panzer::IntegrationRule(cubature_degree, cell_data));
  }

  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY)
  {
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",1);
    pl->set("Y Blocks",1);
    pl->set("X Elements",elemX);
    pl->set("Y Elements",elemY);

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

    return mesh;
  }

  void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb, int integration_order)
  {
    // Physics block
    ipb->setName("test physics");
    {
      Teuchos::ParameterList& p = ipb->sublist("a");
      p.set("Type","Energy");
      p.set("Prefix","");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",integration_order);
    }
    {
      Teuchos::ParameterList& p = ipb->sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","solid");
      p.set("Basis Type","HCurl");
      p.set("Basis Order",1);
      p.set("Integration Order",integration_order);
    }

  }

}
