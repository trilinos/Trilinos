#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Panzer_STK_ProjectField.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Phalanx_FieldManager.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

// TODO BWR MOVE W FUNCS
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_CubeTetMeshFactory.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SquareTriMeshFactory.hpp"
#include "Panzer_STK_LineMeshFactory.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_LagrangianInterpolation.hpp"

typedef Kokkos::DynRankView<double,PHX::Device> DynRankView;

// for making explicit instantiated tests easier
#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(project_field,value,TYPE)

Teuchos::RCP<panzer_stk::STK_Interface> createInlineMesh(Teuchos::RCP<Teuchos::ParameterList> pl);
Teuchos::RCP<panzer::WorksetContainer> createWorksetsForFields(Teuchos::RCP<panzer_stk::STK_Interface> mesh, 
                                                               Teuchos::RCP<panzer::DOFManager> dof_manager, 
                                                               std::map<std::string,panzer::BasisDescriptor> fmap);
DynRankView getCoeffsForAnalyticFunction(Teuchos::RCP<panzer::PureBasis> basis, DynRankView cellNodes, 
                                         Teuchos::RCP<const std::vector<Intrepid2::Orientation> > orts /*, Functor ?? */ );

struct Fun {
double
KOKKOS_INLINE_FUNCTION
operator()(const double& x, const double& y, const double& z) const {
  return x*x + y*y + z*z;
}
};

template<typename DynRankViewType>
class EvalSolFunctor {
public:
  DynRankViewType funAtPoints;
  DynRankViewType points;
  Fun fun;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int elem) const
  {
    for(int i=0;i<static_cast<int>(points.extent(1));i++) {
      auto x = points(elem,i,0), y = points(elem,i,1), z= points(elem,i,2);
      funAtPoints(elem,i) = fun(x,y,z);
    }
  }
};

namespace panzer_stk {

typedef Kokkos::DynRankView<double,PHX::Device> FieldArray;

//**********************************************************************
template<typename EvalT, typename Traits>
class GetCoeffsEvaluator
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    GetCoeffsEvaluator(std::string name, Teuchos::RCP<panzer::PureBasis> basis,
                       DynRankView cellNodes, 
                       Teuchos::RCP<const std::vector<Intrepid2::Orientation> > orts /*, Functor ??*/);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;
    Teuchos::RCP<panzer::PureBasis> basis;
    PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS> coeffs;
    Teuchos::RCP<const std::vector<Intrepid2::Orientation> > orts;
    DynRankView cellNodes;
  
}; // end of class

template<typename EvalT, typename Traits>
GetCoeffsEvaluator<EvalT, Traits>::
GetCoeffsEvaluator(std::string name, Teuchos::RCP<panzer::PureBasis> basis, DynRankView cellNodes, 
                   Teuchos::RCP<const std::vector<Intrepid2::Orientation> > orts /*, Functor ?? */) 
  : basis(basis), orts(orts), cellNodes(cellNodes)
{

  // grab information from quadrature rule
  coeffs = PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS>(name, basis->functional);

  this->addEvaluatedField(coeffs);

  std::string n = "GetCoeffsEvaluator: " + name;
  this->setName(n);

}

template<typename EvalT, typename Traits>
void
GetCoeffsEvaluator<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  /* workset */)
{
  typedef Intrepid2::Experimental::LagrangianInterpolation<PHX::Device> li;
  typedef Intrepid2::CellTools<PHX::Device> ct;

  auto numOwnedElems = cellNodes.extent(0);
  auto dim = basis->dimension();
  auto basisCardinality = basis->cardinality();

  auto dofCoeffs = DynRankView("dofCoeffs",numOwnedElems,basisCardinality);
  auto dofCoords = DynRankView("dofCoords",numOwnedElems,basisCardinality,dim);

  // First, need to copy orientations to device
  auto orts_dev = Kokkos::DynRankView<Intrepid2::Orientation,PHX::Device>("orts_dev",orts->size());
  auto orts_host = Kokkos::create_mirror_view(orts_dev);
  for (size_t i=0; i < orts_host.extent(0); ++i)
    orts_host(i) = orts->at(i);
  Kokkos::deep_copy(orts_dev,orts_host);

  auto it2basis = basis->template getIntrepid2Basis<PHX::exec_space,double,double>();

  // get Dof coordinates and coefficients (reference frame)
  li::getDofCoordsAndCoeffs(dofCoords,dofCoeffs,it2basis.get(),orts_dev);

  auto funAtDofPoints = DynRankView("funAtDofPoints",numOwnedElems, basisCardinality);
    {
      auto physDofPoints = DynRankView("physDofPoints", numOwnedElems, basisCardinality, dim);
      // TODO BWR update this?
      ct::mapToPhysicalFrame(physDofPoints,dofCoords,cellNodes,it2basis->getBaseCellTopology());

      EvalSolFunctor<DynRankView> functor;
      functor.funAtPoints = funAtDofPoints;
      functor.points = physDofPoints;
      Kokkos::parallel_for("loop for evaluating the function at DoF points", numOwnedElems,functor);
      Kokkos::fence(); //make sure that funAtDofPoints has been evaluated
    }

    li::getBasisCoeffs(coeffs.get_view(), funAtDofPoints, dofCoeffs);

  return;
}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(project_field,value,EvalType)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  // build a mesh and worksets
  //////////////////////////////////////////////////////////
  Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
  pl->set("Mesh Dimension",3);
  pl->set("Type","Hex");
  pl->sublist("Mesh Factory Parameter List").set("X Blocks",1);
  pl->sublist("Mesh Factory Parameter List").set("Y Blocks",1);
  pl->sublist("Mesh Factory Parameter List").set("Z Blocks",1);
  pl->sublist("Mesh Factory Parameter List").set("X Elements",2);
  pl->sublist("Mesh Factory Parameter List").set("Y Elements",2);
  pl->sublist("Mesh Factory Parameter List").set("Z Elements",2);

  auto mesh = createInlineMesh(pl);
  Teuchos::RCP<panzer::DOFManager> dof_manager;
  std::map<std::string,panzer::BasisDescriptor> fmap;
  // TODO CHANGE THIS >>>
  fmap["Pressure"] = panzer::BasisDescriptor(2,"HCurl");

  // set up worksets and orientations
  auto workset_container = createWorksetsForFields(mesh,dof_manager,fmap);
  auto orientations = workset_container->getOrientations();
  std::cout << orientations->size() << std::endl;
  auto worksets = workset_container->getWorksets(panzer::WorksetDescriptor("eblock-0_0_0", panzer::WorksetSizeType::ALL_ELEMENTS, true, true));

  auto cellNodes = (*worksets)[0].getCellVertices(); // TODO BWR CHANGE
  auto numCells = cellNodes.extent(0);

  // we will project from a second order basis to a first
  auto topo = mesh->getCellTopology("eblock-0_0_0");
  Teuchos::RCP<panzer::PureBasis> sourceBasis = Teuchos::rcp(new panzer::PureBasis("HGrad",2,numCells,topo));
  Teuchos::RCP<panzer::PureBasis> targetBasis = Teuchos::rcp(new panzer::PureBasis("HGrad",1,numCells,topo));

  Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm
   = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>);

  // need an evaluator to get coeffs of analytic function
  ///////////////////////////////////////////////////
  {
     // fill the basis with some dummy values
     RCP<panzer_stk::GetCoeffsEvaluator<EvalType,panzer::Traits> > e =
       rcp(new panzer_stk::GetCoeffsEvaluator<EvalType,panzer::Traits>(
        "Pressure",sourceBasis,cellNodes.get_view(),orientations)
       );

     fm->registerEvaluator<EvalType>(e);
  }

    // add evaluator under test
    ///////////////////////////////////////////////////

  {
    RCP<panzer_stk::ProjectField<EvalType,panzer::Traits> > e =
      rcp(new panzer_stk::ProjectField<EvalType,panzer::Traits>("Pressure",sourceBasis,targetBasis));

    fm->registerEvaluator<EvalType>(e);

    Teuchos::RCP<PHX::FieldTag> ft = e->evaluatedFields()[0];
    fm->requireField<EvalType>(*ft);
  }

  panzer::Traits::SD setupData;

  setupData.orientations_ = orientations;

  setupData.worksets_ = worksets;
  fm->postRegistrationSetup(setupData);

  panzer::Traits::PED preEvalData;
  fm->preEvaluate<EvalType>(preEvalData);
  // TODO BWR loop thru worksets?
  fm->evaluateFields<EvalType>((*worksets)[0]);
  fm->postEvaluate<EvalType>(0);

  typedef typename EvalType::ScalarT ScalarT;

  typename PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS> s("Pressure",sourceBasis->functional);
  typename PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS> t("Pressure",targetBasis->functional);

  fm->getFieldData<EvalType>(s);
  fm->getFieldData<EvalType>(t);

  auto s_h = Kokkos::create_mirror_view(s.get_view());
  Kokkos::deep_copy(s_h, s.get_view());
  auto t_h = Kokkos::create_mirror_view(t.get_view());
  Kokkos::deep_copy(t_h, t.get_view());

  for (size_t ncell=0;ncell<numCells;++ncell){
    for (size_t idx_dof=0;idx_dof<8;++idx_dof){
      TEST_FLOATING_EQUALITY(s_h(ncell,idx_dof),t_h(ncell,idx_dof),1e-14);
      std::cout << s_h(ncell,idx_dof) << " " << t_h(ncell,idx_dof) << std::endl;
    }
  }

}

// move and wrap inside some test utils namespace?
Teuchos::RCP<panzer_stk::STK_Interface> createInlineMesh(Teuchos::RCP<Teuchos::ParameterList> pl)
{
  std::cout << *pl << std::endl;

  auto dim = pl->get<int>("Mesh Dimension");
  auto type = pl->get<std::string>("Type");

  Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory;

  if (type == "Quad" && dim == 2) {
    mesh_factory = Teuchos::rcp(new panzer_stk::SquareQuadMeshFactory);
  } else if (type == "Tri" && dim == 2) {
    mesh_factory = Teuchos::rcp(new panzer_stk::SquareTriMeshFactory);
  } else if (type == "Tet" && dim == 3) {
    mesh_factory = Teuchos::rcp(new panzer_stk::CubeTetMeshFactory);
  } else if (type == "Hex" && dim == 3) {
    mesh_factory = Teuchos::rcp(new panzer_stk::CubeHexMeshFactory);
  } else if (type == "Line" && dim == 1) {
    mesh_factory = Teuchos::rcp(new panzer_stk::LineMeshFactory);
  } else {
    // Throw an error
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
      "ERROR: Type and/or dimension of inline mesh not valid!");
  }

  // TODO BWR THIS SET PARAM LIST IS NOT WORKING??
  std::cout << pl->sublist("Mesh Factory Parameter List") << std::endl;
  mesh_factory->setParameterList(Teuchos::rcp(new Teuchos::ParameterList(pl->sublist("Mesh Factory Parameter List"))));
  auto mesh = mesh_factory->buildMesh(MPI_COMM_WORLD);

  return mesh;
}

Teuchos::RCP<panzer::WorksetContainer> createWorksetsForFields(Teuchos::RCP<panzer_stk::STK_Interface> mesh, 
                                                               Teuchos::RCP<panzer::DOFManager> dof_manager, 
                                                               std::map<std::string,panzer::BasisDescriptor> fmap)
{

  // Dof manager is, at least, a global indexer
  // and is passed in/out, if needed

  // And a needs map (block to needs)
  std::map<std::string, panzer::WorksetNeeds> needs_map;

  // Assuming one element block for now...
  std::vector<std::string> eblocks;
  mesh->getElementBlockNames(eblocks);
  {

    // Build a connectivity manager for the given mesh
    const auto conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    // Initialize the dof manager with the conn manager
    dof_manager = Teuchos::rcp(new panzer::DOFManager(conn_manager,MPI_COMM_WORLD));

    // Build basis and fields
    const auto & cell_topology = *mesh->getCellTopology(eblocks[0]);

    for (auto & map : fmap)
    {
      auto name = map.first;
      auto bd = map.second;

      auto intrepid_basis = panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>(bd.getType(),bd.getOrder(), cell_topology);
      auto field_pattern = Teuchos::rcp(new panzer::Intrepid2FieldPattern(intrepid_basis));

      // Add field to dof manager to lock in the field pattern
      dof_manager->addField(eblocks[0], name, field_pattern);

      needs_map[eblocks[0]].addBasis(bd);
      // TODO BWR needed?
      //panzer::IntegrationDescriptor integrator(3,panzer::IntegrationDescriptor::VOLUME);
      //needs_map[eblocks[0]].addIntegrator(integrator);
    }
  }

  // Finalize the dof manager
  dof_manager->buildGlobalUnknowns();

  // Build a workset container
  Teuchos::RCP<panzer::WorksetContainer> workset_container;
  {
    // Build workset factory
    auto factory = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh));

    workset_container = Teuchos::rcp(new panzer::WorksetContainer(factory, needs_map));
    workset_container->setGlobalIndexer(dof_manager);
  }

  return workset_container;

}

DynRankView getCoeffsForAnalyticFunction(Teuchos::RCP<panzer::PureBasis> basis, DynRankView cellNodes, 
                                         Teuchos::RCP<const std::vector<Intrepid2::Orientation> > orts /*, Functor ?? */ )
{
  typedef Intrepid2::Experimental::LagrangianInterpolation<PHX::Device> li;
  typedef Intrepid2::CellTools<PHX::Device> ct;

  auto numOwnedElems = cellNodes.extent(0);
  auto dim = basis->dimension();
  auto basisCardinality = basis->cardinality();

  auto coeffs = DynRankView("coeffs",numOwnedElems,basisCardinality);
  auto dofCoeffs = DynRankView("dofCoeffs",numOwnedElems,basisCardinality);
  auto dofCoords = DynRankView("dofCoords",numOwnedElems,basisCardinality,dim);

  // First, need to copy orientations to device
  auto orts_dev = Kokkos::DynRankView<Intrepid2::Orientation,PHX::Device>("orts_dev",orts->size());
  auto orts_host = Kokkos::create_mirror_view(orts_dev);
  for (size_t i=0; i < orts_host.extent(0); ++i)
    orts_host(i) = orts->at(i);
  Kokkos::deep_copy(orts_dev,orts_host);

  auto it2basis = basis->template getIntrepid2Basis<PHX::exec_space,double,double>();

  // get Dof coordinates and coefficients (reference frame)
  li::getDofCoordsAndCoeffs(dofCoords,dofCoeffs,it2basis.get(),orts_dev);

  auto funAtDofPoints = DynRankView("funAtDofPoints",numOwnedElems, basisCardinality);
    {
      auto physDofPoints = DynRankView("physDofPoints", numOwnedElems, basisCardinality, dim);
      // TODO BWR update this?
      ct::mapToPhysicalFrame(physDofPoints,dofCoords,cellNodes,it2basis->getBaseCellTopology());

      EvalSolFunctor<DynRankView> functor;
      functor.funAtPoints = funAtDofPoints;
      functor.points = physDofPoints;
      Kokkos::parallel_for("loop for evaluating the function at DoF points", numOwnedElems,functor);
      Kokkos::fence(); //make sure that funAtDofPoints has been evaluated
    }

    li::getBasisCoeffs(coeffs, funAtDofPoints, dofCoeffs);

    return coeffs;
}


////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////////////

typedef panzer::Traits::Residual ResidualType;

UNIT_TEST_GROUP(ResidualType)