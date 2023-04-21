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
#include "Panzer_DOFManager.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_OrientationsInterface.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_ProjectionTools.hpp"

typedef Kokkos::DynRankView<double,PHX::Device> DynRankView;

// for making explicit instantiated tests easier
#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(project_field,value,TYPE)

struct WorksetsAndOrts {

  Teuchos::RCP<std::vector<panzer::Workset> > worksets;
  Teuchos::RCP<const std::vector<Intrepid2::Orientation> > orientations;

};

Teuchos::RCP<panzer_stk::STK_Interface> createInlineMesh(Teuchos::RCP<Teuchos::ParameterList> pl);
WorksetsAndOrts getWorksetsAndOrtsForFields(
  Teuchos::RCP<panzer_stk::STK_Interface> mesh, 
  std::map<std::string,panzer::BasisDescriptor> fmap);
template <typename EvalType> 
bool checkProjection(Teuchos::RCP<panzer_stk::STK_Interface> mesh, 
                     std::map<std::string,panzer::BasisDescriptor> & fmap,
                     std::string & eShape);

struct Fun {

double
KOKKOS_INLINE_FUNCTION
operator()(const double& x, const double& y, const double& z, const size_t & dim) const {
  double f0 = sin(x*2)*sin(y*2)*sin(z*2)+sin(x*y*z*8);
  double f1 = cos(x*2)*cos(y*2)*cos(z*2);
  double f2 = cos(x*y*z*8);

  switch (dim){
    case 0:
      return f0; 
    case 1:
      return f1; 
    case 2:
      return f2; 
    default:
      return 0;
    }

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
      auto x = points(elem,i,0), y = points(elem,i,1), z = points(elem,i,2);
      if (funAtPoints.rank() == 3) {
        for(int d=0;d<static_cast<int>(funAtPoints.extent(2));d++)
          funAtPoints(elem,i,d) = fun(x,y,z,d); // vector basis
      } else {
        funAtPoints(elem,i) = fun(x,y,z,0); // scalar basis
      }
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

    GetCoeffsEvaluator(std::string & name, Teuchos::RCP<panzer::PureBasis> basis, 
                       std::string & eShape /*, Functor ??*/);

    void
    evaluateFields(
      typename Traits::EvalData d);

    void postRegistrationSetup(typename Traits::SetupData d, PHX::FieldManager<Traits> & fm);

  private:

    using ScalarT = typename EvalT::ScalarT;
    Teuchos::RCP<panzer::PureBasis> basis;
    PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS> coeffs;
    Teuchos::RCP<const std::vector<Intrepid2::Orientation> > orts;
    DynRankView cellNodes;
    std::string eShape;
  
}; // end of class

template<typename EvalT, typename Traits>
GetCoeffsEvaluator<EvalT, Traits>::
GetCoeffsEvaluator(std::string & name, Teuchos::RCP<panzer::PureBasis> basis,
                   std::string & eShape /*, Functor ?? */) 
  : basis(basis), eShape(eShape)
{

  // set up coeffs
  coeffs = PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS>(name, basis->functional);

  this->addEvaluatedField(coeffs);

  std::string n = "GetCoeffsEvaluator: " + name;
  this->setName(n);

}

template<typename EvalT,typename Traits>
void GetCoeffsEvaluator<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData  d, 
		      PHX::FieldManager<Traits>& /* fm */)
{
  orts = d.orientations_;
}

template<typename EvalT, typename Traits>
void
GetCoeffsEvaluator<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  workset)
{
  typedef Intrepid2::CellTools<PHX::Device> ct;
  typedef Intrepid2::FunctionSpaceTools<PHX::Device> fst;
  typedef Intrepid2::Experimental::ProjectionTools<PHX::Device> pts;

  // FYI, this all relies on a first-order mesh
  cellNodes = workset.getCellVertices().get_view(); // TODO BWR UPDATE 
  auto numNodesPerElem = cellNodes.extent(1);

  auto numOwnedElems = cellNodes.extent(0);
  TEUCHOS_ASSERT(numOwnedElems==orts->size());

  auto dim = basis->dimension();

  // First, need to copy orientations to device
  auto orts_dev = Kokkos::DynRankView<Intrepid2::Orientation,PHX::Device>("orts_dev",orts->size());
  auto orts_host = Kokkos::create_mirror_view(orts_dev);
  for (size_t i=0; i < orts_host.extent(0); ++i)
    orts_host(i) = orts->at(i);
  Kokkos::deep_copy(orts_dev,orts_host);

  auto it2basis = basis->template getIntrepid2Basis<PHX::exec_space,double,double>();
  auto functionSpace = it2basis->getFunctionSpace();
  auto cell_topology = it2basis->getBaseCellTopology(); // See note above

  {
    Intrepid2::Experimental::ProjectionStruct<PHX::Device,double> projStruct;
    projStruct.createL2ProjectionStruct(it2basis.get(), 3); // cubature order 3

    int numPoints = projStruct.getNumTargetEvalPoints();
    DynRankView evaluationPoints("evaluationPoints", numOwnedElems, numPoints, dim);

    pts::getL2EvaluationPoints(evaluationPoints,
        orts_dev,
        it2basis.get(),
        &projStruct);

    DynRankView refTargetAtEvalPoints, physTargetAtEvalPoints;
    if(functionSpace == Intrepid2::FUNCTION_SPACE_HCURL || functionSpace == Intrepid2::FUNCTION_SPACE_HDIV) {
      refTargetAtEvalPoints = DynRankView("targetAtEvalPoints", numOwnedElems, numPoints, dim);
      physTargetAtEvalPoints = DynRankView("targetAtEvalPoints", numOwnedElems, numPoints, dim);
    } else {
      refTargetAtEvalPoints = DynRankView("targetAtEvalPoints", numOwnedElems, numPoints);
      physTargetAtEvalPoints = DynRankView("targetAtEvalPoints", numOwnedElems, numPoints);
    }

    DynRankView physEvalPoints("physEvalPoints", numOwnedElems, numPoints, dim);
    {
      DynRankView linearBasisValuesAtEvalPoint("linearBasisValuesAtEvalPoint", numOwnedElems, numNodesPerElem);

      Kokkos::parallel_for(Kokkos::RangePolicy<PHX::Device::execution_space>(0,numOwnedElems),
          KOKKOS_LAMBDA (const int &i) {
        auto basisValuesAtEvalPoint = Kokkos::subview(linearBasisValuesAtEvalPoint,i,Kokkos::ALL());
        for(int j=0; j<numPoints; ++j){
          auto evalPoint = Kokkos::subview(evaluationPoints,i,j,Kokkos::ALL());
          if (eShape == "Hex") {
            Intrepid2::Impl::Basis_HGRAD_HEX_C1_FEM::template Serial<Intrepid2::OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
          } else if (eShape == "Tet") {
            Intrepid2::Impl::Basis_HGRAD_TET_C1_FEM::template Serial<Intrepid2::OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
          } else if (eShape == "Quad") {
            Intrepid2::Impl::Basis_HGRAD_QUAD_C1_FEM::template Serial<Intrepid2::OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
          } else if (eShape == "Tri") {
            Intrepid2::Impl::Basis_HGRAD_TRI_C1_FEM::template Serial<Intrepid2::OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
          }
          for(size_t k=0; k<numNodesPerElem; ++k)
            for(int d=0; d<dim; ++d)
              physEvalPoints(i,j,d) += cellNodes(i,k,d)*basisValuesAtEvalPoint(k);
        }
      });
      Kokkos::fence();
    }

    //transform the target function and its derivative to the reference element (inverse of pullback operator)
    DynRankView jacobian("jacobian", numOwnedElems, numPoints, dim, dim);
    DynRankView jacobian_det("jacobian_det", numOwnedElems, numPoints);
    DynRankView jacobian_inv("jacobian_inv", numOwnedElems, numPoints, dim, dim);
    ct::setJacobian(jacobian, evaluationPoints, cellNodes, cell_topology);
    ct::setJacobianDet (jacobian_det, jacobian);
    ct::setJacobianInv (jacobian_inv, jacobian);

    EvalSolFunctor<DynRankView> functor;
    functor.funAtPoints = physTargetAtEvalPoints;
    functor.points = physEvalPoints;
    Kokkos::parallel_for("loop for evaluating function in phys space", numOwnedElems, functor);
    Kokkos::fence();

    switch (functionSpace) {
    case Intrepid2::FUNCTION_SPACE_HGRAD:
      fst::mapHGradDataFromPhysToRef(refTargetAtEvalPoints,physTargetAtEvalPoints);
      break;
    case Intrepid2::FUNCTION_SPACE_HCURL:
      fst::mapHCurlDataFromPhysToRef(refTargetAtEvalPoints,jacobian,physTargetAtEvalPoints);
      break;
    case Intrepid2::FUNCTION_SPACE_HDIV:
      fst::mapHDivDataFromPhysToRef(refTargetAtEvalPoints,jacobian_inv,jacobian_det,physTargetAtEvalPoints);
      break;
    case Intrepid2::FUNCTION_SPACE_HVOL:
      fst::mapHVolDataFromPhysToRef(refTargetAtEvalPoints,jacobian_det,physTargetAtEvalPoints);
      break;
    default: {}
    }

    pts::getL2BasisCoeffs(coeffs.get_view(),
        refTargetAtEvalPoints,
        evaluationPoints,
        orts_dev,
        it2basis.get(),
        &projStruct);
  }

  return;
}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(project_field,value,EvalType)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  // Quad checks 
  //////////////////////////////////////////////////////////
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    std::string eShape("Quad");
    pl->set("Mesh Dimension",2);
    pl->set("Type","Quad");
    pl->sublist("Mesh Factory Parameter List").set("X Blocks",1);
    pl->sublist("Mesh Factory Parameter List").set("Y Blocks",1);
    pl->sublist("Mesh Factory Parameter List").set("X Elements",2);
    pl->sublist("Mesh Factory Parameter List").set("Y Elements",2);

    auto mesh = createInlineMesh(pl);
    {
      std::map<std::string,panzer::BasisDescriptor> fmap;
  
      fmap["MyField"] = panzer::BasisDescriptor(1,"HGrad");
      // return true if successful
      TEST_ASSERT(checkProjection<EvalType>(mesh,fmap,eShape));
    }
    {
      std::map<std::string,panzer::BasisDescriptor> fmap;
  
      fmap["MyField"] = panzer::BasisDescriptor(2,"HCurl");
      // return true if successful
      TEST_ASSERT(checkProjection<EvalType>(mesh,fmap,eShape));
    }
    {
      std::map<std::string,panzer::BasisDescriptor> fmap;
  
      fmap["MyField"] = panzer::BasisDescriptor(1,"HDiv");
      // return true if successful
      TEST_ASSERT(checkProjection<EvalType>(mesh,fmap,eShape));
    }
  }

  // Hex checks 
  //////////////////////////////////////////////////////////
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    std::string eShape("Hex");
    pl->set("Mesh Dimension",3);
    pl->set("Type","Hex");
    pl->sublist("Mesh Factory Parameter List").set("X Blocks",1);
    pl->sublist("Mesh Factory Parameter List").set("Y Blocks",1);
    pl->sublist("Mesh Factory Parameter List").set("Z Blocks",1);
    pl->sublist("Mesh Factory Parameter List").set("X Elements",2);
    pl->sublist("Mesh Factory Parameter List").set("Y Elements",2);
    pl->sublist("Mesh Factory Parameter List").set("Z Elements",2);

    auto mesh = createInlineMesh(pl);
    {
      std::map<std::string,panzer::BasisDescriptor> fmap;
  
      fmap["MyField"] = panzer::BasisDescriptor(1,"HGrad");
      // return true if successful
      TEST_ASSERT(checkProjection<EvalType>(mesh,fmap,eShape));
    }
    {
      std::map<std::string,panzer::BasisDescriptor> fmap;
  
      fmap["MyField"] = panzer::BasisDescriptor(1,"HCurl");
      // return true if successful
      TEST_ASSERT(checkProjection<EvalType>(mesh,fmap,eShape));
    }
    {
      std::map<std::string,panzer::BasisDescriptor> fmap;
  
      fmap["MyField"] = panzer::BasisDescriptor(1,"HDiv");
      // return true if successful
      TEST_ASSERT(checkProjection<EvalType>(mesh,fmap,eShape));
    }
  }
  
  // Tri checks 
  //////////////////////////////////////////////////////////
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    std::string eShape("Tri");
    pl->set("Mesh Dimension",2);
    pl->set("Type","Quad");
    pl->sublist("Mesh Factory Parameter List").set("X Blocks",1);
    pl->sublist("Mesh Factory Parameter List").set("Y Blocks",1);
    pl->sublist("Mesh Factory Parameter List").set("X Elements",2);
    pl->sublist("Mesh Factory Parameter List").set("Y Elements",2);

    auto mesh = createInlineMesh(pl);
    {
      std::map<std::string,panzer::BasisDescriptor> fmap;
  
      fmap["MyField"] = panzer::BasisDescriptor(1,"HGrad");
      // return true if successful
      TEST_ASSERT(checkProjection<EvalType>(mesh,fmap,eShape));
    }
    {
      std::map<std::string,panzer::BasisDescriptor> fmap;
  
      fmap["MyField"] = panzer::BasisDescriptor(1,"HCurl");
      // return true if successful
      TEST_ASSERT(checkProjection<EvalType>(mesh,fmap,eShape));
    }
    {
      std::map<std::string,panzer::BasisDescriptor> fmap;
  
      fmap["MyField"] = panzer::BasisDescriptor(1,"HDiv");
      // return true if successful
      TEST_ASSERT(checkProjection<EvalType>(mesh,fmap,eShape));
    }
  }

  // Tet checks
  //////////////////////////////////////////////////////////
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    std::string eShape("Tet");
    pl->set("Mesh Dimension",3);
    pl->set("Type","Hex");
    pl->sublist("Mesh Factory Parameter List").set("X Blocks",1);
    pl->sublist("Mesh Factory Parameter List").set("Y Blocks",1);
    pl->sublist("Mesh Factory Parameter List").set("Z Blocks",1);
    pl->sublist("Mesh Factory Parameter List").set("X Elements",2);
    pl->sublist("Mesh Factory Parameter List").set("Y Elements",2);
    pl->sublist("Mesh Factory Parameter List").set("Z Elements",2);

    auto mesh = createInlineMesh(pl);
    {
      std::map<std::string,panzer::BasisDescriptor> fmap;
  
      fmap["MyField"] = panzer::BasisDescriptor(1,"HGrad");
      // return true if successful
      TEST_ASSERT(checkProjection<EvalType>(mesh,fmap,eShape));
    }
    {
      std::map<std::string,panzer::BasisDescriptor> fmap;
  
      fmap["MyField"] = panzer::BasisDescriptor(1,"HCurl");
      // return true if successful
      TEST_ASSERT(checkProjection<EvalType>(mesh,fmap,eShape));
    }
    {
      std::map<std::string,panzer::BasisDescriptor> fmap;
  
      fmap["MyField"] = panzer::BasisDescriptor(1,"HDiv");
      // return true if successful
      TEST_ASSERT(checkProjection<EvalType>(mesh,fmap,eShape));
    }
  }
}

template <typename EvalType> 
bool checkProjection(Teuchos::RCP<panzer_stk::STK_Interface> mesh, 
                     std::map<std::string,panzer::BasisDescriptor> & fmap,
                     std::string & eShape)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  // This can be expanded, if needed but we assume only one field
  // and one block
  TEUCHOS_ASSERT(fmap.size()==1);
  auto fname = fmap.begin()->first;
  auto fbd = fmap.begin()->second; // BasisDescriptor
  auto type = fbd.getType();
  std::vector<std::string> eblocks;
  mesh->getElementBlockNames(eblocks);

  // set up worksets and orientations

  // For now, we make our lives easier and return one workset
  auto wkstsAndOrts = getWorksetsAndOrtsForFields(mesh,fmap);
  auto worksets = wkstsAndOrts.worksets;
  auto orientations = wkstsAndOrts.orientations;

  auto numCells = orientations->size();

  // we will project from first to a second order basis to a first
  // this essentially matches the intrepid2 test, but we've
  // wrapped things in an evaluator, etc.
  auto topo = mesh->getCellTopology(eblocks[0]);
  Teuchos::RCP<panzer::PureBasis> sourceBasis = Teuchos::rcp(new panzer::PureBasis(type,1,numCells,topo));
  Teuchos::RCP<panzer::PureBasis> targetBasis = Teuchos::rcp(new panzer::PureBasis(type,2,numCells,topo));

  Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm
   = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>);

  // need an evaluator to get coeffs of analytic function
  ///////////////////////////////////////////////////
  {
    RCP<panzer_stk::GetCoeffsEvaluator<EvalType,panzer::Traits> > e =
      rcp(new panzer_stk::GetCoeffsEvaluator<EvalType,panzer::Traits>(fname,sourceBasis,eShape)
      );

    fm->registerEvaluator<EvalType>(e);

    Teuchos::RCP<PHX::FieldTag> ft = e->evaluatedFields()[0];
    fm->requireField<EvalType>(*ft);

  }

  // add evaluator we're testing 
  ///////////////////////////////////////////////////
  {
    RCP<panzer_stk::ProjectField<EvalType,panzer::Traits> > e =
      rcp(new panzer_stk::ProjectField<EvalType,panzer::Traits>(fname,sourceBasis,targetBasis));

    fm->registerEvaluator<EvalType>(e);

    Teuchos::RCP<PHX::FieldTag> ft = e->evaluatedFields()[0];
    fm->requireField<EvalType>(*ft);
  }
  {
    std::string suffix("_final"); // otherwise we'd have two fields with same name
    RCP<panzer_stk::ProjectField<EvalType,panzer::Traits> > e =
      rcp(new panzer_stk::ProjectField<EvalType,panzer::Traits>(fname,targetBasis,sourceBasis,suffix));

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
  for ( auto & wkst : *worksets)
    fm->evaluateFields<EvalType>(wkst);
  fm->postEvaluate<EvalType>(0);

  typedef typename EvalType::ScalarT ScalarT;

  PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS> s(fname,sourceBasis->functional);
  PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS> t(fname+"_final",sourceBasis->functional);

  fm->getFieldData<EvalType>(s);
  fm->getFieldData<EvalType>(t);

  auto s_h = Kokkos::create_mirror_view(s.get_view());
  Kokkos::deep_copy(s_h, s.get_view());
  auto t_h = Kokkos::create_mirror_view(t.get_view());
  Kokkos::deep_copy(t_h, t.get_view());

  bool matched = true;
  for (size_t ncell=0;ncell<numCells;++ncell){
    for (int idx_pt=0;idx_pt<sourceBasis->cardinality();++idx_pt){
        if (std::abs(s_h(ncell,idx_pt) - t_h(ncell,idx_pt)) > 1e-14) matched = false;
    }
  }

  return matched;
}


// move and wrap inside some test utils namespace?
Teuchos::RCP<panzer_stk::STK_Interface> createInlineMesh(Teuchos::RCP<Teuchos::ParameterList> pl)
{

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
  } else {
    // Throw an error
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
      "ERROR: Type and/or dimension of inline mesh not valid!");
  }

  mesh_factory->setParameterList(Teuchos::rcp(new Teuchos::ParameterList(pl->sublist("Mesh Factory Parameter List"))));
  auto mesh = mesh_factory->buildMesh(MPI_COMM_WORLD);

  return mesh;
}

WorksetsAndOrts getWorksetsAndOrtsForFields(
  Teuchos::RCP<panzer_stk::STK_Interface> mesh, 
  std::map<std::string,panzer::BasisDescriptor> fmap)
{

  WorksetsAndOrts wksOrts;

  // Dof manager is, at least, a global indexer
  // and is passed in/out, if needed

  // And a needs map (block to needs)
  std::map<std::string, panzer::WorksetNeeds> needs_map;

  // Assuming one element block for now...
  std::vector<std::string> eblocks;
  mesh->getElementBlockNames(eblocks);

  Teuchos::RCP<panzer::DOFManager> dof_manager;
  {

    // Build a connectivity manager for the given mesh
    const auto conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    // Initialize the dof manager with the conn manager
    dof_manager = Teuchos::rcp(new panzer::DOFManager(conn_manager,MPI_COMM_WORLD));

    // Build basis and fields
    const auto & cell_topology = mesh->getCellTopology(eblocks[0]);
    // Set cell data
    needs_map[eblocks[0]].cellData = panzer::CellData(1,cell_topology); // numCells will get overwritten

    for (auto & map : fmap)
    {
      auto name = map.first;
      auto bd = map.second;

      auto intrepid_basis = panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>(bd.getType(),bd.getOrder(),cell_topology);
      auto field_pattern = Teuchos::rcp(new panzer::Intrepid2FieldPattern(intrepid_basis));

      // Add field to dof manager to lock in the field pattern
      dof_manager->addField(eblocks[0], name, field_pattern);

      needs_map[eblocks[0]].addBasis(bd);
    }
  }

  // Finalize the dof manager
  dof_manager->buildGlobalUnknowns();

  // Build workset factory
  auto factory = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh));

  auto workset_container = Teuchos::rcp(new panzer::WorksetContainer(factory, needs_map));
  workset_container->setGlobalIndexer(dof_manager);
  wksOrts.worksets = workset_container->getWorksets(panzer::WorksetDescriptor(eblocks[0],
                                                    panzer::WorksetSizeType::ALL_ELEMENTS, false, true));

  wksOrts.orientations = workset_container->getOrientations();

  return wksOrts;

}

////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////////////

typedef panzer::Traits::Residual ResidualType;

UNIT_TEST_GROUP(ResidualType)