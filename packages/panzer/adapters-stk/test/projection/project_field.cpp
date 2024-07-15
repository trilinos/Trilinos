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
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Panzer_STK_ProjectField.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Phalanx_FieldManager.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_ProjectionTools.hpp"

#include "EvaluatorTestTools.hpp"

typedef Kokkos::DynRankView<double,PHX::Device> DynRankView;

// for making explicit instantiated tests easier
#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(project_field,value,TYPE)

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
      auto x = points(elem,i,0);
      auto y = points(elem,i,1);
      double z = 0.;
      if (points.extent(2) == 3)
        z = points(elem,i,2);
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

//**********************************************************************
template<typename EvalT, typename Traits>
class GetCoeffsEvaluator
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    GetCoeffsEvaluator(std::string & name, Teuchos::RCP<panzer::PureBasis> basis, 
                       std::string & eShape, const size_t numNodesPerElem, const size_t dim);

    void
    evaluateFields(
      typename Traits::EvalData d);

    void postRegistrationSetup(typename Traits::SetupData d, PHX::FieldManager<Traits> & fm);

    enum ElemShape {HEX, TET, QUAD, TRI};

  private:

    using ScalarT = typename EvalT::ScalarT;
    Teuchos::RCP<panzer::PureBasis> basis;
    PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS> coeffs;
    Teuchos::RCP<const std::vector<Intrepid2::Orientation> > orts;
    Kokkos::DynRankView<Intrepid2::Orientation,PHX::Device> local_orts;
    Kokkos::DynRankView<double,PHX::Device> local_nodes;
    DynRankView local_jacobian, local_jacobian_inv, local_jacobian_det;
    DynRankView local_physEvalPoints;
    DynRankView local_refTargetAtEvalPoints, local_physTargetAtEvalPoints;
    ElemShape elemShape;
    Intrepid2::ProjectionStruct<PHX::Device,double> projStruct;
    Teuchos::RCP<Intrepid2::Basis<PHX::Device::execution_space,double,double> > it2basis;
  
}; // end of class

template<typename EvalT, typename Traits>
GetCoeffsEvaluator<EvalT, Traits>::
GetCoeffsEvaluator(std::string & name, Teuchos::RCP<panzer::PureBasis> basis,
                   std::string & eShape, const size_t numNodesPerElem, const size_t dim) 
  : basis(basis)
{

  // set up coeffs
  coeffs = PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS>(name, basis->functional);

  auto workset_size = basis->functional->extent(0); // Max number of cells

  this->addEvaluatedField(coeffs);

  std::string n = "GetCoeffsEvaluator: " + name;
  this->setName(n);

  if (eShape == "Hex") {
    elemShape = HEX;
  } else if (eShape == "Tet") {
    elemShape = TET;
  } else if (eShape == "Quad") {
    elemShape = QUAD;
  } else if (eShape == "Tri") {
    elemShape = TRI;
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
      "ERROR: Element shape not supported!");
  }

  it2basis = basis->template getIntrepid2Basis<PHX::exec_space,double,double>();
  
  // set up projection structure
  projStruct.createL2ProjectionStruct(it2basis.get(), 3); // cubature order 3
  int numPoints = projStruct.getNumTargetEvalPoints();

  // storage for local (to the workset) orientations and cell nodes
  local_orts  = Kokkos::DynRankView<Intrepid2::Orientation,PHX::Device>("orts",workset_size);
  local_nodes = Kokkos::DynRankView<double,PHX::Device>("cellNodes",workset_size,numNodesPerElem,dim);

  // storage for local objects that don't need to know the number of target points
  local_physEvalPoints = DynRankView("physEvalPoints",workset_size,numPoints,dim);
  local_jacobian = DynRankView("jacobian", workset_size, numPoints, dim, dim);
  local_jacobian_det = DynRankView("jacobian_det", workset_size, numPoints);
  local_jacobian_inv = DynRankView("jacobian_inv", workset_size, numPoints, dim, dim);
 
  if (basis->isVectorBasis()) {
    local_refTargetAtEvalPoints = DynRankView("targetAtEvalPoints",workset_size, numPoints,dim);
    local_physTargetAtEvalPoints = DynRankView("targetAtEvalPoints",workset_size, numPoints,dim);
  } else {
    local_refTargetAtEvalPoints = DynRankView("targetAtEvalPoints",workset_size, numPoints);
    local_physTargetAtEvalPoints = DynRankView("targetAtEvalPoints",workset_size, numPoints);
  }

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
  using ct = Intrepid2::CellTools<PHX::Device>;
  using fst = Intrepid2::FunctionSpaceTools<PHX::Device>;
  using pts = Intrepid2::ProjectionTools<PHX::Device>;

  // FYI, this all relies on a first-order mesh
  auto cellNodesAll = workset.getCellNodes();
  size_t numOwnedElems = workset.num_cells;
  TEUCHOS_ASSERT(local_nodes.extent(0)==local_orts.extent(0));

  auto dim = basis->dimension();
  auto numNodesPerElem = local_nodes.extent(1);

  // Get subview of local data in case num_cells doesn't match workset_size,
  // as can happen for the final workset
  const auto cell_range = std::pair<int,int>(0,numOwnedElems);
  auto sub_local_nodes = Kokkos::subview(local_nodes,cell_range,Kokkos::ALL(),Kokkos::ALL());
  auto sub_local_orts  = Kokkos::subview(local_orts,cell_range);

  auto sub_coeffs      = Kokkos::subview(coeffs.get_view(),cell_range,Kokkos::ALL());

  // First, need to copy orientations to device and grab local nodes
  auto orts_host = Kokkos::create_mirror_view(sub_local_orts);
  // subselect
  auto nodes_host = Kokkos::create_mirror_view(sub_local_nodes);
  auto nodesAll_host = Kokkos::create_mirror_view(cellNodesAll.get_view());
  for (size_t i=0; i < numOwnedElems; ++i) {
    orts_host(i) = orts->at(workset.cell_local_ids[i]);
    for (size_t j=0; j < numNodesPerElem; ++j)
      for (int k=0; k < dim; ++k)
        nodes_host(i,j,k) = nodesAll_host(i,j,k);
  }
  Kokkos::deep_copy(sub_local_orts,orts_host);
  Kokkos::deep_copy(sub_local_nodes,nodes_host);

  auto functionSpace = it2basis->getFunctionSpace();
  auto cell_topology = it2basis->getBaseCellTopology(); // See note above

  {

    auto evaluationPoints = projStruct.getAllEvalPoints();

    auto physEvalPoints = Kokkos::subview(local_physEvalPoints,cell_range,Kokkos::ALL(),Kokkos::ALL());
    ct::mapToPhysicalFrame(physEvalPoints,evaluationPoints,sub_local_nodes,it2basis->getBaseCellTopology());

    auto refTargetAtEvalPoints = Kokkos::subview(local_refTargetAtEvalPoints,cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto physTargetAtEvalPoints = Kokkos::subview(local_physTargetAtEvalPoints,cell_range,Kokkos::ALL(),Kokkos::ALL());

    // transform the target function and its derivative to the reference element (inverse of pullback operator)
    auto jacobian = Kokkos::subview(local_jacobian,cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    auto jacobian_det = Kokkos::subview(local_jacobian_det,cell_range,Kokkos::ALL());
    auto jacobian_inv = Kokkos::subview(local_jacobian_inv,cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    ct::setJacobian(jacobian, evaluationPoints, sub_local_nodes, cell_topology);
    ct::setJacobianDet (jacobian_det, jacobian);
    ct::setJacobianInv (jacobian_inv, jacobian);

    EvalSolFunctor<decltype(physTargetAtEvalPoints)> functor;
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

    pts::getL2BasisCoeffs(sub_coeffs,
        refTargetAtEvalPoints,
        sub_local_orts,
        it2basis.get(),
        &projStruct);
  }
  return;
}

template<typename EvalT, typename Traits>
class CheckCoeffsEvaluator
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    CheckCoeffsEvaluator(std::string & name, Teuchos::RCP<panzer::PureBasis> basis,
                         Teuchos::RCP<std::vector<double> > errors);

    void
    evaluateFields(
      typename Traits::EvalData d);

    void postRegistrationSetup(typename Traits::SetupData d, PHX::FieldManager<Traits> & fm);

  private:

    using ScalarT = typename EvalT::ScalarT;
    Teuchos::RCP<panzer::PureBasis> basis;
    PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS> before, after, diff;
    Teuchos::RCP<std::vector<double> > errors;
  
}; // end of class

template<typename EvalT, typename Traits>
CheckCoeffsEvaluator<EvalT, Traits>::
CheckCoeffsEvaluator(std::string & name, Teuchos::RCP<panzer::PureBasis> basis,
Teuchos::RCP<std::vector<double> > errors)
  : basis(basis), errors(errors)
{
  
  // set up fields

  before = PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS>(name,         basis->functional);
  after = PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS>( name+"_final",basis->functional);
  diff = PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS>("diff",basis->functional);
  this->addNonConstDependentField(before);
  this->addNonConstDependentField(after);
  this->addEvaluatedField(diff); // this is hack to get the evaluator to run

  std::string n = "CheckCoeffsEvaluator: " + name;
  this->setName(n);

}

template<typename EvalT,typename Traits>
void CheckCoeffsEvaluator<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData  d, 
		      PHX::FieldManager<Traits>& /* fm */)
{
}

template<typename EvalT, typename Traits>
void
CheckCoeffsEvaluator<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  workset)
{

  // all on host for checking
  auto numOwnedElems = workset.num_cells;

  auto s_h = Kokkos::create_mirror_view(before.get_view());
  Kokkos::deep_copy(s_h, before.get_view());
  auto t_h = Kokkos::create_mirror_view(after.get_view());
  Kokkos::deep_copy(t_h, after.get_view());

  double L2err = 0.;
  for (int ncell=0;ncell<numOwnedElems;++ncell){
    for (int idx_pt=0;idx_pt<basis->cardinality();++idx_pt){
      L2err += ( s_h(ncell,idx_pt) - t_h(ncell,idx_pt) ) * ( s_h(ncell,idx_pt) - t_h(ncell,idx_pt) );
    }
  }

  errors->push_back(std::sqrt(L2err));

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

    auto mesh = EvaluatorTestTools::createInlineMesh(pl);
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

    auto mesh = EvaluatorTestTools::createInlineMesh(pl);
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
    pl->set("Type","Tri");
    pl->sublist("Mesh Factory Parameter List").set("X Blocks",1);
    pl->sublist("Mesh Factory Parameter List").set("Y Blocks",1);
    pl->sublist("Mesh Factory Parameter List").set("X Elements",2);
    pl->sublist("Mesh Factory Parameter List").set("Y Elements",2);

    auto mesh = EvaluatorTestTools::createInlineMesh(pl);
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
    pl->set("Type","Tet");
    pl->sublist("Mesh Factory Parameter List").set("X Blocks",1);
    pl->sublist("Mesh Factory Parameter List").set("Y Blocks",1);
    pl->sublist("Mesh Factory Parameter List").set("Z Blocks",1);
    pl->sublist("Mesh Factory Parameter List").set("X Elements",2);
    pl->sublist("Mesh Factory Parameter List").set("Y Elements",2);
    pl->sublist("Mesh Factory Parameter List").set("Z Elements",2);

    auto mesh = EvaluatorTestTools::createInlineMesh(pl);
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

  // return worksets and orientations and get relevant sizes
  auto wkstsAndOrts = EvaluatorTestTools::getWorksetsAndOrtsForFields(mesh,fmap,2);
  auto worksets = wkstsAndOrts.worksets;
  auto orientations = wkstsAndOrts.orientations;
  auto wkstSize = (*worksets)[0].num_cells; // As long as we don't pick the last workset, this gives the max size 
  auto numNodesPerElem = (*worksets)[0].getCellNodes().extent(1);

  // we will project from first to a second order basis to a first
  // this essentially matches the intrepid2 test, but we've
  // wrapped things in an evaluator, etc.
  auto topo = mesh->getCellTopology(eblocks[0]);
  auto dim = topo->getDimension();

  Teuchos::RCP<panzer::PureBasis> sourceBasis = Teuchos::rcp(new panzer::PureBasis(type,1,wkstSize,topo));
  Teuchos::RCP<panzer::PureBasis> targetBasis = Teuchos::rcp(new panzer::PureBasis(type,2,wkstSize,topo));

  Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm
   = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>);

  // need an evaluator to get coeffs of analytic function
  ///////////////////////////////////////////////////
  {
    RCP<panzer_stk::GetCoeffsEvaluator<EvalType,panzer::Traits> > e =
      rcp(new panzer_stk::GetCoeffsEvaluator<EvalType,panzer::Traits>(fname,sourceBasis,eShape,numNodesPerElem,dim)
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
    std::string outName = fname+"_final"; // otherwise we'd have two fields with same name
    RCP<panzer_stk::ProjectField<EvalType,panzer::Traits> > e =
      rcp(new panzer_stk::ProjectField<EvalType,panzer::Traits>(fname,targetBasis,sourceBasis,outName));

    fm->registerEvaluator<EvalType>(e);

    Teuchos::RCP<PHX::FieldTag> ft = e->evaluatedFields()[0];
    fm->requireField<EvalType>(*ft);
  }

  // need an evaluator to test projection
  ///////////////////////////////////////////////////
  RCP<std::vector<double> > errors = rcp(new std::vector<double>);
  {
    RCP<panzer_stk::CheckCoeffsEvaluator<EvalType,panzer::Traits> > e =
      rcp(new panzer_stk::CheckCoeffsEvaluator<EvalType,panzer::Traits>(fname,sourceBasis,errors)
      );

    fm->registerEvaluator<EvalType>(e);

    // again a bit of a hack
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

  bool matched = true;

  for ( auto & err : *errors ) {
    if ( err > 1e-14 ) matched = false;
  }

  return matched;
}

////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////////////

typedef panzer::Traits::Residual ResidualType;

UNIT_TEST_GROUP(ResidualType)
