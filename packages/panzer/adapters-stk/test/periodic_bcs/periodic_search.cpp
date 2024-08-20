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

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Tuple.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_Utilities.hpp"

// TODO using this ALSO for global comm of info...
// If old functionality is removed, the tests will need to be 
// rewritten
#include "Panzer_STK_PeriodicBC_Matcher.hpp"
#include "Panzer_STK_PeriodicBC_Parser.hpp"
#include "Panzer_STK_PeriodicBC_MatchConditions.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_PauseToAttach.hpp"

#include "Kokkos_DynRankView.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HCURL_QUAD_I1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"

#include <string>

#ifdef PANZER_HAVE_STKSEARCH

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer {

  template <typename Intrepid2Type>
  RCP<const panzer::FieldPattern> buildFieldPattern()
  {
     // build a geometric pattern from a single basis
     RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2Type);
     RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
     return pattern;
  }

  TEUCHOS_UNIT_TEST(periodic_search, fillLocalSearchVector)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
      RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
      pl->set("X Blocks",2);
      pl->set("Y Blocks",1);
      pl->set("X Elements",2);
      pl->set("Y Elements",1);
      mesh_factory.setParameterList(pl);
      mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }
    TEST_ASSERT(mesh!=Teuchos::null);
    TEST_EQUALITY(mesh->getBulkData()->parallel_size(),2);

    auto myrank = mesh->getBulkData()->parallel_rank();

    panzer_stk::CoordMatcher x_matcher(0);
    panzer_stk::periodic_helpers::SphereIdVector coordsIds;
    auto error = x_matcher.getAbsoluteTolerance();
    panzer_stk::periodic_helpers::fillLocalSearchVector(*mesh,coordsIds,error,"top","coord");
   
    if (myrank==0) TEST_EQUALITY(coordsIds.size(),4);
    if (myrank==1) TEST_EQUALITY(coordsIds.size(),1);
  }

  TEUCHOS_UNIT_TEST(periodic_search, fillLocalSearchVector_edge)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
      RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
      pl->set("X Blocks",2);
      pl->set("Y Blocks",1);
      pl->set("X Elements",2);
      pl->set("Y Elements",1);
      mesh_factory.setParameterList(pl);
      mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }
    TEST_ASSERT(mesh!=Teuchos::null);
    TEST_EQUALITY(mesh->getBulkData()->parallel_size(),2);

    auto myrank = mesh->getBulkData()->parallel_rank();

    panzer_stk::CoordMatcher x_matcher(0);
    panzer_stk::periodic_helpers::SphereIdVector coordsIds;
    auto error = x_matcher.getAbsoluteTolerance();
    panzer_stk::periodic_helpers::fillLocalSearchVector(*mesh,coordsIds,error,"top","edge");
   
    if (myrank==0) TEST_EQUALITY(coordsIds.size(),2);
    if (myrank==1) TEST_EQUALITY(coordsIds.size(),2);
  }

  TEUCHOS_UNIT_TEST(periodic_search, fillLocalSearchVector_face)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::CubeHexMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("Z Blocks",1);
       pl->set("X Elements",2);
       pl->set("Y Elements",1);
       pl->set("Z Elements",1);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }
    TEST_ASSERT(mesh!=Teuchos::null);
    TEST_EQUALITY(mesh->getBulkData()->parallel_size(),2);

    auto myrank = mesh->getBulkData()->parallel_rank();

    panzer_stk::CoordMatcher x_matcher(0);
    panzer_stk::periodic_helpers::SphereIdVector coordsIds;
    auto error = x_matcher.getAbsoluteTolerance();
    panzer_stk::periodic_helpers::fillLocalSearchVector(*mesh,coordsIds,error,"top","face");

    if (myrank==0) TEST_EQUALITY(coordsIds.size(),2);
    if (myrank==1) TEST_EQUALITY(coordsIds.size(),2);
  }

  TEUCHOS_UNIT_TEST(periodic_search, computeGlobalCentroid)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::CubeHexMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("Z Blocks",1);
       pl->set("X Elements",2);
       pl->set("Y Elements",1);
       pl->set("Z Elements",1);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }
    TEST_ASSERT(mesh!=Teuchos::null);
    TEST_EQUALITY(mesh->getBulkData()->parallel_size(),2);

    std::vector<double> centroid = panzer_stk::periodic_helpers::computeGlobalCentroid(*mesh,"top");

    double exact[3] = {.5,1.,.5};

    for (size_t d=0; d<3; ++d) TEST_FLOATING_EQUALITY(centroid[d],exact[d],1e-12);

  }

  TEUCHOS_UNIT_TEST(periodic_search, fillLocalSearchVector_removeDuplicates)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::CubeHexMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("Z Blocks",1);
       pl->set("X Elements",2);
       pl->set("Y Elements",1);
       pl->set("Z Elements",1);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }
    TEST_ASSERT(mesh!=Teuchos::null);
    TEST_EQUALITY(mesh->getBulkData()->parallel_size(),2);

    panzer_stk::CoordMatcher x_matcher(0),y_matcher(1),z_matcher(2);
    panzer_stk::periodic_helpers::SphereIdVector topCoordsIds,leftCoordsIds,frontCoordsIds;
    panzer_stk::periodic_helpers::SphereIdVector uniqueLeftCoordsIds,uniqueFrontCoordsIds;
    auto error = x_matcher.getAbsoluteTolerance();

    // first get all the ids on each face
    panzer_stk::periodic_helpers::fillLocalSearchVector(*mesh,topCoordsIds,error,"top","coord");
    panzer_stk::periodic_helpers::fillLocalSearchVector(*mesh,leftCoordsIds,error,"left","coord");
    panzer_stk::periodic_helpers::fillLocalSearchVector(*mesh,frontCoordsIds,error,"front","coord");

    // now only get ids if they have not already been found
    std::vector<std::vector<std::string> > matchedSides(3);
    std::vector<panzer_stk::periodic_helpers::SearchId> doubleRequestsL, doubleRequestsF;
    matchedSides[0].push_back("top");
    panzer_stk::periodic_helpers::fillLocalSearchVector(*mesh,uniqueLeftCoordsIds,error,"left","coord",false,matchedSides[0],doubleRequestsL);
    matchedSides[0].push_back("left");
    panzer_stk::periodic_helpers::fillLocalSearchVector(*mesh,uniqueFrontCoordsIds,error,"front","coord",false,matchedSides[0],doubleRequestsF);

    // ensure that we do have shared ids and that we found them

    std::vector<size_t> topIds,leftIds,frontIds;
    std::vector<size_t> uniqueLeftIds, uniqueFrontIds;
    std::vector<size_t> doubleLfound, doubleFfound;

    // all ids
    for (size_t j=0; j<topCoordsIds.size(); ++j)
      topIds.push_back(topCoordsIds[j].second.id().id());
    for (size_t j=0; j<leftCoordsIds.size(); ++j)
      leftIds.push_back(leftCoordsIds[j].second.id().id());
    for (size_t j=0; j<frontCoordsIds.size(); ++j)
      frontIds.push_back(frontCoordsIds[j].second.id().id());
    // unique ids (when needed)
    for (size_t j=0; j<uniqueLeftCoordsIds.size(); ++j)
      uniqueLeftIds.push_back(uniqueLeftCoordsIds[j].second.id().id());
    for (size_t j=0; j<uniqueFrontCoordsIds.size(); ++j)
      uniqueFrontIds.push_back(uniqueFrontCoordsIds[j].second.id().id());

    // ids requested multiple times
    for (auto & id : doubleRequestsL) doubleLfound.push_back(id.id());
    for (auto & id : doubleRequestsF) doubleFfound.push_back(id.id());

    std::vector<size_t> TL_intersection, TF_intersection, TLF_intersection;
    std::vector<size_t> L_unique_exact, F_unique_exact;

    std::sort(topIds.begin(),topIds.end());
    std::sort(leftIds.begin(),leftIds.end());
    std::sort(frontIds.begin(),frontIds.end());

    std::set_intersection(topIds.begin(),topIds.end(),
                          leftIds.begin(),leftIds.end(),
                          std::back_inserter(TL_intersection));

    std::set_intersection(topIds.begin(),topIds.end(),
                          frontIds.begin(),frontIds.end(),
                          std::back_inserter(TF_intersection));

    int mySize = TL_intersection.size();
    int totalSize;

    Teuchos::reduceAll(*(mesh->getComm()),Teuchos::REDUCE_SUM,1,&mySize,&totalSize);

    TEST_ASSERT(totalSize>0);

    std::sort(doubleLfound.begin(),doubleLfound.end());

    // comparing directly to TL_intersection here is ok
    // ensure we found the nodes in the intersection

    int foundSize = doubleLfound.size();

    TEST_ASSERT(mySize==foundSize);
    size_t errL = 0;
    for (size_t i=0; i<doubleLfound.size(); ++i) {
      if (doubleLfound[i] != TL_intersection[i]) 
        ++errL;
    }
    TEST_EQUALITY(errL,0);

    // now ensure the unique nodes we found are correct
    std::set_difference(leftIds.begin(),leftIds.end(),
                        topIds.begin(),topIds.end(),
                        std::back_inserter(L_unique_exact));
    std::sort(L_unique_exact.begin(),L_unique_exact.end());
    std::sort(uniqueLeftIds.begin(),uniqueLeftIds.end());

    TEST_ASSERT(L_unique_exact.size()==uniqueLeftIds.size());
    size_t errLUnique = 0;
    for (size_t i=0; i<L_unique_exact.size(); ++i) {
      if (L_unique_exact[i] != uniqueLeftIds[i]) 
        ++errLUnique;
    }
    TEST_EQUALITY(errLUnique,0);

    std::sort(TL_intersection.begin(),TL_intersection.end());

    std::set_intersection(TL_intersection.begin(),TL_intersection.end(),
                          frontIds.begin(),frontIds.end(),
                          std::back_inserter(TLF_intersection));

    mySize = TLF_intersection.size();
    totalSize=0;

    Teuchos::reduceAll(*(mesh->getComm()),Teuchos::REDUCE_SUM,1,&mySize,&totalSize);

    TEST_ASSERT(totalSize>0);

    mySize = TF_intersection.size();
    totalSize=0;

    Teuchos::reduceAll(*(mesh->getComm()),Teuchos::REDUCE_SUM,1,&mySize,&totalSize);

    TEST_ASSERT(totalSize>0);

    // doubleFfound contains nodes on the front that have been previously matched
    // so this is the intersection of front with union of left and top

    std::vector<size_t> frontRepeatedNodes,LT_union; 
    std::set_union(leftIds.begin(),leftIds.end(),
                   topIds.begin(),topIds.end(),
                   std::back_inserter(LT_union));

    std::sort(LT_union.begin(),LT_union.end());

    std::set_intersection(LT_union.begin(),LT_union.end(),
                          frontIds.begin(),frontIds.end(),
                          std::back_inserter(frontRepeatedNodes));

    std::sort(doubleFfound.begin(),doubleFfound.end());
    std::sort(frontRepeatedNodes.begin(),frontRepeatedNodes.end());

    TEST_ASSERT(frontRepeatedNodes.size()==doubleFfound.size());
    size_t errF = 0;
    for (size_t i=0; i<doubleFfound.size(); ++i) {
      if (doubleFfound[i] != frontRepeatedNodes[i]) 
        ++errF;
    }
    TEST_EQUALITY(errF,0);

    // ensure the unique nodes we found are correct
    std::set_difference(frontIds.begin(),frontIds.end(),
                        LT_union.begin(),LT_union.end(),
                        std::back_inserter(F_unique_exact));
    std::sort(F_unique_exact.begin(),F_unique_exact.end());
    std::sort(uniqueFrontIds.begin(),uniqueFrontIds.end());

    TEST_ASSERT(F_unique_exact.size()==uniqueFrontIds.size());
    size_t errFUnique = 0;
    for (size_t i=0; i<F_unique_exact.size(); ++i) {
      if (F_unique_exact[i] != uniqueFrontIds[i]) 
        ++errFUnique;
    }
    TEST_EQUALITY(errFUnique,0);

    // now ensure that the final set of unique ids is really unique

    std::vector<size_t> finalIds;

    for (size_t j=0; j<topCoordsIds.size(); ++j)
      finalIds.push_back(topCoordsIds[j].second.id().id());
    for (size_t j=0; j<uniqueLeftCoordsIds.size(); ++j)
      finalIds.push_back(uniqueLeftCoordsIds[j].second.id().id());
    for (size_t j=0; j<uniqueFrontCoordsIds.size(); ++j)
      finalIds.push_back(uniqueFrontCoordsIds[j].second.id().id());

    std::sort(finalIds.begin(),finalIds.end());
    auto finalSize = finalIds.size();

    // remove duplicates, if there are any
    finalIds.erase(std::unique(finalIds.begin(),finalIds.end()),finalIds.end());
    auto uniqueSize = finalIds.size();

    TEST_ASSERT(uniqueSize==finalSize);

  }

  TEUCHOS_UNIT_TEST(periodic_search, transformLocalSearchVector_coordMatcher)
  {

     panzer_stk::CoordMatcher x_matcher(0);
     panzer_stk::CoordMatcher y_matcher(1);

     panzer_stk::periodic_helpers::SphereIdVector bottom, left;

     // create lines of points to be shifted

     size_t nPoints = 5;

     stk::mesh::EntityId id(0); // doesnt matter
     stk::mesh::EntityKey key(stk::topology::NODE_RANK,id); // doesnt matter
     panzer_stk::periodic_helpers::SearchId search_id(key,0); // doesnt matter
     for (size_t n=0; n<nPoints; ++n) {
       stk::search::Point<double> yCenter(0,n,0);
       stk::search::Point<double> xCenter(n,0,0); 
       left.emplace_back( stk::search::Sphere<double>(yCenter,1e-10), search_id );
       bottom.emplace_back( stk::search::Sphere<double>(xCenter,1e-10), search_id );
     }

     // centroids of lines parallel to the data with unity offset
     std::vector<double> rightCentroid = {1.,0.,0.}; 
     std::vector<double> topCentroid = {0.,1.,0.}; 

     panzer_stk::periodic_helpers::transformLocalSearchVector(left,y_matcher,rightCentroid);
     panzer_stk::periodic_helpers::transformLocalSearchVector(bottom,x_matcher,topCentroid);

     for (size_t n=0; n<nPoints; ++n){
       double exact[3] = {1.,(double)n,0.};
       auto ptShifted = left[n].first.center();
       for (size_t i=0; i<3; ++i) TEST_FLOATING_EQUALITY(exact[i],ptShifted[i],1e-14);
     }

     for (size_t n=0; n<nPoints; ++n){
       double exact[3] = {(double)n,1.,0.};
       auto ptShifted = bottom[n].first.center();
       for (size_t i=0; i<3; ++i) TEST_FLOATING_EQUALITY(exact[i],ptShifted[i],1e-14);
     }

  }

  TEUCHOS_UNIT_TEST(periodic_search, transformLocalSearchVector_planeMatcher)
  {

     panzer_stk::PlaneMatcher xy_matcher(0,1);
     panzer_stk::PlaneMatcher yz_matcher(1,2);
     panzer_stk::PlaneMatcher xz_matcher(0,2);

     panzer_stk::periodic_helpers::SphereIdVector xy, yz, xz;

     // create planes of points to be shifted

     size_t nPoints = 5;

     stk::mesh::EntityId id(0); // doesnt matter
     stk::mesh::EntityKey key(stk::topology::NODE_RANK,id); // doesnt matter
     panzer_stk::periodic_helpers::SearchId search_id(key,0); // doesnt matter
     for (size_t i=0; i<nPoints; ++i) {
       for (size_t j=0; j<nPoints; ++j) {
         stk::search::Point<double> xyCenter(i,j,0);
         stk::search::Point<double> yzCenter(0,i,j); 
         stk::search::Point<double> xzCenter(i,0,j); 
         xy.emplace_back( stk::search::Sphere<double>(xyCenter,1e-10), search_id );
         yz.emplace_back( stk::search::Sphere<double>(yzCenter,1e-10), search_id );
         xz.emplace_back( stk::search::Sphere<double>(xzCenter,1e-10), search_id );
       }
     }

     // centroids of planes parallel to the data with unity offset
     std::vector<double> xyCentroid = {2.,2.,1.}; 
     std::vector<double> yzCentroid = {1.,2.,2.}; 
     std::vector<double> xzCentroid = {2.,1.,2.}; 

     panzer_stk::periodic_helpers::transformLocalSearchVector(xy,xy_matcher,xyCentroid);
     panzer_stk::periodic_helpers::transformLocalSearchVector(yz,yz_matcher,yzCentroid);
     panzer_stk::periodic_helpers::transformLocalSearchVector(xz,xz_matcher,xzCentroid);

     size_t n = 0;

     for (size_t i=0; i<nPoints; ++i){
      for (size_t j=0; j<nPoints; ++j){
       double exact[3] = {(double)i,(double)j,1.};
       auto ptShifted = xy[n].first.center();
       ++n;
       for (size_t dim=0; dim<3; ++dim) TEST_FLOATING_EQUALITY(exact[dim],ptShifted[dim],1e-14); 
      }
     }

     n = 0;
     for (size_t i=0; i<nPoints; ++i){
        for (size_t j=0; j<nPoints; ++j){
           double exact[3] = {1.,(double)i,(double)j};
           auto ptShifted = yz[n].first.center();
           ++n;
           for (size_t dim=0; dim<3; ++dim) TEST_FLOATING_EQUALITY(exact[dim],ptShifted[dim],1e-14); 
        }
     }

     n = 0;
     for (size_t i=0; i<nPoints; ++i){
        for (size_t j=0; j<nPoints; ++j){
           double exact[3] = {(double)i,1.,(double)j};
           auto ptShifted = xz[n].first.center();
           ++n;
           for (size_t dim=0; dim<3; ++dim) TEST_FLOATING_EQUALITY(exact[dim],ptShifted[dim],1e-14); 
        }
     }

  }

  TEUCHOS_UNIT_TEST(periodic_search, transformLocalSearchVector_quarterPlaneMatcher)
  {

     // caps indicate the fixed coordinate
     // e.g. xyZ is quarter symmetry xz plane matched with yz plane
     panzer_stk::QuarterPlaneMatcher xyZ_matcher(0,1,2);
     panzer_stk::QuarterPlaneMatcher yzX_matcher(1,2,0);
     panzer_stk::QuarterPlaneMatcher xzY_matcher(0,2,1);
     panzer_stk::QuarterPlaneMatcher yxZ_matcher(1,0,2);

     panzer_stk::periodic_helpers::SphereIdVector yz,zx,zy,xz;

     // create planes of points to be shifted (these are side B's)

     size_t nPoints = 5;

     stk::mesh::EntityId id(0); // doesnt matter
     stk::mesh::EntityKey key(stk::topology::NODE_RANK,id); // doesnt matter
     panzer_stk::periodic_helpers::SearchId search_id(key,0); // doesnt matter
     for (size_t i=0; i<nPoints; ++i) {
       for (size_t j=0; j<nPoints; ++j) {
         stk::search::Point<double> yzCenter(0,i,j); 
         stk::search::Point<double> zxCenter(j,0,i);
         stk::search::Point<double> zyCenter(0,j,i); 
         stk::search::Point<double> xzCenter(i,0,j); 
         yz.emplace_back( stk::search::Sphere<double>(yzCenter,1e-10), search_id );
         zx.emplace_back( stk::search::Sphere<double>(zxCenter,1e-10), search_id );
         zy.emplace_back( stk::search::Sphere<double>(zyCenter,1e-10), search_id );
         xz.emplace_back( stk::search::Sphere<double>(xzCenter,1e-10), search_id );
       }
     }

     // centroids of planes (xz,yx,xy,yz) -- these are side A's
     std::vector<double> xzCentroid = {2.,0.,2.}; 
     std::vector<double> yxCentroid = {2.,2.,0.}; 
     std::vector<double> xyCentroid = {2.,2.,0.}; 
     std::vector<double> yzCentroid = {0.,2.,2.}; 

     panzer_stk::periodic_helpers::transformLocalSearchVector(yz,xyZ_matcher,xzCentroid);
     panzer_stk::periodic_helpers::transformLocalSearchVector(zx,yzX_matcher,yxCentroid);
     panzer_stk::periodic_helpers::transformLocalSearchVector(zy,xzY_matcher,xyCentroid);
     panzer_stk::periodic_helpers::transformLocalSearchVector(xz,yxZ_matcher,yzCentroid);

     size_t n = 0;

     for (size_t i=0; i<nPoints; ++i){
      for (size_t j=0; j<nPoints; ++j){
       // yz becomes xz
       double exact[3] = {(double)i,0.,(double)j};
       auto ptShifted = yz[n].first.center();
       ++n;
       for (size_t dim=0; dim<3; ++dim) TEST_FLOATING_EQUALITY(exact[dim],ptShifted[dim],1e-14); 
      }
     }

     n = 0;
     for (size_t i=0; i<nPoints; ++i){
        for (size_t j=0; j<nPoints; ++j){
           // zx becomes yx
           double exact[3] = {(double)j,(double)i,0.};
           auto ptShifted = zx[n].first.center();
           ++n;
           for (size_t dim=0; dim<3; ++dim) TEST_FLOATING_EQUALITY(exact[dim],ptShifted[dim],1e-14); 
        }
     }

     n = 0;
     for (size_t i=0; i<nPoints; ++i){
        for (size_t j=0; j<nPoints; ++j){
           // zy becomes xy
           double exact[3] = {(double)i,(double)j,0.};
           auto ptShifted = zy[n].first.center();
           ++n;
           for (size_t dim=0; dim<3; ++dim) TEST_FLOATING_EQUALITY(exact[dim],ptShifted[dim],1e-14); 
        }
     }

     n = 0;
     for (size_t i=0; i<nPoints; ++i){
        for (size_t j=0; j<nPoints; ++j){
           // xz becomes yz
           double exact[3] = {0.,(double)i,(double)j};
           auto ptShifted = xz[n].first.center();
           ++n;
           for (size_t dim=0; dim<3; ++dim) TEST_FLOATING_EQUALITY(exact[dim],ptShifted[dim],1e-14); 
        }
     }

  }
    TEUCHOS_UNIT_TEST(periodic_search, transformLocalSearchVector_wedgeMatcher)
  {

     // define wedge matchers based on the mirror plane
     std::vector<std::string> params; // need something... default is 1e-8, 3D
     panzer_stk::WedgeMatcher YZ_matcher(panzer_stk::WedgeMatcher::MirrorPlane::YZ_PLANE,params);
     panzer_stk::WedgeMatcher XZ_matcher(panzer_stk::WedgeMatcher::MirrorPlane::XZ_PLANE,params);

     panzer_stk::periodic_helpers::SphereIdVector YZ_sideB, XZ_sideB;

     // create planes of points to be shifted (these are side B's)

     size_t nPoints = 5;

     stk::mesh::EntityId id(0); // doesnt matter
     stk::mesh::EntityKey key(stk::topology::NODE_RANK,id); // doesnt matter
     panzer_stk::periodic_helpers::SearchId search_id(key,0); // doesnt matter
     // we will create planes with corners (0,0,0) (1,1,0) (1,1,1) (0,0,1) 
     for (size_t i=0; i<nPoints; ++i) {
       for (size_t j=0; j<nPoints; ++j) {
         double a = (double)i/(nPoints-1);
         double b = (double)j/(nPoints-1);
         double c = a;
         stk::search::Point<double> YZCenter(a,c,b); 
         stk::search::Point<double> XZCenter(c,a,b);
         YZ_sideB.emplace_back( stk::search::Sphere<double>(YZCenter,1e-10), search_id );
         XZ_sideB.emplace_back( stk::search::Sphere<double>(XZCenter,1e-10), search_id );
       }
     }

     // centroids of mirrored planes -- these are side A's
     std::vector<double> YZ_sideACentroid = {-0.5,0.5,0.5}; 
     std::vector<double> XZ_sideACentroid = {0.5,-0.5,0.5}; 

     panzer_stk::periodic_helpers::transformLocalSearchVector(YZ_sideB,YZ_matcher,YZ_sideACentroid);
     panzer_stk::periodic_helpers::transformLocalSearchVector(XZ_sideB,XZ_matcher,XZ_sideACentroid);

     size_t n = 0;

     for (size_t i=0; i<nPoints; ++i){
      for (size_t j=0; j<nPoints; ++j){
       // mirror over YZ (x --> -x)
       double a = -(double)i/(nPoints-1);
       double b = (double)j/(nPoints-1);
       double c = -a; // double negative...
       double exact[3] = {a,c,b};
       auto ptShifted = YZ_sideB[n].first.center();
       ++n;
       for (size_t dim=0; dim<3; ++dim) TEST_FLOATING_EQUALITY(exact[dim],ptShifted[dim],1e-14); 
      }
     }

     n = 0;
     for (size_t i=0; i<nPoints; ++i){
        for (size_t j=0; j<nPoints; ++j){
          // mirror over XZ (y --> -y)
          double a = -(double)i/(nPoints-1);
          double b = (double)j/(nPoints-1);
          double c = -a; // double negative...
          double exact[3] = {c,a,b};
          auto ptShifted = XZ_sideB[n].first.center();
          ++n;
          for (size_t dim=0; dim<3; ++dim) TEST_FLOATING_EQUALITY(exact[dim],ptShifted[dim],1e-14); 
        }
     }

  }

  TEUCHOS_UNIT_TEST(periodic_search, getGlobalPairing)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",6);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }
    TEST_ASSERT(mesh!=Teuchos::null);

    auto myrank = mesh->getBulkData()->parallel_rank();

    // run tests
    /////////////////////////////////////////////

    panzer_stk::CoordMatcher matcher(1);

    // next line requires communication
    Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds
          = panzer_stk::periodic_helpers::matchPeriodicSidesSearch("left","right",*mesh,matcher,"coord");
    Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds_edge
          = panzer_stk::periodic_helpers::matchPeriodicSidesSearch("left","right",*mesh,matcher,"edge");

    panzer_stk::periodic_helpers::SphereIdVector leftCoordsIds_edge,rightCoordsIds_edge;
    auto error = matcher.getAbsoluteTolerance();
    panzer_stk::periodic_helpers::fillLocalSearchVector(*mesh,leftCoordsIds_edge,error,"left","edge");
    panzer_stk::periodic_helpers::fillLocalSearchVector(*mesh,rightCoordsIds_edge,error,"right","edge");

    // TODO this is not the best... but using these old routines for global comms
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"left");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_right = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"right");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"left","edge");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_right = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"right","edge");
       std::vector<std::size_t> & sideIds = *idsAndCoords.first;
       std::vector<std::size_t> & sideIds_right = *idsAndCoords_right.first;
       std::vector<std::size_t> & sideIds_edge = *idsAndCoords_edge.first;
       std::vector<std::size_t> & sideIds_edge_right = *idsAndCoords_edge_right.first;
       std::vector<Tuple<double,3> > & sideCoords = *idsAndCoords.second;
       std::vector<Tuple<double,3> > & sideCoords_right = *idsAndCoords_right.second;
       std::vector<Tuple<double,3> > & sideCoords_edge = *idsAndCoords_edge.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_right = *idsAndCoords_edge_right.second;

    if(myrank==0) {
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          TEST_EQUALITY(pair.first,pair.second-12);
       }
       for(std::size_t i=0;i<globallyMatchedIds_edge->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds_edge)[i];

          // Get coordinates for matched edges
          Tuple<double,3> coord_left;
          Tuple<double,3> coord_right;
          int flag0 = 0;
          for(std::size_t j=0;j<sideIds_edge.size();j++){
            if(pair.first == sideIds_edge[j]){
              coord_left = sideCoords_edge[j];
              flag0++;
            }
            if(pair.second == sideIds_edge_right[j]){
              coord_right = sideCoords_edge_right[j];
              flag0++;
            }
          }
          TEST_EQUALITY(flag0,2);

          // Get coordinates of associated nodes
          int left_l  = -1;
          int left_u  = -1;
          int right_l = -1;
          int right_u = -1;
          int flag1   = 0;

          for(std::size_t j=0;j<sideCoords.size();j++) {
            if ((std::abs(sideCoords[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords[j][1] - (coord_left[1]-1.0/8.0)) < 1e-14)){
               left_l = j;
               flag1++;
            }
            if ((std::abs(sideCoords[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords[j][1] - (coord_left[1]+1.0/8.0)) < 1e-14)){
               left_u = j;
               flag1++;
            }
            if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]-1.0/8.0)) < 1e-14)){
               right_l = j;
               flag1++;
            }
            if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]+1.0/8.0)) < 1e-14)){
               right_u = j;
               flag1++;
            }
          }
          TEST_EQUALITY(flag1,4);

          // Test equivalence of node numbers
          TEST_EQUALITY(sideIds[left_l],sideIds_right[right_l]-12);
          TEST_EQUALITY(sideIds[left_u],sideIds_right[right_u]-12);

       }
    }
    else {
       TEST_EQUALITY(globallyMatchedIds->size(),0);
    }
  }

  TEUCHOS_UNIT_TEST(periodic_search, PeriodicBC_Matcher_multi_2x)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup 2D mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",4);
       pl->set("Y Elements",2);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }
    TEST_ASSERT(mesh!=Teuchos::null);

    {
       panzer_stk::CoordMatcher xmatcher(0);
       panzer_stk::CoordMatcher ymatcher(1);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> tb_Match
             = panzer_stk::buildPeriodicBC_Matcher("top","bottom",xmatcher);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> lr_Match
             = panzer_stk::buildPeriodicBC_Matcher("left","right",ymatcher);

       RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds;
       std::vector<std::string> matchedSides;
       globallyMatchedIds = tb_Match->getMatchedPair(*mesh,matchedSides);
       matchedSides.push_back("top");
       globallyMatchedIds = lr_Match->getMatchedPair(*mesh,matchedSides,globallyMatchedIds);

       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          out << pair.first << " " << pair.second << " " << std::endl;
       }

       // check
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];

          // double match here
          if(pair.first==1 || pair.first==19)
          {   TEST_EQUALITY(pair.second,9); }
          // remaining l/r match
          else if(pair.first==10)
          {   TEST_EQUALITY(pair.second,18); }
          // top/bottom spacing is 18
          else
          {   TEST_EQUALITY(pair.second,pair.first-18); }
       }

    }
  }

  TEUCHOS_UNIT_TEST(periodic_search, PeriodicBC_Matcher_multi_3x)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::CubeHexMeshFactory mesh_factory;

    // setup 3D mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("Z Blocks",1);
       pl->set("X Elements",4);
       pl->set("Y Elements",2);
       pl->set("Z Elements",2);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }
    TEST_ASSERT(mesh!=Teuchos::null);

    {
      panzer_stk::PlaneMatcher xy_matcher(0,1);
      panzer_stk::PlaneMatcher yz_matcher(1,2);
      panzer_stk::PlaneMatcher xz_matcher(0,2);

       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> tb_Match
             = panzer_stk::buildPeriodicBC_Matcher("top","bottom",xz_matcher);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> lr_Match
             = panzer_stk::buildPeriodicBC_Matcher("left","right",yz_matcher);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> fb_Match
             = panzer_stk::buildPeriodicBC_Matcher("front","back",xy_matcher);

       RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds;
       std::vector<std::string> matchedSides;
       globallyMatchedIds = tb_Match->getMatchedPair(*mesh,matchedSides);
       matchedSides.push_back("top");
       globallyMatchedIds = lr_Match->getMatchedPair(*mesh,matchedSides,globallyMatchedIds);
       matchedSides.push_back("left");
       globallyMatchedIds = fb_Match->getMatchedPair(*mesh,matchedSides,globallyMatchedIds);

       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          out << pair.first << " " << pair.second << " " << std::endl;
       }

       // check
       // triply periodic node goes from 73 to 55 to 63 to 9
       // doubly periodic edge (top left) goes from (19,46) to (1,28) to (9,36)
       // doubly periodic edge (top front) goes from (74-81) to (56-63) to (2-9)
       // doubly periodic edge (left front) goes from (64,55) to (72,63) to (18,9)
       // top to bottom spacing is -18
       // left to right spacing is 8
       // front to back spacing is -54
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];

          if(pair.first==73 || pair.first==55 || pair.first==63 || pair.first==19)
          {   TEST_EQUALITY(pair.second,9); }
          else if(pair.first==46)
          {   TEST_EQUALITY(pair.second,36); }
          else if(pair.first==64)
          {   TEST_EQUALITY(pair.second,18); }
          else if(pair.first>=74 && pair.first<=81)
          {   TEST_EQUALITY(pair.first-pair.second,72); }
          else
          {   size_t test=0;
              // one and only one should be satisfied
              size_t diff = (pair.first >= pair.second) ? 
                pair.first - pair.second : pair.second - pair.first;
              if ( diff == 18) ++test;
              if ( diff == 8 ) ++test;
              if ( diff == 54) ++test;
              TEST_EQUALITY(test,1); }
       }

    }
  }

  TEUCHOS_UNIT_TEST(periodic_search, matchPeriodicSides)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",6);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }
    TEST_ASSERT(mesh!=Teuchos::null);

    auto myrank = mesh->getBulkData()->parallel_rank();

       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_left = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"left");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_right = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"right");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_left = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"left","edge");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_right = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"right","edge");
       std::vector<std::size_t> & sideIds_left = *idsAndCoords_left.first;
       std::vector<std::size_t> & sideIds_right = *idsAndCoords_right.first;
       std::vector<std::size_t> & sideIds_edge_left = *idsAndCoords_edge_left.first;
       std::vector<std::size_t> & sideIds_edge_right = *idsAndCoords_edge_right.first;
       std::vector<Tuple<double,3> > & sideCoords_left = *idsAndCoords_left.second;
       std::vector<Tuple<double,3> > & sideCoords_right = *idsAndCoords_right.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_left = *idsAndCoords_edge_left.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_right = *idsAndCoords_edge_right.second;

       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_top = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"top");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_bottom = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"bottom");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_top = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"top","edge");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_bottom = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"bottom","edge");
       std::vector<std::size_t> & sideIds_top = *idsAndCoords_top.first;
       std::vector<std::size_t> & sideIds_bottom = *idsAndCoords_bottom.first;
       std::vector<std::size_t> & sideIds_edge_top = *idsAndCoords_edge_top.first;
       std::vector<std::size_t> & sideIds_edge_bottom = *idsAndCoords_edge_bottom.first;
       std::vector<Tuple<double,3> > & sideCoords_top = *idsAndCoords_top.second;
       std::vector<Tuple<double,3> > & sideCoords_bottom = *idsAndCoords_bottom.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_top = *idsAndCoords_edge_top.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_bottom = *idsAndCoords_edge_bottom.second;

    // run tests

    // Nodes
    {
       panzer_stk::CoordMatcher matcher(1);
       Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds
             = panzer_stk::periodic_helpers::matchPeriodicSidesSearch("left","right",*mesh,matcher);

       // match left & right sides
       if(myrank==0) {
          for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
             std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
             TEST_EQUALITY(pair.first,pair.second-12);
          }
       }
       else
          TEST_EQUALITY(globallyMatchedIds->size(),0);
    }

    // Edges
    {
       panzer_stk::CoordMatcher matcher(1);
       Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds_edge
             = panzer_stk::periodic_helpers::matchPeriodicSidesSearch("left","right",*mesh,matcher,"edge");

       // match left & right sides
       if(myrank==0) {
          for(std::size_t i=0;i<globallyMatchedIds_edge->size();i++) {
             std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds_edge)[i];

             // Get coordinates for matched edges
             Tuple<double,3> coord_left;
             Tuple<double,3> coord_right;
             int flag0 = 0;
             for(std::size_t j=0;j<sideIds_edge_left.size();j++){
               if(pair.first == sideIds_edge_left[j]){
                 coord_left = sideCoords_edge_left[j];
                 flag0++;
               }
               if(pair.second == sideIds_edge_right[j]){
                 coord_right = sideCoords_edge_right[j];
                 flag0++;
               }
             }
             TEST_EQUALITY(flag0,2);

             // Get coordinates of associated nodes
             int left_l  = -1;
             int left_u  = -1;
             int right_l = -1;
             int right_u = -1;
             int flag1   = 0;

             for(std::size_t j=0;j<sideCoords_left.size();j++) {
               if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]-1.0/8.0)) < 1e-14)){
                  left_l = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]+1.0/8.0)) < 1e-14)){
                  left_u = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]-1.0/8.0)) < 1e-14)){
                  right_l = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]+1.0/8.0)) < 1e-14)){
                  right_u = j;
                  flag1++;
               }
             }
             TEST_EQUALITY(flag1,4);

             // Test equivalence of node numbers
             TEST_EQUALITY(sideIds_left[left_l],sideIds_right[right_l]-12);
             TEST_EQUALITY(sideIds_left[left_u],sideIds_right[right_u]-12);
          }
       }
       else
          TEST_EQUALITY(globallyMatchedIds_edge->size(),0);
    }

    // Nodes
    {
       panzer_stk::CoordMatcher matcher(0);
       Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds
             = panzer_stk::periodic_helpers::matchPeriodicSidesSearch("top","bottom",*mesh,matcher);

       Teuchos::RCP<std::vector<std::size_t> > locallyRequiredIds = panzer_stk::periodic_helpers::getLocalSideIds(*mesh,"top");

       TEST_EQUALITY(globallyMatchedIds->size(),locallyRequiredIds->size());

       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          TEST_EQUALITY(pair.first,pair.second+52);
       }
    }

    // Edges
    {
       panzer_stk::CoordMatcher matcher(0);
       Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds_edge
             = panzer_stk::periodic_helpers::matchPeriodicSidesSearch("top","bottom",*mesh,matcher,"edge");

       Teuchos::RCP<std::vector<std::size_t> > locallyRequiredIds = panzer_stk::periodic_helpers::getLocalSideIds(*mesh,"top","edge");

       TEST_EQUALITY(globallyMatchedIds_edge->size(),locallyRequiredIds->size());

       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds_edge->size();i++) {
             std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds_edge)[i];

             // Get coordinates for matched edges
             Tuple<double,3> coord_top;
             Tuple<double,3> coord_bottom;
             int flag0 = 0;
             for(std::size_t j=0;j<sideIds_edge_top.size();j++){
               if(pair.first == sideIds_edge_top[j]){
                 coord_top = sideCoords_edge_top[j];
                 flag0++;
               }
               if(pair.second == sideIds_edge_bottom[j]){
                 coord_bottom = sideCoords_edge_bottom[j];
                 flag0++;
               }
             }
             TEST_EQUALITY(flag0,2);

             // Get coordinates of associated nodes
             int top_l  = -1;
             int top_u  = -1;
             int bottom_l = -1;
             int bottom_u = -1;
             int flag1   = 0;

             for(std::size_t j=0;j<sideCoords_top.size();j++) {
               if ((std::abs(sideCoords_top[j][0] - (coord_top[0]-1.0/24.0)) < 1e-14) && (std::abs(sideCoords_top[j][1] - coord_top[1]) < 1e-14)){
                  top_l = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_top[j][0] - (coord_top[0]+1.0/24.0)) < 1e-14) && (std::abs(sideCoords_top[j][1] - coord_top[1]) < 1e-14)){
                  top_u = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_bottom[j][0] - (coord_bottom[0]-1.0/24.0)) < 1e-14) && (std::abs(sideCoords_bottom[j][1] - coord_bottom[1]) < 1e-14)){
                  bottom_l = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_bottom[j][0] - (coord_bottom[0]+1.0/24.0)) < 1e-14) && (std::abs(sideCoords_bottom[j][1] - coord_bottom[1]) < 1e-14)){
                  bottom_u = j;
                  flag1++;
               }
             }
             TEST_EQUALITY(flag1,4);

             // Test equivalence of node numbers
             TEST_EQUALITY(sideIds_top[top_l],sideIds_bottom[bottom_l]+52);
             TEST_EQUALITY(sideIds_top[top_u],sideIds_bottom[bottom_u]+52);
       }
    }


    // test a failure case!
    {
       panzer_stk::CoordMatcher matcherX(0);
       panzer_stk::CoordMatcher matcherY(1);

       TEST_THROW(panzer_stk::periodic_helpers::matchPeriodicSidesSearch("top","right",*mesh,matcherY),std::logic_error);
       TEST_THROW(panzer_stk::periodic_helpers::matchPeriodicSidesSearch("top","right",*mesh,matcherX),std::logic_error);
       TEST_THROW(panzer_stk::periodic_helpers::matchPeriodicSidesSearch("bottom","left",*mesh,matcherY),std::logic_error);
       TEST_THROW(panzer_stk::periodic_helpers::matchPeriodicSidesSearch("bottom","left",*mesh,matcherX),std::logic_error);
       TEST_THROW(panzer_stk::periodic_helpers::matchPeriodicSidesSearch("top","right",*mesh,matcherY,"edge"),std::logic_error);
       TEST_THROW(panzer_stk::periodic_helpers::matchPeriodicSidesSearch("top","right",*mesh,matcherX,"edge"),std::logic_error);
       TEST_THROW(panzer_stk::periodic_helpers::matchPeriodicSidesSearch("bottom","left",*mesh,matcherY,"edge"),std::logic_error);
       TEST_THROW(panzer_stk::periodic_helpers::matchPeriodicSidesSearch("bottom","left",*mesh,matcherX,"edge"),std::logic_error);
    }
  }

  TEUCHOS_UNIT_TEST(periodic_search, PeriodicBC_Matcher)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",6);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }
    TEST_ASSERT(mesh!=Teuchos::null);

       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_top = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"top");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_bottom = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"bottom");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_top = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"top","edge");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_bottom = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"bottom","edge");
       std::vector<std::size_t> & sideIds_top = *idsAndCoords_top.first;
       std::vector<std::size_t> & sideIds_bottom = *idsAndCoords_bottom.first;
       std::vector<std::size_t> & sideIds_edge_top = *idsAndCoords_edge_top.first;
       std::vector<std::size_t> & sideIds_edge_bottom = *idsAndCoords_edge_bottom.first;
       std::vector<Tuple<double,3> > & sideCoords_top = *idsAndCoords_top.second;
       std::vector<Tuple<double,3> > & sideCoords_bottom = *idsAndCoords_bottom.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_top = *idsAndCoords_edge_top.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_bottom = *idsAndCoords_edge_bottom.second;

    // Nodes
    {
       panzer_stk::CoordMatcher matcher(0);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> pMatch
             = panzer_stk::buildPeriodicBC_Matcher("top","bottom",matcher);

       TEST_EQUALITY(pMatch->getLeftSidesetName(),std::string("top"));
       TEST_EQUALITY(pMatch->getRightSidesetName(),std::string("bottom"));
       const auto* check_cast_coord = pMatch->getAs<panzer_stk::PeriodicBC_Matcher<panzer_stk::CoordMatcher>>();
       const auto* check_cast_plane = pMatch->getAs<panzer_stk::PeriodicBC_Matcher<panzer_stk::PlaneMatcher>>();
       const auto* check_cast_wedge = pMatch->getAs<panzer_stk::PeriodicBC_Matcher<panzer_stk::WedgeMatcher>>();
       TEST_ASSERT(check_cast_coord != nullptr);
       TEST_ASSERT(check_cast_plane == nullptr);
       TEST_ASSERT(check_cast_wedge == nullptr);
       TEST_EQUALITY(check_cast_coord->getMatcher().getIndex(),0);

       std::vector<std::string> matchedSides;
       RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds = pMatch->getMatchedPair(*mesh,matchedSides);

       // for testing purposes!
       RCP<std::vector<std::size_t> > locallyRequiredIds = panzer_stk::periodic_helpers::getLocalSideIds(*mesh,"top");

       TEST_EQUALITY(globallyMatchedIds->size(),locallyRequiredIds->size());

       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          TEST_EQUALITY(pair.first,pair.second+52);
       }
    }

    // Edges
    {
       panzer_stk::CoordMatcher matcher(0);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> pMatch
             = panzer_stk::buildPeriodicBC_Matcher("top","bottom",matcher,"edge");

       std::vector<std::string> matchedSides;
       RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds_edge = pMatch->getMatchedPair(*mesh,matchedSides);

       // for testing purposes!
       RCP<std::vector<std::size_t> > locallyRequiredIds = panzer_stk::periodic_helpers::getLocalSideIds(*mesh,"top","edge");

       TEST_EQUALITY(globallyMatchedIds_edge->size(),locallyRequiredIds->size());

       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds_edge->size();i++) {
             std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds_edge)[i];

             // Get coordinates for matched edges
             Tuple<double,3> coord_top;
             Tuple<double,3> coord_bottom;
             int flag0 = 0;
             for(std::size_t j=0;j<sideIds_edge_top.size();j++){
               if(pair.first == sideIds_edge_top[j]){
                 coord_top = sideCoords_edge_top[j];
                 flag0++;
               }
               if(pair.second == sideIds_edge_bottom[j]){
                 coord_bottom = sideCoords_edge_bottom[j];
                 flag0++;
               }
             }
             TEST_EQUALITY(flag0,2);

             // Get coordinates of associated nodes
             int top_l  = -1;
             int top_u  = -1;
             int bottom_l = -1;
             int bottom_u = -1;
             int flag1   = 0;

             for(std::size_t j=0;j<sideCoords_top.size();j++) {
               if ((std::abs(sideCoords_top[j][0] - (coord_top[0]-1.0/24.0)) < 1e-14) && (std::abs(sideCoords_top[j][1] - coord_top[1]) < 1e-14)){
                  top_l = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_top[j][0] - (coord_top[0]+1.0/24.0)) < 1e-14) && (std::abs(sideCoords_top[j][1] - coord_top[1]) < 1e-14)){
                  top_u = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_bottom[j][0] - (coord_bottom[0]-1.0/24.0)) < 1e-14) && (std::abs(sideCoords_bottom[j][1] - coord_bottom[1]) < 1e-14)){
                  bottom_l = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_bottom[j][0] - (coord_bottom[0]+1.0/24.0)) < 1e-14) && (std::abs(sideCoords_bottom[j][1] - coord_bottom[1]) < 1e-14)){
                  bottom_u = j;
                  flag1++;
               }
             }
             TEST_EQUALITY(flag1,4);

             // Test equivalence of node numbers
             TEST_EQUALITY(sideIds_top[top_l],sideIds_bottom[bottom_l]+52);
             TEST_EQUALITY(sideIds_top[top_u],sideIds_bottom[bottom_u]+52);
       }
    }

    // test a failure case!
    {
       panzer_stk::CoordMatcher matcherX(0);
       panzer_stk::CoordMatcher matcherY(1);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> pMatch;
       std::vector<std::string> matchedSides;

       pMatch = panzer_stk::buildPeriodicBC_Matcher("bottom","left",matcherX);
       TEST_THROW(pMatch->getMatchedPair(*mesh,matchedSides),std::logic_error);

       pMatch = panzer_stk::buildPeriodicBC_Matcher("top","right",matcherX);
       TEST_THROW(pMatch->getMatchedPair(*mesh,matchedSides),std::logic_error);

       pMatch = panzer_stk::buildPeriodicBC_Matcher("top","right",matcherY);
       TEST_THROW(pMatch->getMatchedPair(*mesh,matchedSides),std::logic_error);

       pMatch = panzer_stk::buildPeriodicBC_Matcher("bottom","left",matcherY);
       TEST_THROW(pMatch->getMatchedPair(*mesh,matchedSides),std::logic_error);

       pMatch = panzer_stk::buildPeriodicBC_Matcher("bottom","left",matcherX,"edge");
       TEST_THROW(pMatch->getMatchedPair(*mesh,matchedSides),std::logic_error);

       pMatch = panzer_stk::buildPeriodicBC_Matcher("top","right",matcherX,"edge");
       TEST_THROW(pMatch->getMatchedPair(*mesh,matchedSides),std::logic_error);

       pMatch = panzer_stk::buildPeriodicBC_Matcher("top","right",matcherY,"edge");
       TEST_THROW(pMatch->getMatchedPair(*mesh,matchedSides),std::logic_error);

       pMatch = panzer_stk::buildPeriodicBC_Matcher("bottom","left",matcherY,"edge");
       TEST_THROW(pMatch->getMatchedPair(*mesh,matchedSides),std::logic_error);
    }

  }

  //TEUCHOS_UNIT_TEST(periodic_search, PeriodicBC_Matcher_relative)
  //{
  //  using Teuchos::RCP;
  //  using Teuchos::Tuple;


  //  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  //  panzer_stk::SquareQuadMeshFactory mesh_factory;

  //  // setup mesh
  //  /////////////////////////////////////////////
  //  RCP<panzer_stk::STK_Interface> mesh;
  //  {
  //     // make a mesh with small length-scale
  //     RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
  //     pl->set("X Blocks",2);
  //     pl->set("Y Blocks",1);
  //     pl->set("X Elements",6);
  //     pl->set("Y Elements",4);
  //     pl->set("X0",0.0);
  //     pl->set("Xf",1.0e-6);
  //     pl->set("Y0",0.0);
  //     pl->set("Yf",1.0e-6);
  //     mesh_factory.setParameterList(pl);
  //     mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
  //  }

  //     std::pair<RCP<std::vector<std::size_t> >,
  //               RCP<std::vector<Tuple<double,3> > > > idsAndCoords_top = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"top");
  //     std::pair<RCP<std::vector<std::size_t> >,
  //               RCP<std::vector<Tuple<double,3> > > > idsAndCoords_bottom = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"bottom");

  //  // Nodes
  //  {
  //     // set up a matcher with a tolerance of 1e-6
  //     std::vector<std::string> params;
  //     params.push_back("1e-6");
  //     CoordMatcher bad_matcher(0,params);
  //     Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> bad_pMatch
  //           = panzer_stk::buildPeriodicBC_Matcher("top","bottom",bad_matcher);

  //     // matching should fail since the tolerance is larger than the mesh size
  //     TEST_THROW(bad_pMatch->getMatchedPair(*mesh),std::logic_error);

  //     // make the tolerance relative, then matching shouldn't fail
  //     params.push_back("relative");
  //     CoordMatcher matcher(0,params);
  //     Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> pMatch
  //           = panzer_stk::buildPeriodicBC_Matcher("top","bottom",matcher);

  //     RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds = pMatch->getMatchedPair(*mesh);

  //     // for testing purposes!
  //     RCP<std::vector<std::size_t> > locallyRequiredIds = panzer_stk::periodic_helpers::getLocalSideIds(*mesh,"top");

  //     TEST_EQUALITY(globallyMatchedIds->size(),locallyRequiredIds->size());

  //     // match top & bottom sides
  //     for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
  //        std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
  //        TEST_EQUALITY(pair.first,pair.second+52);
  //     }
  //  }
  //}

  TEUCHOS_UNIT_TEST(periodic_search, PeriodicBC_Matcher_multi)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",4);
       pl->set("Y Elements",2);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }
    TEST_ASSERT(mesh!=Teuchos::null);

    {
       panzer_stk::CoordMatcher xmatcher(0);
       panzer_stk::CoordMatcher ymatcher(1);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> tb_Match
             = panzer_stk::buildPeriodicBC_Matcher("top","bottom",xmatcher);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> lr_Match
             = panzer_stk::buildPeriodicBC_Matcher("left","right",ymatcher);

       RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds;
       std::vector<std::string> matchedSides;
       globallyMatchedIds = tb_Match->getMatchedPair(*mesh,matchedSides);
       matchedSides.push_back("top");
       globallyMatchedIds = lr_Match->getMatchedPair(*mesh,matchedSides,globallyMatchedIds);

       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          out << pair.first << " " << pair.second << std::endl;
       }

       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];

          if(pair.first==1 || pair.first==19)
          {   TEST_EQUALITY(pair.second,9); }
          else if(pair.first==10)
          {   TEST_EQUALITY(pair.second,18); }
          else
          {   TEST_EQUALITY(pair.second,pair.first-18); }
       }

    }
  }

  TEUCHOS_UNIT_TEST(periodic_search, PeriodicBC_Matcher_multi_edge)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",4);
       pl->set("Y Elements",2);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }
    TEST_ASSERT(mesh!=Teuchos::null);

       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_left = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"left");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_right = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"right");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_left = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"left","edge");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_right = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"right","edge");
       std::vector<std::size_t> & sideIds_left = *idsAndCoords_left.first;
       std::vector<std::size_t> & sideIds_right = *idsAndCoords_right.first;
       std::vector<std::size_t> & sideIds_edge_left = *idsAndCoords_edge_left.first;
       std::vector<std::size_t> & sideIds_edge_right = *idsAndCoords_edge_right.first;
       std::vector<Tuple<double,3> > & sideCoords_left = *idsAndCoords_left.second;
       std::vector<Tuple<double,3> > & sideCoords_right = *idsAndCoords_right.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_left = *idsAndCoords_edge_left.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_right = *idsAndCoords_edge_right.second;

       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_top = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"top");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_bottom = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"bottom");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_top = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"top","edge");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_bottom = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"bottom","edge");
       std::vector<std::size_t> & sideIds_top = *idsAndCoords_top.first;
       std::vector<std::size_t> & sideIds_bottom = *idsAndCoords_bottom.first;
       std::vector<std::size_t> & sideIds_edge_top = *idsAndCoords_edge_top.first;
       std::vector<std::size_t> & sideIds_edge_bottom = *idsAndCoords_edge_bottom.first;
       std::vector<Tuple<double,3> > & sideCoords_top = *idsAndCoords_top.second;
       std::vector<Tuple<double,3> > & sideCoords_bottom = *idsAndCoords_bottom.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_top = *idsAndCoords_edge_top.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_bottom = *idsAndCoords_edge_bottom.second;

    {
       panzer_stk::CoordMatcher xmatcher(0);
       panzer_stk::CoordMatcher ymatcher(1);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> tb_Match
             = panzer_stk::buildPeriodicBC_Matcher("top","bottom",xmatcher,"edge");
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> lr_Match
             = panzer_stk::buildPeriodicBC_Matcher("left","right",ymatcher,"edge");

       RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds_edge;
       std::vector<std::string> matchedSides;
       globallyMatchedIds_edge = tb_Match->getMatchedPair(*mesh,matchedSides);
       matchedSides.push_back("top");
       globallyMatchedIds_edge = lr_Match->getMatchedPair(*mesh,matchedSides,globallyMatchedIds_edge);

       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds_edge->size();i++) {

          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds_edge)[i];

          // Get coordinates for matched edges on top and bottom
          Tuple<double,3> coord_top;
          Tuple<double,3> coord_bottom;
          int flag_tb = 0;
          for(std::size_t j=0;j<sideIds_edge_top.size();j++){
            if(pair.first == sideIds_edge_top[j]){
              coord_top = sideCoords_edge_top[j];
              flag_tb++;
            }
            if(pair.second == sideIds_edge_bottom[j]){
              coord_bottom = sideCoords_edge_bottom[j];
              flag_tb++;
            }
          }

          // Get coordinates for matched edges on left and right
          Tuple<double,3> coord_left;
          Tuple<double,3> coord_right;
          int flag_lr = 0;
          for(std::size_t j=0;j<sideIds_edge_left.size();j++){
            if(pair.first == sideIds_edge_left[j]){
              coord_left = sideCoords_edge_left[j];
              flag_lr++;
            }
            if(pair.second == sideIds_edge_right[j]){
              coord_right = sideCoords_edge_right[j];
              flag_lr++;
            }
          }
          TEST_EQUALITY(flag_tb+flag_lr,2);

          // If node is on top
          if(flag_tb == 2) {
             int top_l  = -1;
             int top_u  = -1;
             int bottom_l = -1;
             int bottom_u = -1;
             int flag   = 0;

             for(std::size_t j=0;j<sideCoords_top.size();j++) {
               if ((std::abs(sideCoords_top[j][0] - (coord_top[0]-1.0/16.0)) < 1e-14) && (std::abs(sideCoords_top[j][1] - coord_top[1]) < 1e-14)){
                  top_l = j;
                  flag++;
               }
               if ((std::abs(sideCoords_top[j][0] - (coord_top[0]+1.0/16.0)) < 1e-14) && (std::abs(sideCoords_top[j][1] - coord_top[1]) < 1e-14)){
                  top_u = j;
                  flag++;
               }
               if ((std::abs(sideCoords_bottom[j][0] - (coord_bottom[0]-1.0/16.0)) < 1e-14) && (std::abs(sideCoords_bottom[j][1] - coord_bottom[1]) < 1e-14)){
                  bottom_l = j;
                  flag++;
               }
               if ((std::abs(sideCoords_bottom[j][0] - (coord_bottom[0]+1.0/16.0)) < 1e-14) && (std::abs(sideCoords_bottom[j][1] - coord_bottom[1]) < 1e-14)){
                  bottom_u = j;
                  flag++;
               }
             }
             TEST_EQUALITY(flag,4);

             // Test equivalence of node numbers
             TEST_EQUALITY(sideIds_top[top_l],sideIds_bottom[bottom_l]+18);
             TEST_EQUALITY(sideIds_top[top_u],sideIds_bottom[bottom_u]+18);
          }
          // If node is on left
          else if(flag_lr == 2) {
             int left_l  = -1;
             int left_u  = -1;
             int right_l = -1;
             int right_u = -1;
             int flag   = 0;

             for(std::size_t j=0;j<sideCoords_left.size();j++) {
               if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]-1.0/4.0)) < 1e-14)){
                  left_l = j;
                  flag++;
               }
               if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]+1.0/4.0)) < 1e-14)){
                  left_u = j;
                  flag++;
               }
               if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]-1.0/4.0)) < 1e-14)){
                  right_l = j;
                  flag++;
               }
               if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]+1.0/4.0)) < 1e-14)){
                  right_u = j;
                  flag++;
               }
             }
             TEST_EQUALITY(flag,4);

             // Test equivalence of node numbers
             TEST_EQUALITY(sideIds_left[left_l],sideIds_right[right_l]-8);
             TEST_EQUALITY(sideIds_left[left_u],sideIds_right[right_u]-8);
          }

       }

    }
  }

  TEUCHOS_UNIT_TEST(periodic_search, PeriodicBC_Matcher_multi_face)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::CubeHexMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("Z Blocks",1);
       pl->set("X Elements",4);
       pl->set("Y Elements",2);
       pl->set("Z Elements",1);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }
    TEST_ASSERT(mesh!=Teuchos::null);

       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_left = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"left");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_right = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"right");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_face_left = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"left","face");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_face_right = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"right","face");
       std::vector<std::size_t> & sideIds_left = *idsAndCoords_left.first;
       std::vector<std::size_t> & sideIds_right = *idsAndCoords_right.first;
       std::vector<std::size_t> & sideIds_face_left = *idsAndCoords_face_left.first;
       std::vector<std::size_t> & sideIds_face_right = *idsAndCoords_face_right.first;
       std::vector<Tuple<double,3> > & sideCoords_left = *idsAndCoords_left.second;
       std::vector<Tuple<double,3> > & sideCoords_right = *idsAndCoords_right.second;
       std::vector<Tuple<double,3> > & sideCoords_face_left = *idsAndCoords_face_left.second;
       std::vector<Tuple<double,3> > & sideCoords_face_right = *idsAndCoords_face_right.second;

       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_top = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"top");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_bottom = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"bottom");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_face_top = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"top","face");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_face_bottom = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"bottom","face");
       std::vector<std::size_t> & sideIds_top = *idsAndCoords_top.first;
       std::vector<std::size_t> & sideIds_bottom = *idsAndCoords_bottom.first;
       std::vector<std::size_t> & sideIds_face_top = *idsAndCoords_face_top.first;
       std::vector<std::size_t> & sideIds_face_bottom = *idsAndCoords_face_bottom.first;
       std::vector<Tuple<double,3> > & sideCoords_top = *idsAndCoords_top.second;
       std::vector<Tuple<double,3> > & sideCoords_bottom = *idsAndCoords_bottom.second;
       std::vector<Tuple<double,3> > & sideCoords_face_top = *idsAndCoords_face_top.second;
       std::vector<Tuple<double,3> > & sideCoords_face_bottom = *idsAndCoords_face_bottom.second;

    {
       panzer_stk::CoordMatcher xmatcher(0);
       panzer_stk::CoordMatcher ymatcher(1);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> tb_Match
             = panzer_stk::buildPeriodicBC_Matcher("top","bottom",xmatcher,"face");
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> lr_Match
             = panzer_stk::buildPeriodicBC_Matcher("left","right",ymatcher,"face");

       RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds_face;
       std::vector<std::string> matchedSides;
       globallyMatchedIds_face = tb_Match->getMatchedPair(*mesh,matchedSides);
       matchedSides.push_back("top");
       globallyMatchedIds_face = lr_Match->getMatchedPair(*mesh,matchedSides,globallyMatchedIds_face);

       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds_face->size();i++) {

          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds_face)[i];

          // Get coordinates for matched faces on top and bottom
          Tuple<double,3> coord_top;
          Tuple<double,3> coord_bottom;
          int flag_tb = 0;
          for(std::size_t j=0;j<sideIds_face_top.size();j++){
            if(pair.first == sideIds_face_top[j]){
              coord_top = sideCoords_face_top[j];
              flag_tb++;
            }
            if(pair.second == sideIds_face_bottom[j]){
              coord_bottom = sideCoords_face_bottom[j];
              flag_tb++;
            }
          }

          // Get coordinates for matched faces on left and right
          Tuple<double,3> coord_left;
          Tuple<double,3> coord_right;
          int flag_lr = 0;
          for(std::size_t j=0;j<sideIds_face_left.size();j++){
            if(pair.first == sideIds_face_left[j]){
              coord_left = sideCoords_face_left[j];
              flag_lr++;
            }
            if(pair.second == sideIds_face_right[j]){
              coord_right = sideCoords_face_right[j];
              flag_lr++;
            }
          }
          TEST_EQUALITY(flag_tb+flag_lr,2);

          // If node is on top
          if(flag_tb == 2) {
             int top_ll  = -1;
             int top_lu  = -1;
             int top_ul  = -1;
             int top_uu  = -1;
             int bottom_ll = -1;
             int bottom_lu = -1;
             int bottom_ul = -1;
             int bottom_uu = -1;
             int flag   = 0;

             for(std::size_t j=0;j<sideCoords_top.size();j++) {
               if ((std::abs(sideCoords_top[j][0] - (coord_top[0]-1.0/16.0)) < 1e-14) && (std::abs(sideCoords_top[j][1] - coord_top[1]) < 1e-14) && (std::abs(sideCoords_top[j][2] - (coord_top[2]-1.0/2.0)) < 1e-14)){
                  top_ll = j;
                  flag++;
               }
               if ((std::abs(sideCoords_top[j][0] - (coord_top[0]-1.0/16.0)) < 1e-14) && (std::abs(sideCoords_top[j][1] - coord_top[1]) < 1e-14) && (std::abs(sideCoords_top[j][2] - (coord_top[2]+1.0/2.0)) < 1e-14)){
                  top_lu = j;
                  flag++;
               }
               if ((std::abs(sideCoords_top[j][0] - (coord_top[0]+1.0/16.0)) < 1e-14) && (std::abs(sideCoords_top[j][1] - coord_top[1]) < 1e-14) && (std::abs(sideCoords_top[j][2] - (coord_top[2]-1.0/2.0)) < 1e-14)){
                  top_ul = j;
                  flag++;
               }
               if ((std::abs(sideCoords_top[j][0] - (coord_top[0]+1.0/16.0)) < 1e-14) && (std::abs(sideCoords_top[j][1] - coord_top[1]) < 1e-14) && (std::abs(sideCoords_top[j][2] - (coord_top[2]+1.0/2.0)) < 1e-14)){
                  top_uu = j;
                  flag++;
               }
               if ((std::abs(sideCoords_bottom[j][0] - (coord_bottom[0]-1.0/16.0)) < 1e-14) && (std::abs(sideCoords_bottom[j][1] - coord_bottom[1]) < 1e-14) && (std::abs(sideCoords_bottom[j][2] - (coord_bottom[2]-1.0/2.0)) < 1e-14)){
                  bottom_ll = j;
                  flag++;
               }
               if ((std::abs(sideCoords_bottom[j][0] - (coord_bottom[0]-1.0/16.0)) < 1e-14) && (std::abs(sideCoords_bottom[j][1] - coord_bottom[1]) < 1e-14) && (std::abs(sideCoords_bottom[j][2] - (coord_bottom[2]+1.0/2.0)) < 1e-14)){
                  bottom_lu = j;
                  flag++;
               }
               if ((std::abs(sideCoords_bottom[j][0] - (coord_bottom[0]+1.0/16.0)) < 1e-14) && (std::abs(sideCoords_bottom[j][1] - coord_bottom[1]) < 1e-14) && (std::abs(sideCoords_bottom[j][2] - (coord_bottom[2]-1.0/2.0)) < 1e-14)){
                  bottom_ul = j;
                  flag++;
               }
               if ((std::abs(sideCoords_bottom[j][0] - (coord_bottom[0]+1.0/16.0)) < 1e-14) && (std::abs(sideCoords_bottom[j][1] - coord_bottom[1]) < 1e-14) && (std::abs(sideCoords_bottom[j][2] - (coord_bottom[2]+1.0/2.0)) < 1e-14)){
                  bottom_uu = j;
                  flag++;
               }
             }
             TEST_EQUALITY(flag,8);

             // Test equivalence of node numbers
             TEST_EQUALITY(sideIds_top[top_ll],sideIds_bottom[bottom_ll]+18);
             TEST_EQUALITY(sideIds_top[top_lu],sideIds_bottom[bottom_lu]+18);
             TEST_EQUALITY(sideIds_top[top_ul],sideIds_bottom[bottom_ul]+18);
             TEST_EQUALITY(sideIds_top[top_uu],sideIds_bottom[bottom_uu]+18);
          }
          // If node is on left
          else if(flag_lr == 2) {
             int left_ll  = -1;
             int left_lu  = -1;
             int left_ul  = -1;
             int left_uu  = -1;
             int right_ll = -1;
             int right_lu = -1;
             int right_ul = -1;
             int right_uu = -1;
             int flag   = 0;

             for(std::size_t j=0;j<sideCoords_left.size();j++) {
               if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]-1.0/4.0)) < 1e-14) && (std::abs(sideCoords_left[j][2] - (coord_left[2]-1.0/2.0)) < 1e-14)){
                  left_ll = j;
                  flag++;
               }
               if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]-1.0/4.0)) < 1e-14) && (std::abs(sideCoords_left[j][2] - (coord_left[2]+1.0/2.0)) < 1e-14)){
                  left_lu = j;
                  flag++;
               }
               if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]+1.0/4.0)) < 1e-14) && (std::abs(sideCoords_left[j][2] - (coord_left[2]-1.0/2.0)) < 1e-14)){
                  left_ul = j;
                  flag++;
               }
               if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]+1.0/4.0)) < 1e-14) && (std::abs(sideCoords_left[j][2] - (coord_left[2]+1.0/2.0)) < 1e-14)){
                  left_uu = j;
                  flag++;
               }
               if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]-1.0/4.0)) < 1e-14) && (std::abs(sideCoords_right[j][2] - (coord_right[2]-1.0/2.0)) < 1e-14)){
                  right_ll = j;
                  flag++;
               }
               if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]-1.0/4.0)) < 1e-14) && (std::abs(sideCoords_right[j][2] - (coord_right[2]+1.0/2.0)) < 1e-14)){
                  right_lu = j;
                  flag++;
               }
               if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]+1.0/4.0)) < 1e-14) && (std::abs(sideCoords_right[j][2] - (coord_right[2]-1.0/2.0)) < 1e-14)){
                  right_ul = j;
                  flag++;
               }
               if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]+1.0/4.0)) < 1e-14) && (std::abs(sideCoords_right[j][2] - (coord_right[2]+1.0/2.0)) < 1e-14)){
                  right_uu = j;
                  flag++;
               }
             }
             TEST_EQUALITY(flag,8);

             // Test equivalence of node numbers
             TEST_EQUALITY(sideIds_left[left_ll],sideIds_right[right_ll]-8);
             TEST_EQUALITY(sideIds_left[left_lu],sideIds_right[right_lu]-8);
             TEST_EQUALITY(sideIds_left[left_ul],sideIds_right[right_ul]-8);
             TEST_EQUALITY(sideIds_left[left_uu],sideIds_right[right_uu]-8);
          }

       }

    }
  }

  TEUCHOS_UNIT_TEST(periodic_search, PeriodicBC_Matcher_nodes_and_edges)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",4);
       pl->set("Y Elements",2);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }
    TEST_ASSERT(mesh!=Teuchos::null);

       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_left = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"left");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_right = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"right");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_left = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"left","edge");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_right = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"right","edge");
       std::vector<std::size_t> & sideIds_left = *idsAndCoords_left.first;
       std::vector<std::size_t> & sideIds_right = *idsAndCoords_right.first;
       std::vector<std::size_t> & sideIds_edge_left = *idsAndCoords_edge_left.first;
       std::vector<std::size_t> & sideIds_edge_right = *idsAndCoords_edge_right.first;
       std::vector<Tuple<double,3> > & sideCoords_left = *idsAndCoords_left.second;
       std::vector<Tuple<double,3> > & sideCoords_right = *idsAndCoords_right.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_left = *idsAndCoords_edge_left.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_right = *idsAndCoords_edge_right.second;


    {
       panzer_stk::CoordMatcher ymatcher(1);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> node_Match
             = panzer_stk::buildPeriodicBC_Matcher("left","right",ymatcher);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> edge_Match
             = panzer_stk::buildPeriodicBC_Matcher("left","right",ymatcher,"edge");

       RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds;
       std::vector<std::vector<std::string> > matchedSides(3);
       // matchedSides is a vector of vectors here (nodes,edges,faces)
       globallyMatchedIds = node_Match->getMatchedPair(*mesh,matchedSides[0]);
       matchedSides[0].push_back("left");
       globallyMatchedIds = edge_Match->getMatchedPair(*mesh,matchedSides[1],globallyMatchedIds);
       matchedSides[1].push_back("left");

       // match left & right sides
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {

          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];

          // Is this a node or edge pairing?
          if(pair.first < 27){ // Node
            TEST_EQUALITY(pair.first,pair.second-8);
          } else {             // Edge

            // Get coordinates for matched edges on left and right
            Tuple<double,3> coord_left;
            Tuple<double,3> coord_right;
            int flag_lr = 0;
            for(std::size_t j=0;j<sideIds_edge_left.size();j++){
              if(pair.first == sideIds_edge_left[j]){
                coord_left = sideCoords_edge_left[j];
                flag_lr++;
              }
              if(pair.second == sideIds_edge_right[j]){
                coord_right = sideCoords_edge_right[j];
                flag_lr++;
              }
            }
            TEST_EQUALITY(flag_lr,2);

            int left_l  = -1;
            int left_u  = -1;
            int right_l = -1;
            int right_u = -1;
            int flag   = 0;

            for(std::size_t j=0;j<sideCoords_left.size();j++) {
              if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]-1.0/4.0)) < 1e-14)){
                 left_l = j;
                 flag++;
              }
              if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]+1.0/4.0)) < 1e-14)){
                 left_u = j;
                 flag++;
              }
              if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]-1.0/4.0)) < 1e-14)){
                 right_l = j;
                 flag++;
              }
              if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]+1.0/4.0)) < 1e-14)){
                 right_u = j;
                 flag++;
              }
            }
            TEST_EQUALITY(flag,4);
            // Test equivalence of node numbers
            TEST_EQUALITY(sideIds_left[left_l],sideIds_right[right_l]-8);
            TEST_EQUALITY(sideIds_left[left_u],sideIds_right[right_u]-8);
          }
       }

    }

  }

  TEUCHOS_UNIT_TEST(periodic_search, PeriodicBC_Matcher_nodes_edges_and_faces)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::CubeHexMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("Z Blocks",1);
       pl->set("X Elements",4);
       pl->set("Y Elements",2);
       pl->set("Z Elements",1);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }
    TEST_ASSERT(mesh!=Teuchos::null);

       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_left = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"left");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_right = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"right");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_left = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"left","edge");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_right = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"right","edge");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_face_left = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"left","face");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_face_right = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"right","face");
       std::vector<std::size_t> & sideIds_left = *idsAndCoords_left.first;
       std::vector<std::size_t> & sideIds_right = *idsAndCoords_right.first;
       std::vector<std::size_t> & sideIds_edge_left = *idsAndCoords_edge_left.first;
       std::vector<std::size_t> & sideIds_edge_right = *idsAndCoords_edge_right.first;
       std::vector<std::size_t> & sideIds_face_left = *idsAndCoords_face_left.first;
       std::vector<std::size_t> & sideIds_face_right = *idsAndCoords_face_right.first;
       std::vector<Tuple<double,3> > & sideCoords_left = *idsAndCoords_left.second;
       std::vector<Tuple<double,3> > & sideCoords_right = *idsAndCoords_right.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_left = *idsAndCoords_edge_left.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_right = *idsAndCoords_edge_right.second;
       std::vector<Tuple<double,3> > & sideCoords_face_left = *idsAndCoords_face_left.second;
       std::vector<Tuple<double,3> > & sideCoords_face_right = *idsAndCoords_face_right.second;

    {
       panzer_stk::PlaneMatcher ymatcher(1,2);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> node_Match
             = panzer_stk::buildPeriodicBC_Matcher("left","right",ymatcher);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> edge_Match
             = panzer_stk::buildPeriodicBC_Matcher("left","right",ymatcher,"edge");
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> face_Match
             = panzer_stk::buildPeriodicBC_Matcher("left","right",ymatcher,"face");

       RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds;
       std::vector<std::vector<std::string> > matchedSides(3);
       // matchedSides is a vector of vectors here (nodes,edges,faces)
       globallyMatchedIds = node_Match->getMatchedPair(*mesh,matchedSides[0]);
       matchedSides[0].push_back("left");
       globallyMatchedIds = edge_Match->getMatchedPair(*mesh,matchedSides[1],globallyMatchedIds);
       matchedSides[1].push_back("left");
       globallyMatchedIds = face_Match->getMatchedPair(*mesh,matchedSides[2],globallyMatchedIds);
       matchedSides[2].push_back("left");

       // match left & right sides
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {

          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          // Is this a node, edge, or face pairing?
          if(pair.first < 54){ // Node
            TEST_EQUALITY(pair.first,pair.second-8);
          } else if(pair.first < 192){ //Edge

            // Get coordinates for matched edges on left and right
            Tuple<double,3> coord_left;
            Tuple<double,3> coord_right;
            int flag_lr = 0;
            for(std::size_t j=0;j<sideIds_edge_left.size();j++){
              if(pair.first == sideIds_edge_left[j]){
                coord_left = sideCoords_edge_left[j];
                flag_lr++;
              }
              if(pair.second == sideIds_edge_right[j]){
                coord_right = sideCoords_edge_right[j];
                flag_lr++;
              }
            }
            TEST_EQUALITY(flag_lr,2);

            int left_l  = -1;
            int left_u  = -1;
            int right_l = -1;
            int right_u = -1;
            int flag   = 0;

            for(std::size_t j=0;j<sideCoords_left.size();j++) {
              if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]-1.0/4.0)) < 1e-14) && (std::abs(sideCoords_left[j][2] - coord_left[2]) < 1e-14)){
                 left_l = j;
                 flag++;
              }
              if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - coord_left[1]) < 1e-14) && (std::abs(sideCoords_left[j][2] - (coord_left[2]-1.0/2.0)) < 1e-14)){
                 left_l = j;
                 flag++;
              }
              if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]+1.0/4.0)) < 1e-14) && (std::abs(sideCoords_left[j][2] - coord_left[2]) < 1e-14)){
                 left_u = j;
                 flag++;
              }
              if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - coord_left[1]) < 1e-14) && (std::abs(sideCoords_left[j][2] - (coord_left[2]+1.0/2.0)) < 1e-14)){
                 left_u = j;
                 flag++;
              }
              if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]-1.0/4.0)) < 1e-14) && (std::abs(sideCoords_right[j][2] - coord_right[2]) < 1e-14)){
                 right_l = j;
                 flag++;
              }
              if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - coord_right[1]) < 1e-14) && (std::abs(sideCoords_right[j][2] - (coord_right[2]-1.0/2.0)) < 1e-14)){
                 right_l = j;
                 flag++;
              }
              if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]+1.0/4.0)) < 1e-14) && (std::abs(sideCoords_right[j][2] - coord_right[2]) < 1e-14)){
                 right_u = j;
                 flag++;
              }
              if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - coord_right[1]) < 1e-14) && (std::abs(sideCoords_right[j][2] - (coord_right[2]+1.0/2.0)) < 1e-14)){
                 right_u = j;
                 flag++;
              }
            }
            TEST_EQUALITY(flag,4);
            // Test equivalence of node numbers
            TEST_EQUALITY(sideIds_left[left_l],sideIds_right[right_l]-8);
            TEST_EQUALITY(sideIds_left[left_u],sideIds_right[right_u]-8);
          } else { //Face

            // Get coordinates for matched faces on left and right
            Tuple<double,3> coord_left;
            Tuple<double,3> coord_right;
            int flag_lr = 0;
            for(std::size_t j=0;j<sideIds_face_left.size();j++){
              if(pair.first == sideIds_face_left[j]){
                coord_left = sideCoords_face_left[j];
                flag_lr++;
              }
              if(pair.second == sideIds_face_right[j]){
                coord_right = sideCoords_face_right[j];
                flag_lr++;
              }
            }
            TEST_EQUALITY(flag_lr,2);

             int left_ll  = -1;
             int left_lu  = -1;
             int left_ul  = -1;
             int left_uu  = -1;
             int right_ll = -1;
             int right_lu = -1;
             int right_ul = -1;
             int right_uu = -1;
             int flag   = 0;

             for(std::size_t j=0;j<sideCoords_left.size();j++) {
               if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]-1.0/4.0)) < 1e-14) && (std::abs(sideCoords_left[j][2] - (coord_left[2]-1.0/2.0)) < 1e-14)){
                  left_ll = j;
                  flag++;
               }
               if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]-1.0/4.0)) < 1e-14) && (std::abs(sideCoords_left[j][2] - (coord_left[2]+1.0/2.0)) < 1e-14)){
                  left_lu = j;
                  flag++;
               }
               if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]+1.0/4.0)) < 1e-14) && (std::abs(sideCoords_left[j][2] - (coord_left[2]-1.0/2.0)) < 1e-14)){
                  left_ul = j;
                  flag++;
               }
               if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]+1.0/4.0)) < 1e-14) && (std::abs(sideCoords_left[j][2] - (coord_left[2]+1.0/2.0)) < 1e-14)){
                  left_uu = j;
                  flag++;
               }
               if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]-1.0/4.0)) < 1e-14) && (std::abs(sideCoords_right[j][2] - (coord_right[2]-1.0/2.0)) < 1e-14)){
                  right_ll = j;
                  flag++;
               }
               if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]-1.0/4.0)) < 1e-14) && (std::abs(sideCoords_right[j][2] - (coord_right[2]+1.0/2.0)) < 1e-14)){
                  right_lu = j;
                  flag++;
               }
               if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]+1.0/4.0)) < 1e-14) && (std::abs(sideCoords_right[j][2] - (coord_right[2]-1.0/2.0)) < 1e-14)){
                  right_ul = j;
                  flag++;
               }
               if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]+1.0/4.0)) < 1e-14) && (std::abs(sideCoords_right[j][2] - (coord_right[2]+1.0/2.0)) < 1e-14)){
                  right_uu = j;
                  flag++;
               }
             }
             TEST_EQUALITY(flag,8);

             // Test equivalence of node numbers
             TEST_EQUALITY(sideIds_left[left_ll],sideIds_right[right_ll]-8);
             TEST_EQUALITY(sideIds_left[left_lu],sideIds_right[right_lu]-8);
             TEST_EQUALITY(sideIds_left[left_ul],sideIds_right[right_ul]-8);
             TEST_EQUALITY(sideIds_left[left_uu],sideIds_right[right_uu]-8);
          }
       }

    }

  }

  TEUCHOS_UNIT_TEST(periodic_search, MatcherQueries)
  {
    using namespace panzer_stk;
    using namespace Teuchos;

    // Box mesh in 2D (match one coordinate on opposite side)
    CoordMatcher cm_2d_periodic_in_x_direction(1);
    TEST_EQUALITY(cm_2d_periodic_in_x_direction.getPeriodicDirection(),0);
    CoordMatcher cm_2d_periodic_in_y_direction(0);
    TEST_EQUALITY(cm_2d_periodic_in_y_direction.getPeriodicDirection(),1);
    CoordMatcher cm_2d_periodic_in_z_direction(2);
    TEST_THROW(cm_2d_periodic_in_z_direction.getPeriodicDirection(),std::runtime_error);

    // Box mesh in 3D (match planes on opposite sides)
    RCP<PlaneMatcher> pm;
    pm = rcp(new PlaneMatcher(1,2));
    TEST_EQUALITY(pm->getPeriodicDirection(),0);
    pm = rcp(new PlaneMatcher(2,1));
    TEST_EQUALITY(pm->getPeriodicDirection(),0);
    pm = rcp(new PlaneMatcher(0,2));
    TEST_EQUALITY(pm->getPeriodicDirection(),1);
    pm = rcp(new PlaneMatcher(2,0));
    TEST_EQUALITY(pm->getPeriodicDirection(),1);
    pm = rcp(new PlaneMatcher(0,1));
    TEST_EQUALITY(pm->getPeriodicDirection(),2);
    pm = rcp(new PlaneMatcher(1,0));
    TEST_EQUALITY(pm->getPeriodicDirection(),2);

    // Wedge in 2D
    RCP<WedgeMatcher> wm;
    std::vector<std::string> params;
    params.push_back("1.0e-8");
    params.push_back("2D");
    wm = rcp(new WedgeMatcher(WedgeMatcher::MirrorPlane::XZ_PLANE,params));
    TEST_ASSERT(wm->getMirrorPlane() == WedgeMatcher::MirrorPlane::XZ_PLANE);
    TEST_ASSERT(!wm->isThreeD());
    TEST_EQUALITY(wm->getIndex(),1);
    wm = rcp(new WedgeMatcher(WedgeMatcher::MirrorPlane::YZ_PLANE,params));
    TEST_ASSERT(wm->getMirrorPlane() == WedgeMatcher::MirrorPlane::YZ_PLANE);
    TEST_ASSERT(!wm->isThreeD());
    TEST_EQUALITY(wm->getIndex(),0);

    // Wedge in 3D
    params[1] = "3D";
    wm = rcp(new WedgeMatcher(WedgeMatcher::MirrorPlane::XZ_PLANE,params));
    TEST_ASSERT(wm->getMirrorPlane() == WedgeMatcher::MirrorPlane::XZ_PLANE);
    TEST_ASSERT(wm->isThreeD());
    TEST_EQUALITY(wm->getIndex(),1);
    wm = rcp(new WedgeMatcher(WedgeMatcher::MirrorPlane::YZ_PLANE,params));
    TEST_ASSERT(wm->getMirrorPlane() == WedgeMatcher::MirrorPlane::YZ_PLANE);
    TEST_ASSERT(wm->isThreeD());
    TEST_EQUALITY(wm->getIndex(),0);
  }
}
#endif
