#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_REFINEMENTFIXTURE_TET_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_REFINEMENTFIXTURE_TET_HPP_

#include <Akri_Unit_RefinementFixture.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_SecondOrderMeshSpecs.hpp>

namespace krino {

class RegularTetRefinement : public RefinementFixture<RegularTet>
{
public:
  RegularTetRefinement()
  {
    set_valid_proc_sizes_for_test({1});
    StkMeshTetFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
  }
  stk::mesh::Entity get_element() { return get_owned_elements()[0]; }
  Edge get_edge(unsigned edgeOrdinal) { std::vector<Edge> elemEdges; fill_entity_edges(mMesh, get_element(), elemEdges); return elemEdges[edgeOrdinal]; }
};

class RegularTet10Refinement : public RefinementFixture<RegularTet10>
{
public:
  RegularTet10Refinement()
  {
    set_valid_proc_sizes_for_test({1});
    StkMeshTet10Fixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
  }
  stk::mesh::Entity get_element() { return get_owned_elements()[0]; }
  Edge get_edge(unsigned edgeOrdinal) { std::vector<Edge> elemEdges; fill_entity_edges(mMesh, get_element(), elemEdges); return elemEdges[edgeOrdinal]; }
};


class UMRRegularTetRefinement : public RefinementFixture<UMRRegularTet>
{
public:
  UMRRegularTetRefinement()
  {
    set_valid_proc_sizes_for_test({1,2,4,8});
    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1,1,1,1,1}, {0,1,1,0,0,0,1,1}); // Balanced for refinement along x=0
    else if(stk::parallel_machine_size(mComm) == 4)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1,1,1,1,1}, {0,1,2,3,0,0,1,1}); // Balanced for refinement along x=0
    else if(stk::parallel_machine_size(mComm) == 8)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1,1,1,1,1}, {0,1,2,3,4,5,6,7});
  }
protected:
};

class UMRRegularTet10Refinement : public RefinementFixture<UMRRegularTet10>
{
public:
  UMRRegularTet10Refinement()
  {
    set_valid_proc_sizes_for_test({1,2,4,8});
    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1,1,1,1,1}, {0,1,1,0,0,0,1,1}); // Balanced for refinement along x=0
    else if(stk::parallel_machine_size(mComm) == 4)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1,1,1,1,1}, {0,1,2,3,0,0,1,1}); // Balanced for refinement along x=0
    else if(stk::parallel_machine_size(mComm) == 8)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1,1,1,1,1}, {0,1,2,3,4,5,6,7});
  }
protected:
};

class FourRightTetsRefinement : public RefinementFixture<FourRightTets>
{
public:
  FourRightTetsRefinement()
  {
    set_valid_proc_sizes_for_test({1,2,3,4});
    if(stk::parallel_machine_size(mComm) <= 4)
    {
      if(stk::parallel_machine_size(mComm) == 1)
        this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
      else if(stk::parallel_machine_size(mComm) == 2)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, {0,0,1,1});
      else if(stk::parallel_machine_size(mComm) == 3)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, {0,1,2,2});
      else if(stk::parallel_machine_size(mComm) == 4)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, {0,1,2,3});
    }
  }
protected:
};

class RightTetSurroundedByEdgeTetsRefinement : public RefinementFixture<RightTetSurroundedByEdgeTets>
{
public:
  RightTetSurroundedByEdgeTetsRefinement()
  {
    set_valid_proc_sizes_for_test({1,2,3,4});
    if(stk::parallel_machine_size(mComm) <= 4)
    {
      if(stk::parallel_machine_size(mComm) == 1)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2,2,2,2,2,2}, {0,0,0,0,0,0,0});
      else if(stk::parallel_machine_size(mComm) == 2)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2,2,2,2,2,2}, {0,1,0,1,0,1,0});
      else if(stk::parallel_machine_size(mComm) == 3)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2,2,2,2,2,2}, {0,1,2,0,1,2,0});
      else if(stk::parallel_machine_size(mComm) == 4)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2,2,2,2,2,2}, {0,1,2,3,0,1,2});
    }
  }
protected:
};

class RightTetSurroundedByFaceTetsRefinement : public RefinementFixture<RightTetSurroundedByFaceTets>
{
public:
  RightTetSurroundedByFaceTetsRefinement()
  {
    set_valid_proc_sizes_for_test({1,2,3,4});
    if(stk::parallel_machine_size(mComm) <= 4)
    {
      if(stk::parallel_machine_size(mComm) == 1)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2,2,2,2}, {0,0,0,0,0});
      else if(stk::parallel_machine_size(mComm) == 2)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2,2,2,2}, {0,1,0,1,0});
      else if(stk::parallel_machine_size(mComm) == 3)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2,2,2,2}, {0,1,2,0,1});
      else if(stk::parallel_machine_size(mComm) == 4)
        this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2,2,2,2}, {0,1,2,3,0});
    }
  }
protected:
};

}

#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_REFINEMENTFIXTURE_TET_HPP_ */
