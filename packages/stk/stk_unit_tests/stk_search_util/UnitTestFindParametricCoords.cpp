/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <gtest/gtest.h>
#include <stk_util/stk_config.h>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_search_util/MasterElementProviderIntrepid2.hpp>
#include <stk_search_util/Intrepid2_HasParametricDistance.hpp>
#include <stk_unit_test_utils/TransferTesters.hpp>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

std::shared_ptr<stk::search::MasterElementProviderInterface>
create_MasterElementProvider(bool useCompositeTet10 = false)
{
  return std::make_shared<stk::search::MasterElementProviderIntrepid2>(useCompositeTet10);
}

TEST(TransferUtil, hex8_findParamCoords)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8"
        "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1";

  std::vector<double> pointInside = {0.5, 0.5, 0.5};
  std::vector<double> pointOutside = {1.5, 1.5, 1.5};

  double parametricDistanceInside = 0.0;
  double parametricDistanceOutside = 2.0;

  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointInside, parametricDistanceInside);
  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOutside, parametricDistanceOutside);
}

TEST(TransferUtil, tet4_findParamCoords)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::string meshDesc = "0,1,TET_4,1,2,3,4"
        "|coordinates:   0,0,0, 1,0,0, 0,1,0, 0,0,1";

  std::vector<double> pointInside = {0.25, 0.25, 0.25};
  std::vector<double> pointOutside = {1.25, 1.25, 1.25};

  double parametricDistanceInside = 0.0;
  double parametricDistanceOutside = 12.0;

  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointInside, parametricDistanceInside);
  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOutside, parametricDistanceOutside);
}

TEST(TransferUtil, tet10_quadratic_findParamCoords)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::string meshDesc = "0,1,TET_10,1,2,3,4,5,6,7,8,9,10"
        "|coordinates:   0,0,0, 1,0,0, 0,1,0, 0,0,1,"
                        " 0.5,0,0, 0.5,0.5,0, 0,0.5,0,"
                        " 0,0,0.5, 0.5,0,0.5, 0,0.5,0.5";

  std::vector<double> pointInside = {0.25, 0.25, 0.25};
  std::vector<double> pointOnEdge = {0.0, 0.0, 0.5};
  std::vector<double> pointOutside = {1.25, 1.25, 1.25};

  double parametricDistanceInside = 0.0;
  double parametricDistanceOnEdge = 1.0;
  double parametricDistanceOutside = 12.0;

  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointInside, parametricDistanceInside);
  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOnEdge, parametricDistanceOnEdge);
  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOutside, parametricDistanceOutside);
}

TEST(TransferUtil, tet10_composite_findParamCoords)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::string meshDesc = "0,1,TET_10,1,2,3,4,5,6,7,8,9,10"
        "|coordinates:   0,0,0, 1,0,0, 0,1,0, 0,0,1,"
                        " 0.5,0,0, 0.5,0.5,0, 0,0.5,0,"
                        " 0,0,0.5, 0.5,0,0.5, 0,0.5,0.5";

  std::vector<double> pointInside = {0.25, 0.25, 0.25};
  std::vector<double> pointOnEdge = {0.0, 0.0, 0.5};
  std::vector<double> pointOutside = {1.25, 1.25, 1.25};

  double parametricDistanceInside = 0.0;
  double parametricDistanceOnEdge = 1.0;
  double parametricDistanceOutside = 222.0;

  const bool useCompositeTet10 = true;
  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(useCompositeTet10), pointInside, parametricDistanceInside);
  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(useCompositeTet10), pointOnEdge, parametricDistanceOnEdge);
  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(useCompositeTet10), pointOutside, parametricDistanceOutside);
}

TEST(TransferUtil, deformed_tet10_findParamCoords)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::string meshDesc = "0,1,TET_10,1,2,3,4,5,6,7,8,9,10"
        "|coordinates:   0,0,0, 1,0,0, 0,1,0, 0,0,1,"
                        "0.5,0,0, 0.5,0.5,0, 0,0.5,0,"
                        "0.125,0.0,0.5, 0.5,0,0.5, 0,0.5,0.5";

  std::vector<double> pointAtMidEdgeNode = {0.125, 0.0, 0.5};

  std::vector<double> pointOnEdgeBetweenNodesComposite = {0.111, 0.0, 0.25}; //?????
  //I don't think this point is actually on the edge... But this is the point
  //on the z=0.25 line, where parametric-distance is 1.0
  
  std::vector<double> pointOnEdgeBetweenNodesQuadratic = {0.095, 0.0, 0.25};

  double parametricDistanceAtMidEdgeNodeComposite = 1.0;
  double parametricDistanceOnEdgeComposite = 1.0;
  double parametricDistanceAtMidEdgeNodeQuadratic = 1.0;
  double parametricDistanceOnEdgeQuadratic = 1.0;

  bool useCompositeTet10 = true;
  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(useCompositeTet10), pointAtMidEdgeNode, parametricDistanceAtMidEdgeNodeComposite);
  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(useCompositeTet10), pointOnEdgeBetweenNodesComposite, parametricDistanceOnEdgeComposite);
  useCompositeTet10 = false;
  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(useCompositeTet10), pointAtMidEdgeNode, parametricDistanceAtMidEdgeNodeQuadratic);
  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(useCompositeTet10), pointOnEdgeBetweenNodesQuadratic, parametricDistanceOnEdgeQuadratic);
}

TEST(TransferUtil, pyramid5_findParamCoords)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::string meshDesc = "0,1,PYRAMID_5,1,2,3,4,5"
        "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,1,0, 0.5,0.5,1";

  std::vector<double> pointInside = {0.5, 0.5, 0.2};
  std::vector<double> pointOutside = {1.0, 1.5, 1.0};

  double parametricDistanceInside = 0.4;
  double parametricDistanceOutside = 4.0;

  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointInside, parametricDistanceInside);
  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOutside, parametricDistanceOutside);
}

TEST(TransferUtil, wedge6_findParamCoords)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::string meshDesc = "0,1,WEDGE_6,1,2,3,4,5,6"
        "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,1,1, 1,1,1";

  std::vector<double> pointInside = {0.5, 2.0 / 3.0, 1.0 / 3.0};
  std::vector<double> pointOutside = {1.5, 1.0, 1.0};

  double parametricDistanceInside = 0.0;
  double parametricDistanceOutside = 4.0;

  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointInside, parametricDistanceInside);
  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOutside, parametricDistanceOutside);
}

TEST(TransferUtil, beam2_findParamCoords)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::string meshDesc = "0,1,BEAM_2,1,2|dimension:3"
        "|coordinates:  0,0,0, 1,0,0";

  std::vector<double> pointOnLine = {0.5, 0.0, 0.0};
  std::vector<double> pointOutside = {1.2, 0.0, 0.0};

  double parametricDistanceInside = 0.0;
  double parametricDistanceOutside = 1.4;

  EXPECT_ANY_THROW(stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOnLine, parametricDistanceInside));
  EXPECT_ANY_THROW(stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOutside, parametricDistanceOutside));

}

TEST(TransferUtil, quad4_2d_findParamCoords)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::string meshDesc =
      "0,1,QUAD_4_2D,1,2,3,4,block_1\n"
      "|coordinates:   0,0, 1,0, 1,1, 0,1"
      "|dimension:2";

  std::vector<double> pointOnQuad = {0.5, 0.5};
  std::vector<double> pointOutside = {1.2, 0.0};

  double parametricDistanceInside = 0.0;
  double parametricDistanceOutside = 1.4;

  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOnQuad, parametricDistanceInside);
  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOutside, parametricDistanceOutside);
}

TEST(TransferUtil, tri3_2d_findParamCoords)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  std::string meshDesc = "0,1,TRI_3_2D,1,2,3"
      "|coordinates:   1,0, 0,1, 0,0"
      "|dimension:2";

  std::vector<double> pointOnTri = {1/3.0, 1/3.0};

  const double delta = 0.01;
  std::vector<double> pointOutside = {1.0 + delta, 0.0 + delta};

  double parametricDistanceInside = 0.0;
  double parametricDistanceOutside = 1.06;

  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOnTri, parametricDistanceInside);
  stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOutside, parametricDistanceOutside);
}

TEST(TransferUtil, shell_quad4_findParamCoords)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  if constexpr (stk::search::intrepid2_has_param_dist) {
    std::string meshDesc = "0,1,SHELL_QUAD_4,1,2,3,4"
          "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,1,0";

    std::vector<double> pointOnShell = {0.5, 0.5, 0.0};
    std::vector<double> pointOutside = {0.5, 1.2, 0.0};

    double parametricDistanceInside = 0.0;
    double parametricDistanceOutside = 1.4;


    stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOnShell, parametricDistanceInside);
    stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOutside, parametricDistanceOutside);
  }
  else {
    std::cout<<"skipping because Intrepid2 doesn't have parametric-distance"<<std::endl;
  }
}

TEST(TransferUtil, shell_tri3_findParamCoords)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) { GTEST_SKIP(); }

  if constexpr (stk::search::intrepid2_has_param_dist) {
    std::string meshDesc = "0,1,SHELL_TRI_3,1,2,3"
          "|coordinates:   1,0,0, 0,1,0, 0,0,1";
  
    // Shell triangle: x + y + z = 1 with unit normal [1,1,1]
    const double delta = 0.001;
    std::vector<double> pointOnShell  = {1/3.0, 1/3.0, 1/3.0};
    std::vector<double> pointOutside  = {1/3.0 + delta, 1/3.0 + delta, 1/3.0 + delta};
    std::vector<double> pointOutside2 = {1.0 + 2*delta, 0.0 + delta, 0.0 + delta};
  
    double parametricDistanceInside   = 0.0;
    double parametricDistanceOutside  = 0.0;
    double parametricDistanceOutside2 = 1.0+delta;
  
    stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOnShell, parametricDistanceInside);
    stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOutside, parametricDistanceOutside);
    stk::unit_test_util::test_find_parametric_coordinates(meshDesc, stk::topology::ELEM_RANK, 1, create_MasterElementProvider(), pointOutside2, parametricDistanceOutside2);
  }
  else {
    std::cout<<"skipping because Intrepid2 doesn't have parametric-distance"<<std::endl;
  }
}

} // anonymous namespace
