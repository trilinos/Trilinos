//
//  Intrepid2_ESEAS_Interface.hpp
//  QuadraturePerformance
//
//  Created by Roberts, Nathan V on 9/3/19.
//

#ifndef Intrepid2_ESEAS_Interface_h
#define Intrepid2_ESEAS_Interface_h

#include "Intrepid2_TestUtils.hpp"

static const int N = 8;
static const int VALUE_LENGTH_1D = N+1;
static const int VALUE_LENGTH_TRI_H1_BASIS    = (N+1)*(N+2)/2;
static const int VALUE_LENGTH_TRI_HCURL_BASIS = N*3 + N*(N-1); // N*3 edge functions, plus two families of face functions, each with N*(N-1)/2 members.
static const int VALUE_LENGTH_TRI_HDIV_BASIS  = VALUE_LENGTH_TRI_HCURL_BASIS; // 2D: H(div) is rotated H(curl)
static const int VALUE_LENGTH_TRI_HVOL_BASIS  = N*(N+1)/2;
static const int VALUE_LENGTH_TET_H1_BASIS    = (N+1)*(N+2)*(N+3)/6;
static const int VALUE_LENGTH_WEDGE_H1_BASIS  = VALUE_LENGTH_TRI_H1_BASIS * VALUE_LENGTH_1D;
static const int VALUE_LENGTH_PYRA_H1_BASIS   = N*N*N + 3*N + 1;

static const int VALUE_LENGTH_WEDGE_HCURL_FAMILY_I_BASIS = VALUE_LENGTH_TRI_HCURL_BASIS * VALUE_LENGTH_1D;
static const int VALUE_LENGTH_WEDGE_HDIV_FAMILY_I_BASIS  = VALUE_LENGTH_TRI_HDIV_BASIS  * (VALUE_LENGTH_1D - 1); // H(div,tri) x H(vol,line)

/*
 For all the extern methods defined below:
 - C mangling on my setup prepends a _
 - Fortran mangling is all lower-case, and both prepends a _ and appends a _
 - So, to specify the externally-defined Fortran subroutine, take its name as lowercase, and append a _
 */
// size[0] should be N
// size[1] should be N+1 (VALUE_LENGTH_1D)
// H^1
extern "C" void shape1dhseg_(const double &xi, const int &order, const int size[2],
                             int &dofCount, double values[VALUE_LENGTH_1D], double values_dx[VALUE_LENGTH_1D]);
// L^2
extern "C" void shape1dqseg_(const double &xi, const int &order, const int size[2],
                             int &dofCount, double values[VALUE_LENGTH_1D]);

// quadrilaterals
// "order" has five entries -- still have to figure out what each of these means
// orientations[4] has the four edge orientations
// size[0] should be the maximum polynomial order (N)
// size[1] should be VALUE_LENGTH_1D*VALUE_LENGTH_1D
// values_grad is a 2 by size[1] fortran array.
// H^1
extern "C" void shape2dhquad_(const double xi[2], const int order[5], const int orientations[4], const int size[2],
                              int &dofCount, double values[VALUE_LENGTH_1D*VALUE_LENGTH_1D], double values_grad[VALUE_LENGTH_1D*VALUE_LENGTH_1D*2]);
// H(curl)
extern "C" void shape2dequad_(const double xi[2], const int order[5], const int orientations[4], const int size[2],
                              int &dofCount, double values[VALUE_LENGTH_1D*VALUE_LENGTH_1D*2], double values_curl[VALUE_LENGTH_1D*VALUE_LENGTH_1D]);
// H(div)
extern "C" void shape2dvquad_(const double xi[2], const int order[5], const int orientations[4], const int size[2],
                              int &dofCount, double values[VALUE_LENGTH_1D*VALUE_LENGTH_1D*2], double values_curl[VALUE_LENGTH_1D*VALUE_LENGTH_1D]);
// H(vol)
extern "C" void shape2dqquad_(const double xi[2], const int order[5], const int size[2],
                              int &dofCount, double values[VALUE_LENGTH_1D*VALUE_LENGTH_1D]);

// triangles
// H(grad)
// shape2dhtri(Xi,Norder,NoriE,nsize, NrdofH,ShapH,GradH)
extern "C" void shape2dhtri_(const double xi[2], const int order[4], const int orientations[3], const int size[2],
                             int &dofCount, double values[VALUE_LENGTH_TRI_H1_BASIS], double values_grad[VALUE_LENGTH_TRI_H1_BASIS*2]);

// H(curl)
// shape2detri
// shape2DETri(X,Nord,NoriE,Nsize, NrdofE,ShapE,CurlE)
extern "C" void shape2detri_(const double xi[2], const int order[4], const int orientations[3], const int size[2],
                             int &dofCount, double values[VALUE_LENGTH_TRI_HCURL_BASIS*2], double values_grad[VALUE_LENGTH_TRI_HCURL_BASIS]);

// H(div)
// shape2dvtri
extern "C" void shape2dvtri_(const double xi[2], const int order[4], const int orientations[3], const int size[2],
                             int &dofCount, double values[VALUE_LENGTH_TRI_HDIV_BASIS*2], double values_grad[VALUE_LENGTH_TRI_HDIV_BASIS]);

// H(vol)
// shape2dqtri
// shape2DQTri(X,Nord,Nsize, NrdofQ,ShapQ)
extern "C" void shape2dqtri_(const double xi[2], const int order[4], const int size[2],
                             int &dofCount, double values[VALUE_LENGTH_TRI_HVOL_BASIS]);

// hexahedra
// H(grad)
extern "C" void shape3dhhexa_(const double xi[3], const int order[19], const int edgeOrientations[12], const int faceOrientations[6], const int size[2],
                              int &dofCount, double values[VALUE_LENGTH_1D*VALUE_LENGTH_1D*VALUE_LENGTH_1D], double values_grad[VALUE_LENGTH_1D*VALUE_LENGTH_1D*VALUE_LENGTH_1D*3]);
// H(curl)
extern "C" void shape3dehexa_(const double xi[3], const int order[19], const int edgeOrientations[12], const int faceOrientations[6], const int size[2],
                              int &dofCount, double values[VALUE_LENGTH_1D*VALUE_LENGTH_1D*VALUE_LENGTH_1D*3], double values_curl[VALUE_LENGTH_1D*VALUE_LENGTH_1D*VALUE_LENGTH_1D*3]);
// H(div)
extern "C" void shape3dvhexa_(const double xi[3], const int order[19], const int faceOrientations[6], const int size[2],
                              int &dofCount, double values[VALUE_LENGTH_1D*VALUE_LENGTH_1D*VALUE_LENGTH_1D*3], double values_div[VALUE_LENGTH_1D*VALUE_LENGTH_1D*VALUE_LENGTH_1D]);
// H(vol)
extern "C" void shape3dqhexa_(const double xi[3], const int order[19], const int size[2],
                              int &dofCount, double values[VALUE_LENGTH_1D*VALUE_LENGTH_1D*VALUE_LENGTH_1D]);

// tetrahedra
// H(grad)
extern "C" void shape3dhtet_(const double xi[3], const int order[11], const int edgeOrientations[6], const int faceOrientations[4], const int size[2],
                             int &dofCount, double values[VALUE_LENGTH_TET_H1_BASIS], double values_grad[VALUE_LENGTH_TET_H1_BASIS*3]);

// H(curl)
extern "C" void shape3detet_(const double xi[3], const int order[11], const int edgeOrientations[6], const int faceOrientations[4], const int size[2],
                             int &dofCount, double values[VALUE_LENGTH_TET_H1_BASIS*3], double values_curl[VALUE_LENGTH_TET_H1_BASIS*3]);

// H(div)
extern "C" void shape3dvtet_(const double xi[3], const int order[11], const int faceOrientations[4], const int size[2],
                             int &dofCount, double values[VALUE_LENGTH_TET_H1_BASIS*3], double values_div[VALUE_LENGTH_TET_H1_BASIS]);

// H(vol)
extern "C" void shape3dqtet_(const double xi[3], const int order[11], const int size[2],
                             int &dofCount, double values[VALUE_LENGTH_TET_H1_BASIS]);

// wedges (a.k.a. prisms)
// H(grad)
extern "C" void shape3dhpris_(const double xi[3], const int order[15], const int edgeOrientations[9], const int faceOrientations[5], const int size[2],
                              int &dofCount, double values[VALUE_LENGTH_WEDGE_H1_BASIS], double values_grad[VALUE_LENGTH_WEDGE_H1_BASIS*3]);

// H(curl)
extern "C" void shape3depris_(const double xi[3], const int order[15], const int edgeOrientations[9], const int faceOrientations[5], const int size[2],
                              int &dofCount, double values[VALUE_LENGTH_WEDGE_H1_BASIS*3], double values_curl[VALUE_LENGTH_WEDGE_H1_BASIS*3]);

// H(div)
extern "C" void shape3dvpris_(const double xi[3], const int order[15], const int faceOrientations[5], const int size[2],
                              int &dofCount, double values[VALUE_LENGTH_WEDGE_H1_BASIS*3], double values_div[VALUE_LENGTH_WEDGE_H1_BASIS]);

// H(vol)
extern "C" void shape3dqpris_(const double xi[3], const int order[15], const int size[2],
                              int &dofCount, double values[VALUE_LENGTH_WEDGE_H1_BASIS]);

// Pyramids
// H(grad)
extern "C" void shape3dhpyra_(const double xi[3], const int order[14], const int edgeOrientations[8], const int faceOrientations[5], const int size[2],
                              int &dofCount, double values[VALUE_LENGTH_WEDGE_H1_BASIS], double values_grad[VALUE_LENGTH_WEDGE_H1_BASIS*3]);

// H(curl)
extern "C" void shape3depris_(const double xi[3], const int order[15], const int edgeOrientations[9], const int faceOrientations[5], const int size[2],
                              int &dofCount, double values[VALUE_LENGTH_WEDGE_H1_BASIS*3], double values_curl[VALUE_LENGTH_WEDGE_H1_BASIS*3]);

// H(div)
extern "C" void shape3dvpris_(const double xi[3], const int order[15], const int faceOrientations[5], const int size[2],
                              int &dofCount, double values[VALUE_LENGTH_WEDGE_H1_BASIS*3], double values_div[VALUE_LENGTH_WEDGE_H1_BASIS]);

// H(vol)
extern "C" void shape3dqpris_(const double xi[3], const int order[15], const int size[2],
                              int &dofCount, double values[VALUE_LENGTH_WEDGE_H1_BASIS]);

// helper functions, etc. for setting up ESEAS inputs
static const int MODORDER = 100; // used for integer packing in ESEAS

void getQuadPolyOrderESEAS(int polyOrder[5], int edge1_polyOrder, int edge2_polyOrder, int edge3_polyOrder, int edge4_polyOrder, int x_polyOrder, int y_polyOrder);

void getQuadPolyOrderESEAS(int polyOrderArray[5], int polyOrder);

void getHexPolyOrderESEAS(int polyOrder[19], std::vector<int> twelveEdgePolyOrders, std::vector<std::pair<int,int>> sixFacePolyOrderPairs,
                          std::vector<int> threeHexPolyOrders);

void getHexPolyOrderESEAS  (int polyOrderArray[19], int polyOrder);
void getTetPolyOrderESEAS  (int polyOrderArray[11], int polyOrder);
void getTriPolyOrderESEAS  (int polyOrderArray[4],  int polyOrder);
void getWedgePolyOrderESEAS(int polyOrderArray[15], int polyOrder);

KOKKOS_INLINE_FUNCTION int fieldOrdinal_QuadH1(int fieldOrdinal_x, int fieldOrdinal_y, int polyOrder_x)
{
  return fieldOrdinal_y * (polyOrder_x + 1) + fieldOrdinal_x;
}

// column-major ordering for spaceDim x something array (e.g. gradients in ESEAS)
KOKKOS_INLINE_FUNCTION int flatVectorIndex(int spaceDim, int entryOrdinal, int spatialComponent)
{
  return entryOrdinal * spaceDim + spatialComponent;
}

template<typename ExecutionSpace, typename OutputScalar, typename PointScalar>
std::map<int,int> getESEASOrdinalMap(Intrepid2::Basis<ExecutionSpace,OutputScalar,PointScalar> basis)
{
  // ESEAS numbers bases according to the subcell topologies that they're associated with.
  // Our tensor-product implementations use a tensor-product ordering.
  // This method returns a map that takes as key our basis ordinal, and returns the corresponding ESEAS ordinal.
  shards::CellTopology cellTopo = basis.getBaseCellTopology();
  int spaceDim = cellTopo.getDimension();
  std::map<int,int> dofMap;
  int dofOrdinalESEAS = 0;
  for (int d=0; d<=spaceDim; d++)
  {
    int subcellCount = cellTopo.getSubcellCount(d);
    for (int subcell_ESEAS=0; subcell_ESEAS<subcellCount; subcell_ESEAS++)
    {
      // shards subcell numberings match ESEAS for lines, quadrilaterals
      // for hexahedra, the vertex and edge numberings match, but the face numberings differ
      // I don't yet know whether numberings match on other cell topologies
      int subcell_shards = subcell_ESEAS;
      if (d==2)
      {
        if (cellTopo.getBaseKey() == shards::Hexahedron<>::key)
        {
          // shards has order bottom, right, top, left, front, back
          // ESEAS  has order front, back, bottom, right, top, left
          unsigned shardsOrdinal[6] = {4,5,0,1,2,3};
          subcell_shards = shardsOrdinal[subcell_ESEAS];
        }
        else if (cellTopo.getBaseKey() == shards::Triangle<>::key)
        {
          // these agree
        }
        else if (cellTopo.getBaseKey() == shards::Wedge<>::key)
        {
          // vertices agree, and I *think* edges do, too.  But faces differ.
          const unsigned shardsOrdinal[5] {3,4,0,1,2}; // note that the bottom face has order (0,2,1) in Shards, compared with (0,1,2) in ESEAS (presumably; have not confirmed this yet).  And the hypotenuse quad face is ordered as (0,3,5,2) in Shards, compared with (2,0,3,5) in ESEAS (presumably).
          subcell_shards = shardsOrdinal[subcell_ESEAS];
        }
        else if (cellTopo.getBaseKey() == shards::Pyramid<>::key)
        {
          // I believe these agree
//          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Still need to check ESEAS subcell ordering matches shards for Pyramid");
        }
      }
      int subcellDofCount = basis.getDofCount(d, subcell_shards);
      for (int subcellDofOrdinal=0; subcellDofOrdinal<subcellDofCount; subcellDofOrdinal++)
      {
        int intrepid2Ordinal = basis.getDofOrdinal(d, subcell_shards, subcellDofOrdinal);
        dofMap[intrepid2Ordinal] = dofOrdinalESEAS++;
      }
    }
  }
  return dofMap;
}

std::map<int,int> getESEASOrdinalMap(shards::CellTopology cellTopo, Intrepid2::EFunctionSpace fs, int polyOrder_x, int polyOrder_y, int polyOrder_z = -1);

#endif /* Intrepid2_ESEAS_Interface_h */
