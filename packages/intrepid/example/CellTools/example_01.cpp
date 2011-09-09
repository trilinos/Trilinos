// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file
\brief  Example of the CellTools class.
\author Created by P. Bochev and D. Ridzal
*/
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Shards_CellTopology.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace std;
using namespace Intrepid;
using namespace shards;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  std::cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                   Example use of the CellTools class                        |\n" \
  << "|                                                                             |\n" \
  << "|     1) Using shards::CellTopology to get cell types and topology            |\n" \
  << "|     2) Using CellTools to get cell Jacobians and their inverses and dets    |\n" \
  << "|     3) Testing points for inclusion in reference and physical cells         |\n" \
  << "|     4) Mapping points to and from reference cells with base and extended    |\n" \
  << "|        topologies.                                                          |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov)                      |\n" \
  << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Shards's website:   http://trilinos.sandia.gov/packages/shards             |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| EXAMPLE 1: Query of cell types and topology                                 |\n"\
  << "===============================================================================\n";
  
  
  typedef CellTools<double>       CellTools;
  typedef RealSpaceTools<double>  RealSpaceTools;
  typedef shards::CellTopology    CellTopology;
  
  
  // Vector to hold cell topologies
  std::vector<CellTopology> shardsTopologies;

  
  // All 2D cell topologies in Shards:
  int cellDim = 2;
  shards::getTopologies(shardsTopologies, cellDim);
  std::cout << "Number of all " << cellDim 
    << "D Shards cell topologies = " << shardsTopologies.size() << "\n\n";
  
  for(unsigned i = 0; i < shardsTopologies.size(); i++){
    std::cout << shardsTopologies[i].getName() << "\n"; 
  }
  std::cout <<"\n";

  
  // All standard 2D cell topologies in Shards:
  shards::getTopologies(shardsTopologies, cellDim, 
                        shards::STANDARD_CELL);
  std::cout << "Number of all " 
    << cellDim << "D standard Shards cell topologies = " << shardsTopologies.size() << "\n\n";
  
  for(unsigned i = 0; i < shardsTopologies.size(); i++){
    std::cout << shardsTopologies[i].getName() << "\n"; 
  }
  std::cout <<"\n";
  
  
  // All standard 2D cells with base topologies in Shards:
  shards::getTopologies(shardsTopologies, cellDim, 
                        shards::STANDARD_CELL, 
                        shards::BASE_TOPOLOGY);
  std::cout << "Number of all " << cellDim 
    << "D standard Shards cells with base topologies = " << shardsTopologies.size() << "\n\n";
  
  for(unsigned i = 0; i < shardsTopologies.size(); i++){
    std::cout << shardsTopologies[i].getName() << "\n"; 
  }
  std::cout <<"\n";
  
  
  // All standard 2D cells with extended topologies in Shards:
  shards::getTopologies(shardsTopologies, cellDim, 
                        shards::STANDARD_CELL, 
                        shards::EXTENDED_TOPOLOGY);
  std::cout << "Number of all " << cellDim 
    << "D standard Shards cells with extended topologies = " << shardsTopologies.size() << "\n\n";
  
  for(unsigned i = 0; i < shardsTopologies.size(); i++){
    std::cout << shardsTopologies[i].getName() << "\n"; 
  }
  std::cout <<"\n";

  
  // All non-standard 2D cells with base topologies in Shards:
  shards::getTopologies(shardsTopologies, cellDim, 
                        shards::NONSTANDARD_CELL, 
                        shards::BASE_TOPOLOGY);
  std::cout << "Number of all " << cellDim 
    << "D non-standard Shards cells with base topologies = " << shardsTopologies.size() << "\n\n";
  
  for(unsigned i = 0; i < shardsTopologies.size(); i++){
    std::cout << shardsTopologies[i].getName() << "\n"; 
  }
  std::cout <<"\n";
  
 
  // All non-standard 2D cells with extended topologies in Shards:
  shards::getTopologies(shardsTopologies, cellDim, 
                        shards::NONSTANDARD_CELL, 
                        shards::EXTENDED_TOPOLOGY);
  std::cout << "Number of all " << cellDim 
    << "D non-standard Shards cells with extended topologies = " << shardsTopologies.size() << "\n\n";
  
  for(unsigned i = 0; i < shardsTopologies.size(); i++){
    std::cout << shardsTopologies[i].getName() << "\n"; 
  }
  std::cout <<"\n";
  
  // Finally, get all shards cell topologies and print them!
  shards::getTopologies(shardsTopologies);
  std::cout << "Number of all Shards cell topologies = " << shardsTopologies.size() << "\n\n";
  for(unsigned i = 0; i < shardsTopologies.size(); i++){
    std::cout << shardsTopologies[i] << "\n"; 
  }
  std::cout <<"\n";
  
  
  
cout \
<< "===============================================================================\n"\
<< "| EXAMPLE 2: Using CellTools to get Jacobian, Jacobian inverse & Jacobian det |\n"\
<< "===============================================================================\n";

// 4 triangles with basic triangle topology: number of nodes = number of vertices
CellTopology triangle_3(shards::getCellTopologyData<Triangle<3> >() );
int numCells = 4;
int numNodes = triangle_3.getNodeCount();
int spaceDim = triangle_3.getDimension();

// Rank-3 array with dimensions (C,N,D) for the node coordinates of 3 traingle cells
FieldContainer<double> triNodes(numCells, numNodes, spaceDim);

// Initialize node data: accessor is (cellOrd, nodeOrd, coordinateOrd)
triNodes(0, 0, 0) = 0.0;  triNodes(0, 0, 1) = 0.0;    // 1st triangle =  the reference tri
triNodes(0, 1, 0) = 1.0;  triNodes(0, 1, 1) = 0.0;
triNodes(0, 2, 0) = 0.0;  triNodes(0, 2, 1) = 1.0;

triNodes(1, 0, 0) = 1.0;  triNodes(1, 0, 1) = 1.0;    // 2nd triangle = 1st shifted by (1,1)
triNodes(1, 1, 0) = 2.0;  triNodes(1, 1, 1) = 1.0;
triNodes(1, 2, 0) = 1.0;  triNodes(1, 2, 1) = 2.0;

triNodes(2, 0, 0) = 0.0;  triNodes(2, 0, 1) = 0.0;    // 3rd triangle = flip 1st vertical
triNodes(2, 1, 0) =-1.0;  triNodes(2, 1, 1) = 0.0;
triNodes(2, 2, 0) = 0.0;  triNodes(2, 2, 1) = 1.0;

triNodes(3, 0, 0) = 2.0;  triNodes(3, 0, 1) = 1.0;    // 4th triangle = just a triangle
triNodes(3, 1, 0) = 3.0;  triNodes(3, 1, 1) = 0.5;
triNodes(3, 2, 0) = 3.5;  triNodes(3, 2, 1) = 2.0;


// Rank-2 array with dimensions (P,D) for some points on the reference triangle
int numRefPoints = 2;
FieldContainer<double> refPoints(numRefPoints, spaceDim);
refPoints(0,0) = 0.0;   refPoints(0,1) = 0.0;
refPoints(1,0) = 0.5;   refPoints(1,1) = 0.5;


// Rank-4 array (C,P,D,D) for the Jacobian and its inverse and Rank-2 array (C,P) for its determinant
FieldContainer<double> triJacobian(numCells, numRefPoints, spaceDim, spaceDim);
FieldContainer<double> triJacobInv(numCells, numRefPoints, spaceDim, spaceDim);
FieldContainer<double> triJacobDet(numCells, numRefPoints);

// Rank-4 and Rank-2 auxiliary arrays
FieldContainer<double> rank4Aux (numCells, numRefPoints, spaceDim, spaceDim);
FieldContainer<double> rank2Aux (numCells, numRefPoints);


// Methods to compute cell Jacobians, their inverses and their determinants
CellTools::setJacobian(triJacobian, refPoints, triNodes, triangle_3); 
CellTools::setJacobianInv(triJacobInv, triJacobian );
CellTools::setJacobianDet(triJacobDet, triJacobian );

// Checks: compute det(Inv(DF)) and Inv(Inv(DF))
RealSpaceTools::det(rank2Aux, triJacobInv);
RealSpaceTools::inverse(rank4Aux, triJacobInv);

// Print data
std::cout 
<< std::scientific<< std::setprecision(4)
<< std::right << std::setw(16) << "DF(P)" 
<< std::right << std::setw(30) << "Inv(DF(P))"
<< std::right << std::setw(30) << "Inv(Inv(DF(P)))\n";

for(int cellOrd = 0; cellOrd < numCells; cellOrd++){
  std::cout 
  << "===============================================================================\n"
  << "Cell " << cellOrd << "\n";
  for(int pointOrd = 0; pointOrd < numRefPoints; pointOrd++){
    std::cout << "Point = ("
    << std::setw(4) << std::right << refPoints(pointOrd, 0) << ","<< std::setw(4) << std::right<< refPoints(pointOrd,1) << ")\n";
    for(int row = 0; row < spaceDim; row++){
      std::cout 
      << std::setw(11) << std::right << triJacobian(cellOrd, pointOrd, row, 0) << " "
      << std::setw(11) << std::right << triJacobian(cellOrd, pointOrd, row, 1) 
      //
      << std::setw(16) << std::right << triJacobInv(cellOrd, pointOrd, row, 0) << " "
      << std::setw(11) << std::right << triJacobInv(cellOrd, pointOrd, row, 1)
      //
      << std::setw(16) << std::right << rank4Aux(cellOrd, pointOrd, row, 0) << " "
      << std::setw(11) << std::right << rank4Aux(cellOrd, pointOrd, row, 1) << "\n";
    }
    std::cout 
      << setw(5)<<std::left<< "Determinant:\n"
      << std::setw(11) << std::right << triJacobDet(cellOrd, pointOrd)
      << std::setw(28) << std::right << rank2Aux(cellOrd, pointOrd)
      << std::setw(28) << std::right << " product = " << triJacobDet(cellOrd, pointOrd)*rank2Aux(cellOrd, pointOrd);
    std::cout<< "\n\n";
  }
}


// 2 Quadrilateral cells with base topology 
CellTopology quad_4(shards::getCellTopologyData<Quadrilateral<4> >() );
numCells = 2;
numNodes = quad_4.getNodeCount();
spaceDim = quad_4.getDimension();

FieldContainer<double> quadNodes(numCells, numNodes, spaceDim);
// 1st QUAD
quadNodes(0,0,0) = 1.00;  quadNodes(0,0,1) = 1.00;               
quadNodes(0,1,0) = 2.00;  quadNodes(0,1,1) = 0.75;
quadNodes(0,2,0) = 1.75;  quadNodes(0,2,1) = 2.00;  
quadNodes(0,3,0) = 1.25;  quadNodes(0,3,1) = 2.00, 
// 2ND QUAD
quadNodes(1,0,0) = 2.00;  quadNodes(1,0,1) = 0.75;               
quadNodes(1,1,0) = 3.00;  quadNodes(1,1,1) = 1.25;
quadNodes(1,2,0) = 2.75;  quadNodes(1,2,1) = 2.25;
quadNodes(1,3,0) = 1.75;  quadNodes(1,3,1) = 2.00;



// 1 Hexahedron cell with base topology: number of nodes = number of vertices
CellTopology hex_8(shards::getCellTopologyData<Hexahedron<8> >() );
numCells = 1;
numNodes = hex_8.getNodeCount();
spaceDim = hex_8.getDimension();

FieldContainer<double> hexNodes(numCells, numNodes, spaceDim);
hexNodes(0,0,0) =
// bottom face vertices
hexNodes(0,0,0) = 1.00;   hexNodes(0,0,1) = 1.00;   hexNodes(0,0,2) = 0.00;          
hexNodes(0,1,0) = 2.00;   hexNodes(0,1,1) = 0.75;   hexNodes(0,1,2) =-0.25;
hexNodes(0,2,0) = 1.75;   hexNodes(0,2,1) = 2.00;   hexNodes(0,2,2) = 0.00;
hexNodes(0,3,0) = 1.25;   hexNodes(0,3,1) = 2.00;   hexNodes(0,3,1) = 0.25;
// top face vertices
hexNodes(0,4,0) = 1.25;   hexNodes(0,4,1) = 0.75;   hexNodes(0,4,2) = 0.75;          
hexNodes(0,5,0) = 1.75;   hexNodes(0,5,1) = 1.00;   hexNodes(0,5,2) = 1.00;
hexNodes(0,6,0) = 2.00;   hexNodes(0,6,1) = 2.00;   hexNodes(0,6,2) = 1.25;
hexNodes(0,7,0) = 1.00;   hexNodes(0,7,1) = 2.00;   hexNodes(0,7,2) = 1.00;



std::cout << std::setprecision(16) << "\n" \
<< "===============================================================================\n"\
<< "| EXAMPLE 3: Using single point inclusion test method                         |\n"\
<< "===============================================================================\n";

// Define cell topologies
CellTopology edge3(shards::getCellTopologyData<Line<3> >() );
CellTopology tri6 (shards::getCellTopologyData<Triangle<6> >() );
CellTopology quad9(shards::getCellTopologyData<Quadrilateral<9> >() );
CellTopology tet4 (shards::getCellTopologyData<Tetrahedron<> >() );
CellTopology hex27(shards::getCellTopologyData<Hexahedron<27> >() );
CellTopology wedge(shards::getCellTopologyData<Wedge<> >() );
CellTopology pyr  (shards::getCellTopologyData<Pyramid<> >() );

// Points that are close to the boundaries of their reference cells
double point_in_edge[1]    = {1.0-INTREPID_EPSILON};
double point_in_quad[2]    = {1.0,                  1.0-INTREPID_EPSILON};
double point_in_tri[2]     = {0.5-INTREPID_EPSILON, 0.5-INTREPID_EPSILON};
double point_in_tet[3]     = {0.5-INTREPID_EPSILON, 0.5-INTREPID_EPSILON, 2.0*INTREPID_EPSILON};
double point_in_hex[3]     = {1.0-INTREPID_EPSILON, 1.0-INTREPID_EPSILON, 1.0-INTREPID_EPSILON};
double point_in_wedge[3]   = {0.5,                  0.25,                 1.0-INTREPID_EPSILON};
double point_in_pyramid[3] = {-INTREPID_EPSILON,    INTREPID_EPSILON,     1.0-INTREPID_EPSILON};

// Run the inclusion test for each point and print results
int in_edge     = CellTools::checkPointInclusion( point_in_edge,    1, edge3, INTREPID_THRESHOLD );
int in_quad     = CellTools::checkPointInclusion( point_in_quad,    2, quad9 );
int in_tri      = CellTools::checkPointInclusion( point_in_tri,     2, tri6 );
int in_tet      = CellTools::checkPointInclusion( point_in_tet,     3, tet4 );
int in_hex      = CellTools::checkPointInclusion( point_in_hex,     3, hex27 );
int in_prism    = CellTools::checkPointInclusion( point_in_wedge,   3, wedge );
int in_pyramid  = CellTools::checkPointInclusion( point_in_pyramid, 3, pyr );

if(in_edge) {
  std::cout << "(" << point_in_edge[0] << ")" 
  << " is inside reference Line " << endl;
}
if(in_quad) {
  std::cout << "(" << point_in_quad[0] << "," << point_in_quad[1]  << ")" 
  << " is inside reference Quadrilateral " << endl;
}
if(in_tri) {
  std::cout << "(" << point_in_tri[0] << "," << point_in_tri[1] << ")"  
  << " is inside reference Triangle " << endl;
}
if(in_tet) {
  std::cout << "(" << point_in_tet[0] << "," << point_in_tet[1] << "," << point_in_tet[2]<<")" 
  << " is inside reference Tetrahedron " << endl;
}
if(in_hex) {
  std::cout << "(" << point_in_hex[0] << "," << point_in_hex[1] << "," << point_in_hex[2]<<")" 
  << " is inside reference Hexahedron " << endl;
}
if(in_prism) {
  std::cout << "(" << point_in_wedge[0] << "," << point_in_wedge[1] << "," << point_in_wedge[2]<<")" 
  << " is inside reference Wedge " << endl;
}
if(in_pyramid) {
  std::cout << "(" << point_in_pyramid[0] << "," << point_in_pyramid[1] << "," << point_in_pyramid[2]<<")" 
  << " is inside reference Pyramid " << endl;
}

// Change the points to be outside their reference cells.
double small = 2.0*INTREPID_THRESHOLD;
point_in_edge[0] += small;

point_in_tri[0] += small;
point_in_tri[1] += small;

point_in_pyramid[0] += small;
point_in_pyramid[1] += small;
point_in_pyramid[2] += small;

in_edge     = CellTools::checkPointInclusion(point_in_edge,    1,   edge3);
in_tri      = CellTools::checkPointInclusion(point_in_tri,     2,   tri6);
in_pyramid  = CellTools::checkPointInclusion(point_in_pyramid, 3,   pyr);

std::cout << "\nChecking if perturbed Points belong to reference cell: " << endl;
if(!in_edge) {
  std::cout << "(" << point_in_edge[0] << ")" << " is NOT inside reference Line " << endl;
}
if(!in_tri) {
  std::cout << "(" << point_in_tri[0] << "," << point_in_tri[1] << ")"  << " is NOT inside reference Triangle " << endl;
}
if(!in_pyramid) {
  std::cout << "(" << point_in_pyramid[0] << "," << point_in_pyramid[1] << "," << point_in_pyramid[2]<<")" 
  << " is NOT inside reference Pyramid " << endl;
}



std::cout << std::setprecision(6) << "\n" \
<< "===============================================================================\n"\
<< "| EXAMPLE 4-A: Using pointwise inclusion test method for reference cells      |\n"\
<< "===============================================================================\n";


// Rank-1 array for one 2D reference point and rank-1 array for the test result
FieldContainer<double> onePoint(2);
FieldContainer<int> testOnePoint(1);

onePoint(0) = 0.2;   onePoint(1) = 0.3;

std::cout <<"\t Pointwise inclusion test for Triangle<6>: rank-1 array with a single 2D point: \n";

CellTools::checkPointwiseInclusion(testOnePoint, onePoint, tri6);  

std::cout << "point(" 
<< std::setw(13) << std::right << onePoint(0) << "," 
<< std::setw(13) << std::right << onePoint(1) << ") ";
if( testOnePoint(0) ) {
  std::cout << " is inside. \n";
}
else{
  std::cout << " is not inside. \n";
}
std::cout << "\n";


// Rank-2 array for 4 2D reference points (vector of points) and rank-1 array for the test result
FieldContainer<double>  fourPoints(4, 2);
FieldContainer<int> testFourPoints(4);

fourPoints(0,0) = 0.5;   fourPoints(0,1) = 0.5;
fourPoints(1,0) = 1.0;   fourPoints(1,1) = 1.1;
fourPoints(2,0) =-1.0;   fourPoints(2,1) =-1.1;
fourPoints(3,0) =-1.0;   fourPoints(3,1) = 0.5;

std::cout <<"\t  Pointwise inclusion test for Quadrilateral<9>: rank-2 array with 4 2D points: \n";

CellTools::checkPointwiseInclusion(testFourPoints, fourPoints, quad9);

for(int i1 = 0; i1 < fourPoints.dimension(0); i1++) {
  std::cout << " point(" << i1 << ") = (" 
  << std::setw(13) << std::right << fourPoints(i1, 0) << "," 
  << std::setw(13) << std::right << fourPoints(i1, 1) << ") ";
  if( testFourPoints(i1) ) {
    std::cout << " is inside. \n";
  }
  else{
    std::cout << " is not inside. \n";
  }
}
std::cout << "\n";


// Rank-3 array for 6 2D points and rank-2 array for the test result
FieldContainer<double>  sixPoints(2, 3, 2);
FieldContainer<int> testSixPoints(2, 3);

sixPoints(0,0,0) = -1.0;   sixPoints(0,0,1) =  1.0;
sixPoints(0,1,0) =  1.0;   sixPoints(0,1,1) =  0.0;
sixPoints(0,2,0) =  0.0;   sixPoints(0,2,1) =  1.0;
sixPoints(1,0,0) = -1.0;   sixPoints(1,0,1) = -1.0;
sixPoints(1,1,0) =  0.1;   sixPoints(1,1,1) =  0.2;
sixPoints(1,2,0) =  0.2;   sixPoints(1,2,1) =  0.3;

std::cout <<"\t  Pointwise inclusion test for Triangle<6>: rank-3 array with six 2D points: \n";

CellTools::checkPointwiseInclusion(testSixPoints, sixPoints, tri6);

for(int i0 = 0; i0 < sixPoints.dimension(0); i0++){
  for(int i1 = 0; i1 < sixPoints.dimension(1); i1++) {
    std::cout << " point(" << i0 << "," << i1 << ") = (" 
    << std::setw(13) << std::right << sixPoints(i0, i1, 0) << "," 
    << std::setw(13) << std::right << sixPoints(i0, i1, 1) << ") ";
    if( testSixPoints(i0, i1) ) {
      std::cout << " is inside. \n";
    }
    else{
      std::cout << " is not inside. \n";
    }
    
  }
}
std::cout << "\n";


// Rank-3 array for 6 3D reference points and rank-2 array for the test results
FieldContainer<double> six3DPoints(2, 3, 3);
FieldContainer<int> testSix3DPoints(2, 3);

six3DPoints(0,0,0) = -1.0;   six3DPoints(0,0,1) =  1.0;   six3DPoints(0,0,2) =  1.0;
six3DPoints(0,1,0) =  1.0;   six3DPoints(0,1,1) =  1.0;   six3DPoints(0,1,2) = -1.0;
six3DPoints(0,2,0) =  0.0;   six3DPoints(0,2,1) =  1.1;   six3DPoints(0,2,2) =  1.0;
six3DPoints(1,0,0) = -1.1;   six3DPoints(1,0,1) = -1.0;   six3DPoints(1,0,2) = -1.0;
six3DPoints(1,1,0) =  0.1;   six3DPoints(1,1,1) =  0.2;   six3DPoints(1,1,2) =  0.2;
six3DPoints(1,2,0) =  1.1;   six3DPoints(1,2,1) =  0.3;   six3DPoints(1,2,2) =  0.3;

std::cout <<"\t  Pointwise inclusion test for Hexahedron<27>: rank-3 array with six 3D points: \n";

CellTools::checkPointwiseInclusion(testSix3DPoints, six3DPoints, hex27);



for(int i0 = 0; i0 < six3DPoints.dimension(0); i0++){
  for(int i1 = 0; i1 < six3DPoints.dimension(1); i1++) {
    std::cout << " point(" << i0 << "," << i1 << ") = (" 
    << std::setw(13) << std::right << six3DPoints(i0, i1, 0) << "," 
    << std::setw(13) << std::right << six3DPoints(i0, i1, 1) << "," 
    << std::setw(13) << std::right << six3DPoints(i0, i1, 2) << ") ";
    if( testSix3DPoints(i0, i1) ) {
      std::cout << " is inside. \n";
    }
    else{
      std::cout << " is not inside. \n";
    }
  }
}


std::cout << std::setprecision(6) << "\n" \
<< "===============================================================================\n"\
<< "| EXAMPLE 4-B: Pointwise inclusion test for a physical cell in cell workset   |\n"\
<< "===============================================================================\n";

// Rank-2 array for a single set of 5 2D physical points and rank-1 array for the test results
FieldContainer<double>  fivePoints(5,2);
FieldContainer<int> testFivePoints(5);

// These points will be tested for inclusion in the last Triangle cell specified by triNodes
fivePoints(0, 0) = 2.1 ;   fivePoints(0, 1) = 1.0 ;       // in
fivePoints(1, 0) = 3.0 ;   fivePoints(1, 1) = 0.75;       // in
fivePoints(2, 0) = 3.5 ;   fivePoints(2, 1) = 1.9 ;       // out
fivePoints(3, 0) = 2.5 ;   fivePoints(3, 1) = 1.0 ;       // in
fivePoints(4, 0) = 2.75;   fivePoints(4, 1) = 2.0 ;       // out

CellTools::checkPointwiseInclusion(testFivePoints, fivePoints, triNodes, triangle_3,  3);

std::cout << " Vertices of Triangle #3: \n" 
<< "\t(" << triNodes(3, 0, 0) << ", " << triNodes(3, 0, 1) << ")\n"
<< "\t(" << triNodes(3, 1, 0) << ", " << triNodes(3, 1, 1) << ")\n"
<< "\t(" << triNodes(3, 1, 0) << ", " << triNodes(3, 1, 1) << ")\n"
<< " Inclusion test results for the physical points: \n\n";

for(int i1 = 0; i1 < fivePoints.dimension(0); i1++) {
  std::cout << " point(" << i1 << ") = (" 
  << std::setw(13) << std::right << fivePoints(i1, 0) << "," 
  << std::setw(13) << std::right << fivePoints(i1, 1) << ") ";
  if( testFivePoints(i1) ) {
    std::cout << " is inside. \n";
  }
  else{
    std::cout << " is not inside. \n";
  }
}
std::cout << "\n";



std::cout << std::setprecision(6) << "\n" \
<< "===============================================================================\n"\
<< "| EXAMPLE 4-C: Pointwise inclusion test for all physical cells in cell workset|\n"\
<< "===============================================================================\n";

// Rank-3 array for 4 sets of 2 2D physical points and rank-2 array for the test results
FieldContainer<double>  fourPointSets(4, 2, 2);
FieldContainer<int> testFourSets(4, 2);

// 1st point set - will be tested for inclusion in Triangle #0
fourPointSets(0, 0, 0) = 0.25;     fourPointSets(0, 0, 1) = 0.75;       // in
fourPointSets(0, 1, 0) = 0.00;     fourPointSets(0, 1, 1) =-0.01;       // out

// 2nd point set - will be tested for inclusion in Triangle #1
fourPointSets(1, 0, 0) = 1.50;     fourPointSets(1, 0, 1) = 1.50;       // in
fourPointSets(1, 1, 0) = 0.99;     fourPointSets(1, 1, 1) = 1.50;       // out

// 3rd point set - will be tested for inclusion in Triangle #2
fourPointSets(2, 0, 0) =-0.25;     fourPointSets(2, 0, 1) = 0.70;       // in
fourPointSets(2, 1, 0) = 0.0001;   fourPointSets(2, 1, 1) = 0.50;       // out

// 4th point set - will be tested for inclusion in Triangle #3
fourPointSets(3, 0, 0) = 3.00;     fourPointSets(3, 0, 1) = 1.00;       // in
fourPointSets(3, 1, 0) = 3.50;     fourPointSets(3, 1, 1) = 0.50;       // out

CellTools::checkPointwiseInclusion(testFourSets, fourPointSets, triNodes, triangle_3);

for(int cell = 0; cell < triNodes.dimension(0); cell++){

  std::cout << " Testing point set inclusion for Triangle #" << cell << " with vertices \n" 
  << "\t(" << triNodes(cell, 0, 0) << ", " << triNodes(cell, 0, 1) << ")\n"
  << "\t(" << triNodes(cell, 1, 0) << ", " << triNodes(cell, 1, 1) << ")\n"
  << "\t(" << triNodes(cell, 1, 0) << ", " << triNodes(cell, 1, 1) << ")\n"
  << " Results for physical point set indexed by the cell ordinal  " << cell << " \n\n";
  
  for(int i1 = 0; i1 < fourPointSets.dimension(1); i1++) {
    std::cout << " point(" << i1 << ") = (" 
    << std::setw(13) << std::right << fourPointSets(cell, i1, 0) << "," 
    << std::setw(13) << std::right << fourPointSets(cell, i1, 1) << ") ";
    if( testFourSets(cell, i1) ) {
      std::cout << " is inside  Triangle #" << cell << " \n";
    }
    else{
      std::cout << " is outside Triangle #" << cell << " \n";
    }
  }
  if(cell < triNodes.dimension(0) - 1) {
    std::cout << "-------------------------------------------------------------------------------\n";
  }
  else{
    std::cout <<" \n"; 
  }
}// cell



std::cout << "\n" \
<< "===============================================================================\n"\
<< "| EXAMPLE 5: Using point set inclusion test method                            |\n"\
<< "===============================================================================\n";

std::cout <<"\t  Point set inclusion test for Triangle<6>: rank-2 array with four 2D point: \n";

if( CellTools::checkPointsetInclusion(fourPoints, tri6) ) {
  std::cout << "\t - All points are inside the reference Triangle<6> cell. \n\n ";
}
else{
  std::cout << "\t - At least one point is not inside the reference Triangle<6> cell. \n\n";
}

std::cout <<"\t  Point set inclusion test for Hexahedron<27>: rank-3 array with six 3D point: \n";

if( CellTools::checkPointsetInclusion(six3DPoints, hex27) ) {
  std::cout << "\t - All points are inside the reference Hexahedron<27> cell. \n\n ";
}
else{
  std::cout << "\t - At least one point is not inside the reference Hexahedron<27> cell. \n\n";
}



std::cout << std::setprecision(4) << "\n" \
<< "===============================================================================\n"\
<< "| EXAMPLE 6: mapping a single point set to physical cells with base topology  |\n"\
<< "===============================================================================\n";

// Rank-3 array with dimensions (P, D) for points on the reference triangle
FieldContainer<double> refTriPoints(3,triangle_3.getDimension() );
refTriPoints(0,0) = 0.2;  refTriPoints(0,1) = 0.0;      // on edge 0
refTriPoints(1,0) = 0.4;  refTriPoints(1,1) = 0.6;      // on edge 1
refTriPoints(2,0) = 0.0;  refTriPoints(2,1) = 0.8;      // on edge 2

// Reference points will be mapped to physical cells with vertices in triNodes: define the appropriate array
FieldContainer<double> physTriPoints(triNodes.dimension(0),      // cell count
                                     refTriPoints.dimension(0),     // point count
                                     refTriPoints.dimension(1));    // point dimension (=2)

CellTools::mapToPhysicalFrame(physTriPoints, refTriPoints, triNodes, triangle_3);

for(int cell = 0; cell < triNodes.dimension(0); cell++){
  std::cout << "====== Triangle " << cell << " ====== \n";
  for(int pt = 0; pt < refTriPoints.dimension(0); pt++){
    std::cout 
    <<  "(" 
    << std::setw(13) << std::right << refTriPoints(pt,0) << "," 
    << std::setw(13) << std::right << refTriPoints(pt,1) << ") -> "
    <<  "(" 
    << std::setw(13) << std::right << physTriPoints(cell, pt, 0) << "," 
    << std::setw(13) << std::right << physTriPoints(cell, pt, 1) << ") \n"; 
  }
  std::cout << "\n";
}


std::cout << std::setprecision(4) << "\n" \
<< "===============================================================================\n"\
<< "| EXAMPLE 7: mapping a single point set to physical cells with ext. topology  |\n"\
<< "===============================================================================\n";
/*
 *  This example illustrates reference-to-physical mapping for a triangle with curved sides. The
 *  physical curved triangle is defined as follows:
 *
 *    - edge 0: starts at (1, -1/2) ends at (1,1)     and is a vertical line segment
 *    - edge 1: starts at (1,1)     ends at (0,0)     and is parabolic segment lying on y=x^2
 *    - edge 2: starts at (0,0)     ends at (1,-1/2)  and is parabolic segment lying on y = -1/2 x^2
 * 
 *  The triangle is uniquely specified by its 3 vertices and 3 edge points. Therefore, to compute the
 *  map from reference to physical coordinates we specify Triangle<6> topology.
 */

// A (C,V,D) array with the 6 nodes of the physical Triangle<6>
FieldContainer<double> tri6Nodes(1,6,2);
tri6Nodes(0,0,0) = 1.0;    tri6Nodes(0,0,1) = -0.5;
tri6Nodes(0,1,0) = 1.0;    tri6Nodes(0,1,1) =  1.0;
tri6Nodes(0,2,0) = 0.0;    tri6Nodes(0,2,1) =  0.0;

tri6Nodes(0,3,0) = 1.0;    tri6Nodes(0,3,1) =  0.0;
tri6Nodes(0,4,0) = 0.5;    tri6Nodes(0,4,1) =  0.25;
tri6Nodes(0,5,0) = 0.5;    tri6Nodes(0,5,1) = -0.125;

// A (P, D) array with the 6 nodes of the reference Triangle<6> plus some extra points
FieldContainer<double> refTri6Points(9,2);
refTri6Points(0,0) = 0.0;    refTri6Points(0,1) =  0.0;
refTri6Points(1,0) = 1.0;    refTri6Points(1,1) =  0.0;
refTri6Points(2,0) = 0.0;    refTri6Points(2,1) =  1.0;

refTri6Points(3,0) = 0.5;    refTri6Points(3,1) =  0.0;
refTri6Points(4,0) = 0.5;    refTri6Points(4,1) =  0.5;
refTri6Points(5,0) = 0.0;    refTri6Points(5,1) =  0.5;

refTri6Points(6,0) = 0.75;   refTri6Points(6,1) =  0.0;               // on ref edge 0
refTri6Points(7,0) = 0.25;   refTri6Points(7,1) =  0.75;              // on ref edge 1
refTri6Points(8,0) = 0.00;   refTri6Points(8,1) =  0.25;              // on ref edge 2

// A (C,P,D) array for the images of the reference points in physical frame
FieldContainer<double> physTri6Points(tri6Nodes.dimension(0),             // cell count
                                      refTri6Points.dimension(0),         // point count
                                      refTri6Points.dimension(1));        // point dimension (=2)

// Define the cell topology: use the extended Triangle<6> topology to fit the curved edges
CellTopology triangle_6(shards::getCellTopologyData<Triangle<6> >() );

// 
CellTools::mapToPhysicalFrame(physTri6Points, refTri6Points, tri6Nodes, triangle_6);

for(int cell = 0; cell < tri6Nodes.dimension(0); cell++){
  std::cout << "====== Triangle " << cell << " ====== \n";
  for(int pt = 0; pt < refTri6Points.dimension(0); pt++){
    std::cout 
    <<  "(" 
    << std::setw(13) << std::right << refTri6Points(pt,0) << "," 
    << std::setw(13) << std::right << refTri6Points(pt,1) << ") -> "
    <<  "(" 
    << std::setw(13) << std::right << physTri6Points(cell, pt, 0) << "," 
    << std::setw(13) << std::right << physTri6Points(cell, pt, 1) << ") \n"; 
  }
  std::cout << "\n";
}



std::cout << "\n" \
<< "===============================================================================\n"\
<< "| EXAMPLE 8: mapping a single physical point set to reference frame           |\n"\
<< "===============================================================================\n";
/*
 * This example shows use of mapToReferenceFrame with rank-2 array (P,D) of points in physical frame.
 * Points are mapped back to reference frame using the mapping for one of the physical cells whose
 * nodes are passed as an argument. Therefore, this use case requires a valid cell ordinal (relative
 * to the nodes array) to be specified. The output array for the images of the points is also rank-2 (P,D).
 *
 */
// Rank-2 arrays with dimensions (P, D) for physical points and their preimages
FieldContainer<double> physPoints(5, triangle_3.getDimension() ); 
FieldContainer<double> preImages (5, triangle_3.getDimension() ); 

// First 3 points are the vertices of the last triangle (ordinal = 3) stored in triNodes
physPoints(0,0) = triNodes(3, 0, 0);   physPoints(0,1) = triNodes(3, 0, 1) ;
physPoints(1,0) = triNodes(3, 1, 0);   physPoints(1,1) = triNodes(3, 1, 1);
physPoints(2,0) = triNodes(3, 2, 0);   physPoints(2,1) = triNodes(3, 2, 1);

// last 2 points are just some arbitrary points contained in the last triangle
physPoints(3,0) = 3.0;                    physPoints(3,1) = 1.0;
physPoints(4,0) = 2.5;                    physPoints(4,1) = 1.1;


// Map physical points from triangle 3 to the reference frame
int whichCell = 3;
CellTools::mapToReferenceFrame(preImages, physPoints, triNodes, triangle_3, whichCell  );

std::cout << " Mapping from Triangle #"<< whichCell << " with vertices: \n" 
<< "\t(" << triNodes(whichCell, 0, 0) << ", " << triNodes(whichCell, 0, 1) << ")\n"
<< "\t(" << triNodes(whichCell, 1, 0) << ", " << triNodes(whichCell, 1, 1) << ")\n"
<< "\t(" << triNodes(whichCell, 1, 0) << ", " << triNodes(whichCell, 1, 1) << ")\n\n"
<< " Physical points and their reference cell preimages: \n";

for(int pt = 0; pt < physPoints.dimension(0); pt++){
  std::cout 
  <<  "(" << std::setw(13) << std::right << physPoints(pt,0) << "," 
  << std::setw(13) << std::right << physPoints(pt,1) << ") -> "
  <<  "(" 
  << std::setw(13) << std::right << preImages(pt, 0) << "," 
  << std::setw(13) << std::right << preImages(pt, 1) << ") \n"; 
}
std::cout << "\n";

// As a check, map pre-images back to Triangle #3
FieldContainer<double> images(5, triangle_3.getDimension() );
CellTools::mapToPhysicalFrame(images, preImages, triNodes, triangle_3, whichCell);

std::cout << " Check: map preimages back to Triangle #3: \n"; 
for(int pt = 0; pt < images.dimension(0); pt++){
  std::cout 
  <<  "(" << std::setw(13) << std::right << preImages(pt,0) << "," 
  << std::setw(13) << std::right << preImages(pt,1) << ") -> "
  <<  "(" 
  << std::setw(13) << std::right << images(pt, 0) << "," 
  << std::setw(13) << std::right << images(pt, 1) << ") \n"; 
}
std::cout << "\n";


std::cout << "\n" \
<< "===============================================================================\n"\
<< "| EXAMPLE 9: mapping physical point sets on a cell workset to reference frame |\n"\
<< "===============================================================================\n";
/*
 * This example shows use of mapToReferenceFrame with rank-3 array (C,P,D) of points in physical frame.
 * For each cell ordinal (relative to the nodes array), the associated point set is mapped back to 
 * reference frame using the mapping corresponding to the physical cell with that ordinal. This use case
 * requires the default value -1 for the cell ordinal. The output array for the images of the points 
 * is also rank-3 (C,P,D).
 *
 */
// Rank-3 arrays with dimensions (C, P, D) for physical points and their preimages
FieldContainer<double> physPointSets(triNodes.dimension(0), 2, triangle_3.getDimension() ); 
FieldContainer<double> preImageSets (triNodes.dimension(0), 2, triangle_3.getDimension() ); 

// Point set on Triangle #0
physPointSets(0,0,0) = 0.25; physPointSets(0,0,1) = 0.25; 
physPointSets(0,1,0) = 0.50; physPointSets(0,1,1) = 0.50; 

// Point set on Triangle #1
physPointSets(1,0,0) = 1.00; physPointSets(1,0,1) = 1.00; 
physPointSets(1,1,0) = 1.25; physPointSets(1,1,1) = 1.75; 

// Point set on Triangle #2
physPointSets(2,0,0) =-0.25; physPointSets(2,0,1) = 0.25; 
physPointSets(2,1,0) =-0.50; physPointSets(2,1,1) = 0.50; 

// Point set on Triangle #3
physPointSets(3,0,0) = 3.00; physPointSets(3,0,1) = 1.00; 
physPointSets(3,1,0) = 2.00; physPointSets(3,1,1) = 1.00; 


// Map physical point sets to reference frame: requires default value of cell ordinal (skip last arg)
CellTools::mapToReferenceFrame(preImageSets, physPointSets, triNodes, triangle_3);

for(int cell = 0; cell < triNodes.dimension(0); cell++){
  std::cout << "====== Triangle " << cell << " ====== \n";
  for(int pt = 0; pt < physPointSets.dimension(1); pt++){
    std::cout 
    <<  "(" 
    << std::setw(13) << std::right << physPointSets(cell, pt, 0) << "," 
    << std::setw(13) << std::right << physPointSets(cell, pt, 1) << ") -> "
    <<  "(" 
    << std::setw(13) << std::right << preImageSets(cell, pt, 0) << "," 
    << std::setw(13) << std::right << preImageSets(cell, pt, 1) << ") \n"; 
  }
  std::cout << "\n";
}
  
// As a check, map preImageSets back to physical frame
FieldContainer<double> postImageSets(triNodes.dimension(0), 2, triangle_3.getDimension() );

CellTools::mapToPhysicalFrame(postImageSets, preImageSets, triNodes, triangle_3);

std::cout << " Check: map preimages back to Triangles: \n"; 
for(int cell = 0; cell < triNodes.dimension(0); cell++){
  std::cout << "====== Triangle " << cell << " ====== \n";
  for(int pt = 0; pt < preImageSets.dimension(1); pt++){
    std::cout 
    <<  "(" 
    << std::setw(13) << std::right << preImageSets(cell, pt, 0) << "," 
    << std::setw(13) << std::right << preImageSets(cell, pt, 1) << ") -> "
    <<  "(" 
    << std::setw(13) << std::right << postImageSets(cell, pt, 0) << "," 
    << std::setw(13) << std::right << postImageSets(cell, pt, 1) << ") \n"; 
  }
  std::cout << "\n";
}

return 0;
}




























