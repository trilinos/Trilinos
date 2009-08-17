// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER


/** \file
    \brief  Test of the CellTools class.
    \author Created by P. Bochev, D. Ridzal and K. Peterson
*/
#include "Intrepid_CellTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ScalarTraits.hpp"

using namespace std;
using namespace Intrepid;


/** \brief  Maps the vertices of the subcell parametrization domain to that subcell. 
            
            Parametrization tests check if the vertices of the parametrization domain are properly 
            mapped to vertices of the resepective reference subcell. Because the parametrization map 
            is a polynomial whose degree depends on the number of vertices, if all vertices are 
            mapped correctly this is sufficient to assert that the parametrization map is correct.
  
            To test reference cells with two different kinds of faces, there are two argument slots
            to pass vertices of parametrization domains "A" and "B". When testing edges always pass
            the same argument because edges have the same parametrization domain [-1,1] (1-cube)

    \param  errorFlag       [out] - counts number of errors
    \param  parentCell      [in]  - topology of the reference cell whose 1 and 2-subcells are parametrized
    \param  subcParamVert_A [in]  - vertex coordinates of parametrization domain "A" (2-simplex for faces)
    \param  subcParamVert_A [in]  - vertex coordinates of parametrization domain "B" (2-cube for faces)
    \param  subcDim         [in]  - dimension of the subcells whose parametrizations are tested
    \param  outStream       [in]  - output stream to write
*/
void testSubcellParametrizations(int&                               errorFlag,
                                 const shards::CellTopology&        parentCell,
                                 const FieldContainer<double>&      subcParamVert_A,
                                 const FieldContainer<double>&      subcParamVert_B,
                                 const int                          subcDim,
                                 const Teuchos::RCP<std::ostream>&  outStream);
  

int main(int argc, char *argv[]) {
 
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  typedef CellTools<double>       CellTools;
  typedef shards::CellTopology    CellTopology;
  
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|                              Unit Test CellTools                            |\n" \
    << "|                                                                             |\n" \
    << "|     1) Edge parametrizations                                                |\n" \
    << "|     2) Face parametrizations                                                |\n" \
    << "|     3) Edge tangents                                                        |\n" \
    << "|     4) Face tangents and normals                                            |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  
  int errorFlag  = 0;

    
  // Vertices of the parametrization domain for 1-subcells: standard 1-cube [-1,1]
  FieldContainer<double> cube_1(2, 1);
  cube_1(0,0) = -1.0; 
  cube_1(1,0) = 1.0;

  
  // Vertices of the parametrization domain for triangular faces: the standard 2-simplex
  FieldContainer<double> simplex_2(3, 2);
  simplex_2(0, 0) = 0.0;   simplex_2(0, 1) = 0.0;
  simplex_2(1, 0) = 1.0;   simplex_2(1, 1) = 0.0;
  simplex_2(2, 0) = 0.0;   simplex_2(2, 1) = 1.0;
  
  
  // Vertices of the parametrization domain for quadrilateral faces: the standard 2-cube
  FieldContainer<double> cube_2(4, 2);
  cube_2(0, 0) =  -1.0;    cube_2(0, 1) =  -1.0;
  cube_2(1, 0) =   1.0;    cube_2(1, 1) =  -1.0;
  cube_2(2, 0) =   1.0;    cube_2(2, 1) =   1.0;
  cube_2(3, 0) =  -1.0;    cube_2(3, 1) =   1.0;

  
  // Pull all available topologies from Shards
  std::vector<shards::CellTopology> allTopologies;
  shards::getTopologies(allTopologies);
  
  
  // Set to 1 for edge and 2 for face tests
  int subcDim;

  try{
    
    *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| Test 1: edge parametrizations:                                              |\n"\
    << "===============================================================================\n\n";
    
    subcDim      = 1;
        
    // Loop over the cell topologies
    for(int topoOrd = 0; topoOrd < (int)allTopologies.size(); topoOrd++){
            
      // Test only 2D and 3D topologies that have reference cells, e.g., exclude Line, Pentagon, etc.
      if(allTopologies[topoOrd].getDimension() > 1 && CellTools::hasReferenceCell(allTopologies[topoOrd]) ){
        *outStream << " Testing edge parametrization for " <<  allTopologies[topoOrd].getName() <<"\n";
        testSubcellParametrizations(errorFlag,
                                    allTopologies[topoOrd],
                                    cube_1,
                                    cube_1,
                                    subcDim,
                                    outStream);
      }
    }

    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| Test 2: face parametrizations:                                              |\n"\
      << "===============================================================================\n\n";
    
    subcDim      = 2;
    
    // Loop over the cell topologies
    for(int topoOrd = 0; topoOrd < (int)allTopologies.size(); topoOrd++){
      
      // Test only 3D topologies that have reference cells
      if(allTopologies[topoOrd].getDimension() > 2 && CellTools::hasReferenceCell(allTopologies[topoOrd]) ){
        *outStream << " Testing face parametrization for cell topology " <<  allTopologies[topoOrd].getName() <<"\n";
        testSubcellParametrizations(errorFlag,
                                    allTopologies[topoOrd],
                                    simplex_2,
                                    cube_2,
                                    subcDim,
                                    outStream);
      }
    }
    
    
    /***********************************************************************************************
      *
      * Common for test 3 and 4: edge tangents and face normals for standard cells with base topo
      *
      **********************************************************************************************/
    
    // Allocate storage and extract all standard cells with base topologies
    std::vector<shards::CellTopology> standardBaseTopologies;    
    shards::getTopologies(standardBaseTopologies, 4, shards::STANDARD_CELL, shards::BASE_TOPOLOGY);

    // Define topologies for the edge and face parametrization domains. (faces are Tri or Quad)
    CellTopology paramEdge    (shards::getCellTopologyData<shards::Line<2> >() );
    CellTopology paramTriFace (shards::getCellTopologyData<shards::Triangle<3> >() );
    CellTopology paramQuadFace(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
    
    // Define CubatureFactory:
    DefaultCubatureFactory<double>  cubFactory;   
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| Test 3: edge tangents/normals for stand. cells with base topologies:        |\n"\
      << "===============================================================================\n\n";
    // This test loops over standard cells with base topologies, creates a set of nodes and tests tangents/normals 
    std::vector<shards::CellTopology>::iterator cti;
    
    // Define cubature on the edge parametrization domain:
    Teuchos::RCP<Cubature<double> > edgeCubature = cubFactory.create(paramEdge, 6); 
    int cubDim       = edgeCubature -> getDimension();
    int numCubPoints = edgeCubature -> getNumPoints();

    // Allocate storage for cubature points and weights on edge parameter domain and fill with points:
    FieldContainer<double> paramEdgePoints(numCubPoints, cubDim);
    FieldContainer<double> paramEdgeWeights(numCubPoints);
    edgeCubature -> getCubature(paramEdgePoints, paramEdgeWeights);
    

    // Loop over admissible topologies 
    for(cti = standardBaseTopologies.begin(); cti !=standardBaseTopologies.end(); ++cti){
      
      // Exclude 0D (node), 1D (Line) and Pyramid<5> cells
      if( ( (*cti).getDimension() >= 2) && ( (*cti).getKey() != shards::Pyramid<5>::key) ){ 
        
        int cellDim = (*cti).getDimension();
        int vCount  = (*cti).getVertexCount();
        FieldContainer<double> refCellVertices(vCount, cellDim);
        CellTools::getReferenceSubcellVertices(refCellVertices, cellDim, 0, (*cti) );
        
        *outStream << " Testing edge tangents";
          if(cellDim == 2) { *outStream << " and normals"; }          
        *outStream <<" for cell topology " <<  (*cti).getName() <<"\n";
        
        
        // Array for physical cell vertices ( must have rank 3 for cetJacobians)
        FieldContainer<double> physCellVertices(1, vCount, cellDim);

        // Randomize reference cell vertices by moving them up to +/- (1/8) units along their
        // coordinate axis. Guaranteed to be non-degenerate for standard cells with base topology 
        for(int v = 0; v < vCount; v++){
          for(int d = 0; d < cellDim; d++){
            double delta = Teuchos::ScalarTraits<double>::random()/8.0;
            physCellVertices(0, v, d) = refCellVertices(v, d) + delta;
          } //for d
        }// for v     
        
        // Allocate storage for cub. points on a ref. edge; Jacobians, phys. edge tangents/normals and 
        // benchmark tangents. The benchmark tangents are stored in array sized for a general case of  
        // non-affine edges to make it easier to add these tests later.
        FieldContainer<double> refEdgePoints(numCubPoints, cellDim);        
        FieldContainer<double> edgePointsJacobians(1, numCubPoints, cellDim, cellDim);
        FieldContainer<double> edgePointTangents(1, numCubPoints, cellDim);
        FieldContainer<double> edgePointNormals(1, numCubPoints, cellDim);
        FieldContainer<double> edgeBenchmarkTangents(1, numCubPoints, cellDim);
        

        // Loop over edges:
        for(int edgeOrd = 0; edgeOrd < (*cti).getEdgeCount(); edgeOrd++){
          
          /* 
           * Compute tangents on the specified edge using CellTools:
           *    1. Map points from edge parametrization domain to ref. edge with specified ordinal
           *    2. Compute parent cell Jacobians at ref. edge points
           *    3. Compute physical edge tangents
           */
          CellTools::mapToReferenceSubcell(refEdgePoints, paramEdgePoints, 1, edgeOrd, (*cti) );
          CellTools::setJacobian(edgePointsJacobians, refEdgePoints, physCellVertices, (*cti) );
          CellTools::getPhysicalEdgeTangents(edgePointTangents, edgePointsJacobians, edgeOrd, (*cti)); 
          
          /*
           * Compute tangents on the specified edge directly and compare with CellTools tangents.
           *    1. Get edge vertices
           *    2. For affine edges tangent coordinates are given by s'(t) = (v1-v0)/2
           *       (for now we only test affine edges, but later we will test edges for cells 
           *        with extended topologies.)
           */
          int v0ord = (*cti).getNodeMap(1, edgeOrd, 0);
          int v1ord = (*cti).getNodeMap(1, edgeOrd, 1);
          
          for(int pt = 0; pt < numCubPoints; pt++){
            for(int d = 0; d < cellDim; d++){
              edgeBenchmarkTangents(0, pt, d) = (physCellVertices(0, v1ord, d) - physCellVertices(0, v0ord, d))/2.0;
              
              // Compare with d-component of edge tangent by CellTools
              if( abs(edgeBenchmarkTangents(0, pt, d) - edgePointTangents(0, pt, d)) > INTREPID_THRESHOLD ){
                errorFlag++;
                *outStream
                  << std::setw(70) << "^^^^----FAILURE!" << "\n"
                  << " Edge tangent computation by CellTools failed for: \n"
                  << "       Cell Topology = " << (*cti).getName() << "\n"
                  << "        Edge ordinal = " << edgeOrd << "\n"
                  << "  Tangent coordinate = " << d << " edge point number = " << pt << "\n"
                  << "     CellTools value = " <<  edgePointTangents(0, pt, d)
                  << "     Benchmark value = " <<  edgeBenchmarkTangents(0, pt, d) << "\n\n";
              }
            } // for d
            
            // Test side normals for 2D cells only: edge normal has coordinates (t1, -t0)
            if(cellDim == 2) {
              CellTools::getPhysicalSideNormals(edgePointNormals, edgePointsJacobians, edgeOrd, (*cti));
              if( abs(edgeBenchmarkTangents(0, pt, 1) - edgePointNormals(0, pt, 0)) > INTREPID_THRESHOLD ){
                errorFlag++;
                *outStream
                  << std::setw(70) << "^^^^----FAILURE!" << "\n"
                  << " Edge Normal computation by CellTools failed for: \n"
                  << "       Cell Topology = " << (*cti).getName() << "\n"
                  << "        Edge ordinal = " << edgeOrd << "\n"
                  << "  Normal coordinate = " << 0 << " edge point number = " << pt << "\n"
                  << "     CellTools value = " <<  edgePointNormals(0, pt, 0)
                  << "     Benchmark value = " <<  edgeBenchmarkTangents(0, pt, 1) << "\n\n";
              }
              if( abs(edgeBenchmarkTangents(0, pt, 0) + edgePointNormals(0, pt, 1)) > INTREPID_THRESHOLD ){
                errorFlag++;
                *outStream
                  << std::setw(70) << "^^^^----FAILURE!" << "\n"
                  << " Edge Normal computation by CellTools failed for: \n"
                  << "       Cell Topology = " << (*cti).getName() << "\n"
                  << "        Edge ordinal = " << edgeOrd << "\n"
                  << "  Normal coordinate = " << 1 << " edge point number = " << pt << "\n"
                  << "     CellTools value = " <<  edgePointNormals(0, pt, 1)
                  << "     Benchmark value = " << -edgeBenchmarkTangents(0, pt, 0) << "\n\n";
              }
            } // edge normals
            
            
          } // for pt
        }// for edgeOrd
      }// if admissible cell
    }// for cti
    
    
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| Test 4: face/side normals for stand. 3D cells with base topologies:         |                                                      |\n"\
      << "===============================================================================\n\n";
    // This test loops over standard 3D cells with base topologies, creates a set of nodes and tests normals 

    // Define cubature on the edge parametrization domain:
    Teuchos::RCP<Cubature<double> > triFaceCubature  = cubFactory.create(paramTriFace, 6); 
    Teuchos::RCP<Cubature<double> > quadFaceCubature = cubFactory.create(paramQuadFace, 6); 
    
    int faceCubDim           = triFaceCubature -> getDimension();
    int numTriFaceCubPoints  = triFaceCubature -> getNumPoints();
    int numQuadFaceCubPoints = quadFaceCubature -> getNumPoints();    
    
    // Allocate storage for cubature points and weights on face parameter domain and fill with points:
    FieldContainer<double> paramTriFacePoints(numTriFaceCubPoints, faceCubDim);
    FieldContainer<double> paramTriFaceWeights(numTriFaceCubPoints);
    FieldContainer<double> paramQuadFacePoints(numQuadFaceCubPoints, faceCubDim);
    FieldContainer<double> paramQuadFaceWeights(numQuadFaceCubPoints);
    
    triFaceCubature -> getCubature(paramTriFacePoints, paramTriFaceWeights);
    quadFaceCubature -> getCubature(paramQuadFacePoints, paramQuadFaceWeights);
    
    
  }// try
  
  //============================================================================================//
  // Wrap up test: check if the test broke down unexpectedly due to an exception                //
  //============================================================================================//
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  };
  
  
  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  
  return errorFlag;
}



void testSubcellParametrizations(int&                               errorFlag,
                                 const shards::CellTopology&        parentCell,
                                 const FieldContainer<double>&      subcParamVert_A,
                                 const FieldContainer<double>&      subcParamVert_B,
                                 const int                          subcDim,
                                 const Teuchos::RCP<std::ostream>&  outStream){
  
  // Get cell dimension and subcell count
  int cellDim      = parentCell.getDimension();
  int subcCount    = parentCell.getSubcellCount(subcDim);
  
  
  // Loop over subcells of the specified dimension
  for(int subcOrd = 0; subcOrd < subcCount; subcOrd++){
    int subcVertexCount = parentCell.getVertexCount(subcDim, subcOrd);
    
    
    // Storage for correct reference subcell vertices and for the images of the parametrization domain points
    FieldContainer<double> refSubcellVertices(subcVertexCount, cellDim);
    FieldContainer<double> mappedParamVertices(subcVertexCount, cellDim);
    
    
    // Retrieve correct reference subcell vertices
    CellTools<double>::getReferenceSubcellVertices(refSubcellVertices, subcDim, subcOrd, parentCell);
    
    
    // Map vertices of the parametrization domain to 1 or 2-subcell with ordinal subcOrd
    // For edges parametrization domain is always 1-cube passed as "subcParamVert_A"
    if(subcDim == 1) {
      CellTools<double>::mapToReferenceSubcell(mappedParamVertices,
                                               subcParamVert_A,
                                               subcDim,
                                               subcOrd,
                                               parentCell);
    }
    // For faces need to treat Triangle and Quadrilateral faces separately
    else if(subcDim == 2) {
      
      // domain "subcParamVert_A" is the standard 2-simplex  
      if(subcVertexCount == 3){
        CellTools<double>::mapToReferenceSubcell(mappedParamVertices,
                                                 subcParamVert_A,
                                                 subcDim,
                                                 subcOrd,
                                                 parentCell);
      }
      // Domain "subcParamVert_B" is the standard 2-cube
      else if(subcVertexCount == 4){
        CellTools<double>::mapToReferenceSubcell(mappedParamVertices,
                                                 subcParamVert_B,
                                                 subcDim,
                                                 subcOrd,
                                                 parentCell);
      }
    }
    
    // Compare the images of the parametrization domain vertices with the true vertices.
    for(int subcVertOrd = 0; subcVertOrd < subcVertexCount; subcVertOrd++){
      for(int dim = 0; dim <  cellDim; dim++){
        
        if(mappedParamVertices(subcVertOrd, dim) != refSubcellVertices(subcVertOrd, dim) ) {
          errorFlag++; 
          *outStream 
            << std::setw(70) << "^^^^----FAILURE!" << "\n"
            << " Cell Topology = " << parentCell.getName() << "\n"
            << " Parametrization of subcell " << subcOrd << " which is "
            << parentCell.getName(subcDim,subcOrd) << " failed for vertex " << subcVertOrd << ":\n"
            << " parametrization map fails to map correctly coordinate " << dim << " of that vertex\n\n";
          
        }//if
      }// for dim 
    }// for subcVertOrd      
  }// for subcOrd
  
}






