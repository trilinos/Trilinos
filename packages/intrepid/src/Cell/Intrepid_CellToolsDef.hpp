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


/** \file   Intrepid_CellToolsDef.hpp
    \brief  Definition file for the Intrepid::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

  template<class Scalar>
  const FieldContainer<double>& CellTools<Scalar>::getSubcellParametrization(const int subcellDim,
                                                                             const int subcellOrd,
                                                                             const shards::CellTopology& parentCell){
    
    // Coefficients of the coordinate functions defining the parametrization maps are stored in 
    // rank-3 arrays with dimensions (SC, PCD, COEF) where:
    //  - SC    is the subcell count of subcells with the specified dimension in the parent cell
    //  - PCD   is Parent Cell Dimension, which gives the number of coordinate functions in the map:
    //          PCD = 2 for standard 2D cells and non-standard 2D cells: shell line and beam
    //          PCD = 3 for standard 3D cells and non-standard 3D cells: shell Tri and Quad
    //  - COEF  is number of coefficients needed to specify a coordinate function:
    //          COEFF = 2 for edge parametrizations
    //          COEFF = 4 for triangular and quadrilateral face parametrizations
    // Arrays are sized and filled only when parametrization of a particular subcell is requested
    // by setSubcellParametrization.
    
    // Edge maps for 2D non-standard cells: ShellLine and Beam
    static FieldContainer<double> lineEdges;   static int lineEdgesSet = 0;
    
    // Edge maps for 2D standard cells: Triangle and Quadrilateral
    static FieldContainer<double> triEdges;         static int triEdgesSet  = 0;
    static FieldContainer<double> quadEdges;        static int quadEdgesSet = 0;

    // Edge maps for 3D non-standard cells: Shell Tri and Quad
    static FieldContainer<double> shellTriEdges;    static int shellTriEdgesSet  = 0;
    static FieldContainer<double> shellQuadEdges;   static int shellQuadEdgesSet = 0;
    
    // Edge maps for 3D standard cells:
    static FieldContainer<double> tetEdges;         static int tetEdgesSet = 0;
    static FieldContainer<double> hexEdges;         static int hexEdgesSet = 0;
    static FieldContainer<double> pyrEdges;         static int pyrEdgesSet = 0;
    static FieldContainer<double> wedgeEdges;       static int wedgeEdgesSet = 0;


    // Face maps for 3D non-standard cells: Shell Triangle and Quadrilateral
    static FieldContainer<double> shellTriFaces;    static int shellTriFacesSet  = 0;
    static FieldContainer<double> shellQuadFaces;   static int shellQuadFacesSet = 0;

    // Face maps for 3D standard cells:
    static FieldContainer<double> tetFaces;         static int tetFacesSet = 0;
    static FieldContainer<double> hexFaces;         static int hexFacesSet = 0;
    static FieldContainer<double> pyrFaces;         static int pyrFacesSet = 0;
    static FieldContainer<double> wedgeFaces;       static int wedgeFacesSet = 0;

    // Select subcell parametrization according to its parent cell type
    switch(parentCell.getKey() ) {
      
      // Tet cells
      case shards::Tetrahedron<4>::key:
      case shards::Tetrahedron<8>::key:
      case shards::Tetrahedron<10>::key:
        if(subcellDim == 2) {
          if(!tetFacesSet){
            setSubcellParametrization(tetFaces, subcellDim, parentCell);
            tetFacesSet = 1;
          }
          return tetFaces;
        }
        else if(subcellDim == 1) {
          if(!tetEdgesSet){
            setSubcellParametrization(tetEdges, subcellDim, parentCell);
            tetEdgesSet = 1;
          }
          return tetEdges;
        }
        else{
          TEST_FOR_EXCEPTION( (subcellDim != 1 || subcellDim != 2), std::invalid_argument, 
                              ">>> ERROR (Intrepid::CellTools::getSubcellParametrization): Tet parametrizations defined for 1 and 2-subcells only");
        }
        break;
        
      // Hex cells
      case shards::Hexahedron<8>::key:
      case shards::Hexahedron<20>::key:
      case shards::Hexahedron<27>::key:
        if(subcellDim == 2) {
          if(!hexFacesSet){
            setSubcellParametrization(hexFaces, subcellDim, parentCell);
            hexFacesSet = 1;
          }
          return hexFaces;
        }
        else if(subcellDim == 1) {
          if(!hexEdgesSet){
            setSubcellParametrization(hexEdges, subcellDim, parentCell);
            hexEdgesSet = 1;
          }
          return hexEdges;
        }
        else{
          TEST_FOR_EXCEPTION( (subcellDim != 1 || subcellDim != 2), std::invalid_argument, 
                              ">>> ERROR (Intrepid::CellTools::getSubcellParametrization): Hex parametrizations defined for 1 and 2-subcells only");
        }
        break;
        
      // Pyramid cells
      case shards::Pyramid<5>::key:
      case shards::Pyramid<13>::key:
      case shards::Pyramid<14>::key:
        if(subcellDim == 2) {
          if(!pyrFacesSet){
            setSubcellParametrization(pyrFaces, subcellDim, parentCell);
            pyrFacesSet = 1;
          }
          return pyrFaces;
        }
        else if(subcellDim == 1) {
          if(!pyrEdgesSet){
            setSubcellParametrization(pyrEdges, subcellDim, parentCell);
            pyrEdgesSet = 1;
          }
          return pyrEdges;
        }
        else {
          TEST_FOR_EXCEPTION( (subcellDim != 1 || subcellDim != 2), std::invalid_argument, 
                              ">>> ERROR (Intrepid::CellTools::getSubcellParametrization): Pyramid parametrizations defined for 1 and 2-subcells only");
        }
        break;
        
      // Wedge cells
      case shards::Wedge<6>::key:
      case shards::Wedge<15>::key:
      case shards::Wedge<18>::key:
        if(subcellDim == 2) {
          if(!wedgeFacesSet){
            setSubcellParametrization(wedgeFaces, subcellDim, parentCell);
            wedgeFacesSet = 1;
          }
          return wedgeFaces;
        }
        else if(subcellDim == 1) {
          if(!wedgeEdgesSet){
            setSubcellParametrization(wedgeEdges, subcellDim, parentCell);
            wedgeEdgesSet = 1;
          }
          return wedgeEdges;
        }
        else {
          TEST_FOR_EXCEPTION( (subcellDim != 1 || subcellDim != 2), std::invalid_argument, 
                              ">>> ERROR (Intrepid::CellTools::getSubcellParametrization): Wedge parametrization defined for 1 and 2-subcells only");
        }
        break;
      //
      // Standard 2D cells have only 1-subcells
      //
      case shards::Triangle<3>::key:
      case shards::Triangle<4>::key:
      case shards::Triangle<6>::key:
        if(subcellDim == 1) {
          if(!triEdgesSet){
            setSubcellParametrization(triEdges, subcellDim, parentCell);
            triEdgesSet = 1;
          }
          return triEdges;
        }
        else{
          TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                              ">>> ERROR (Intrepid::CellTools::getSubcellParametrization): Triangle parametrizations defined for 1-subcells only");
        }
        break;
                  
      case shards::Quadrilateral<4>::key:
      case shards::Quadrilateral<8>::key:
      case shards::Quadrilateral<9>::key:
        if(subcellDim == 1) {
          if(!quadEdgesSet){
            setSubcellParametrization(quadEdges, subcellDim, parentCell);
            quadEdgesSet = 1;
          }
          return quadEdges;
        }
        else{
          TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                              ">>> ERROR (Intrepid::CellTools::getSubcellParametrization): Quad parametrizations defined for 1-subcells only");
        }
        break;        
      //
      // Non-standard 3D Shell Tri and Quad cells have 1 and 2-subcells. Because they are 3D cells
      // can't reuse edge parametrization arrays for 2D Triangle and Quadrilateral.
      //
      case shards::ShellTriangle<3>::key:
      case shards::ShellTriangle<6>::key:
        if(subcellDim == 2) {
          if(!shellTriFacesSet){
            setSubcellParametrization(shellTriFaces, subcellDim, parentCell);
            shellTriFacesSet = 1;
          }
          return shellTriFaces;
        }
        else if(subcellDim == 1) {
          if(!shellTriEdgesSet){
            setSubcellParametrization(shellTriEdges, subcellDim, parentCell);
            shellTriEdgesSet = 1;
          }
          return shellTriEdges;
        }
        else if( subcellDim != 1 || subcellDim != 2){
          TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                              ">>> ERROR (Intrepid::CellTools::getSubcellParametrization): Shell Triangle parametrizations defined for 1 and 2-subcells only");
        }
        break;
        
      case shards::ShellQuadrilateral<4>::key:
      case shards::ShellQuadrilateral<8>::key:
      case shards::ShellQuadrilateral<9>::key:
        if(subcellDim == 2) {
          if(!shellQuadFacesSet){
            setSubcellParametrization(shellQuadFaces, subcellDim, parentCell);
            shellQuadFacesSet = 1;
          }
          return shellQuadFaces;
        }
        else if(subcellDim == 1) {
          if(!shellQuadEdgesSet){
            setSubcellParametrization(shellQuadEdges, subcellDim, parentCell);
            shellQuadEdgesSet = 1;
          }
          return shellQuadEdges;
        }
        else if( subcellDim != 1 || subcellDim != 2){
          TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                              ">>> ERROR (Intrepid::CellTools::getSubcellParametrization): Shell Quad parametrizations defined for 1 and 2-subcells only");
        }
        break;
        
      //
      // Non-standard 2D cells: Shell Lines and Beams have 1-subcells
      //
      case shards::ShellLine<2>::key:
      case shards::ShellLine<3>::key:
      case shards::Beam<2>::key:
      case shards::Beam<3>::key:
        if(subcellDim == 1) {
          if(!lineEdgesSet){
            setSubcellParametrization(lineEdges, subcellDim, parentCell);
            lineEdgesSet = 1;
          }
          return lineEdges;
        }
        else{
          TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                              ">>> ERROR (Intrepid::CellTools::getSubcellParametrization): shell line/beam parametrizations defined for 1-subcells only");
        }
        break;        

      default:
        TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                            ">>> ERROR (Intrepid::CellTools::getSubcellParametrization): invalid cell topology.");
    }//cell key
    // To disable compiler warning, should never be reached
    return lineEdges;
  }
  
  
  
  template<class Scalar>
  void CellTools<Scalar>::setSubcellParametrization(FieldContainer<double>&     subcellParametrization,
                                                    const int                   subcellDim,
                                                    const shards::CellTopology& parentCell)
  {
    // subcellParametrization is rank-3 FieldContainer with dimensions (SC, PCD, COEF) where:
    //  - SC    is the subcell count of subcells with the specified dimension in the parent cell
    //  - PCD   is Parent Cell Dimension, which gives the number of coordinate functions in the map
    //          PCD = 2 for standard 2D cells and non-standard 2D cells: shell line and beam
    //          PCD = 3 for standard 3D cells and non-standard 3D cells: shell Tri and Quad
    //  - COEF  is number of coefficients needed to specify a coordinate function:
    //          COEFF = 2 for edge parametrizations
    //          COEFF = 3 for both Quad and Tri face parametrizations. Because all Quad reference faces
    //          are affine, the coefficient of the bilinear term u*v is zero and is not stored, i.e.,
    //          3 coefficients are sufficient to store Quad face parameterization maps.
    //  
    // Edge parametrization maps [-1,1] to edge defined by (v0, v1)
    // Face parametrization maps [-1,1]^2 to quadrilateral face (v0, v1, v2, v3), or
    // standard 2-simplex  {(0,0),(1,0),(0,1)} to traingle face (v0, v1, v2).
    // This defines orientation-preserving parametrizations with respect to reference edge and
    // face orientations induced by their vertex order. 

    // get subcellParametrization dimensions: (sc, pcd, coeff)
    unsigned sc    = parentCell.getSubcellCount(subcellDim);
    unsigned pcd   = parentCell.getDimension();   
    unsigned coeff = (subcellDim == 1) ? 2 : 3;
    
    
    // Resize container
    subcellParametrization.resize(sc, pcd, coeff);
    
    // Edge parametrizations of 2D and 3D cells (shell lines and beams are 2D cells with edges)
    if(subcellDim == 1){      
      for(unsigned subcellOrd = 0; subcellOrd < sc; subcellOrd++){
        
        int v0ord = parentCell.getNodeMap(subcellDim, subcellOrd, 0);
        int v1ord = parentCell.getNodeMap(subcellDim, subcellOrd, 1);
        
        // vertexK[0] = x_k; vertexK[1] = y_k; vertexK[2] = z_k; z_k = 0 for 2D cells
        // Note that ShellLine and Beam are 2D cells!
        const double* v0 = getReferenceVertex(parentCell, v0ord);
        const double* v1 = getReferenceVertex(parentCell, v1ord);
        
        // x(t) = (x0 + x1)/2 + t*(x1 - x0)/2 
        subcellParametrization(subcellOrd, 0, 0) = (v0[0] + v1[0])/2.0;
        subcellParametrization(subcellOrd, 0, 1) = (v1[0] - v0[0])/2.0;
        
        // y(t) = (y0 + y1)/2 + t*(y1 - y0)/2 
        subcellParametrization(subcellOrd, 1, 0) = (v0[1] + v1[1])/2.0;
        subcellParametrization(subcellOrd, 1, 1) = (v1[1] - v0[1])/2.0;
        
        if( pcd == 3 ) {
          // z(t) = (z0 + z1)/2 + t*(z1 - z0)/2 
          subcellParametrization(subcellOrd, 2, 0) = (v0[2] + v1[2])/2.0;
          subcellParametrization(subcellOrd, 2, 1) = (v1[2] - v0[2])/2.0;
        }
      }// for loop over 1-subcells
    }
      
      // Face parametrizations of 2D and 3D cells: (shell Tri and Quad are 3D cells with faces)
      else if(subcellDim == 2) {
        for(unsigned subcellOrd = 0; subcellOrd < sc; subcellOrd++){
          
          switch(parentCell.getKey(subcellDim,subcellOrd)){
            
            // Admissible triangular faces for 3D cells in Shards:
            case shards::Triangle<3>::key:
            case shards::Triangle<4>::key:
            case shards::Triangle<6>::key: 
              {
                int v0ord = parentCell.getNodeMap(subcellDim, subcellOrd, 0);
                int v1ord = parentCell.getNodeMap(subcellDim, subcellOrd, 1);
                int v2ord = parentCell.getNodeMap(subcellDim, subcellOrd, 2);
                const double* v0 = getReferenceVertex(parentCell, v0ord);
                const double* v1 = getReferenceVertex(parentCell, v1ord);
                const double* v2 = getReferenceVertex(parentCell, v2ord);
                
                // x(u,v) = x0 + (x1 - x0)*u + (x2 - x0)*v
                subcellParametrization(subcellOrd, 0, 0) = v0[0];
                subcellParametrization(subcellOrd, 0, 1) = v1[0] - v0[0];
                subcellParametrization(subcellOrd, 0, 2) = v2[0] - v0[0];
                  
                // y(u,v) = y0 + (y1 - y0)*u + (y2 - y0)*v
                subcellParametrization(subcellOrd, 1, 0) = v0[1];
                subcellParametrization(subcellOrd, 1, 1) = v1[1] - v0[1];
                subcellParametrization(subcellOrd, 1, 2) = v2[1] - v0[1];
                
                // z(u,v) = z0 + (z1 - z0)*u + (z2 - z0)*v
                subcellParametrization(subcellOrd, 2, 0) = v0[2];
                subcellParametrization(subcellOrd, 2, 1) = v1[2] - v0[2];
                subcellParametrization(subcellOrd, 2, 2) = v2[2] - v0[2];
              }
              break;
              
            // Admissible quadrilateral faces for 3D cells in Shards:
            case shards::Quadrilateral<4>::key:
            case shards::Quadrilateral<8>::key:
            case shards::Quadrilateral<9>::key:
              {
                int v0ord = parentCell.getNodeMap(subcellDim, subcellOrd, 0);
                int v1ord = parentCell.getNodeMap(subcellDim, subcellOrd, 1);
                int v2ord = parentCell.getNodeMap(subcellDim, subcellOrd, 2);
                int v3ord = parentCell.getNodeMap(subcellDim, subcellOrd, 3);
                const double* v0 = getReferenceVertex(parentCell, v0ord);
                const double* v1 = getReferenceVertex(parentCell, v1ord);
                const double* v2 = getReferenceVertex(parentCell, v2ord);
                const double* v3 = getReferenceVertex(parentCell, v3ord);
                
                // x(u,v) = (x0+x1+x2+x3)/4+u*(-x0+x1+x2-x3)/4+v*(-x0-x1+x2+x3)/4+uv*(0=x0-x1+x2-x3)/4 
                subcellParametrization(subcellOrd, 0, 0) = ( v0[0] + v1[0] + v2[0] + v3[0])/4.0;
                subcellParametrization(subcellOrd, 0, 1) = (-v0[0] + v1[0] + v2[0] - v3[0])/4.0;
                subcellParametrization(subcellOrd, 0, 2) = (-v0[0] - v1[0] + v2[0] + v3[0])/4.0;
                
                // y(u,v) = (y0+y1+y2+y3)/4+u*(-y0+y1+y2-y3)/4+v*(-y0-y1+y2+y3)/4+uv*(0=y0-y1+y2-y3)/4 
                subcellParametrization(subcellOrd, 1, 0) = ( v0[1] + v1[1] + v2[1] + v3[1])/4.0;
                subcellParametrization(subcellOrd, 1, 1) = (-v0[1] + v1[1] + v2[1] - v3[1])/4.0;
                subcellParametrization(subcellOrd, 1, 2) = (-v0[1] - v1[1] + v2[1] + v3[1])/4.0;
                
                // z(u,v) = (z0+z1+z2+z3)/4+u*(-z0+z1+z2-z3)/4+v*(-z0-z1+z2+z3)/4+uv*(0=z0-z1+z2-z3)/4 
                subcellParametrization(subcellOrd, 2, 0) = ( v0[2] + v1[2] + v2[2] + v3[2])/4.0;
                subcellParametrization(subcellOrd, 2, 1) = (-v0[2] + v1[2] + v2[2] - v3[2])/4.0;
                subcellParametrization(subcellOrd, 2, 2) = (-v0[2] - v1[2] + v2[2] + v3[2])/4.0;
              }                
              break;
            default:
              TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                                  ">>> ERROR (Intrepid::CellTools::setSubcellParametrization): parametrization not defined for the specified face topology.");
             
          }// switch face topology key
        }// for subcellOrd
      }
  }
  
  
  
  template<class Scalar>
  const double* CellTools<Scalar>::getReferenceVertex(const shards::CellTopology& cell,
                                                      const int                   vertexOrd){
    
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( !(hasReferenceCell(cell) ), std::invalid_argument, 
                        ">>> ERROR (Intrepid::CellTools::getReferenceVertex): the specified cell does not have a reference cell.");
#endif
    
    // Cartesian coordinates of supported reference cell vertices, padded to three-dimensions.
    // Vertex order follows cell topology definition in Shards
    static const double line[2][3] =
      {
        {-1, 0, 0}, { 1, 0, 0} 
      };
    
    static const double triangle[3][3] =
      {
        { 0, 0, 0}, { 1, 0, 0}, { 0, 1, 0} 
      };
    
    static const double quadrilateral[4][3] =
      {
        {-1,-1, 0}, { 1,-1, 0}, { 1, 1, 0}, {-1, 1, 0}
      };
    
    static const double tetrahedron[4][3] =
      {
        { 0, 0, 0}, { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}
      };
    
    static const double hexahedron[8][3] = 
      {
        {-1, -1, -1}, {1, -1, -1}, {1,  1, -1}, {-1,  1, -1},
        {-1, -1,  1}, {1, -1,  1}, {1,  1,  1}, {-1,  1,  1}
      };
    
    static const double pyramid[5][3] =
      {
        {-1,-1, 0}, { 1,-1, 0}, { 1, 1, 0}, {-1, 1, 0}, { 0, 0, 1}
      };
    
    static const double wedge[6][3] =
      {
        { 0, 0,-1}, { 1, 0,-1}, { 0, 1,-1}, { 0, 0, 1}, { 1, 0, 1}, { 0, 1, 1} 
      };
    
    
    switch(cell.getKey() ) {
      
      // Line topologies
      case shards::Line<2>::key:
      case shards::Line<3>::key:
      case shards::ShellLine<2>::key:
      case shards::ShellLine<3>::key:
      case shards::Beam<2>::key:
      case shards::Beam<3>::key:
#ifdef HAVE_INTREPID_DEBUG
        TEST_FOR_EXCEPTION( !( (0 <= vertexOrd) && (vertexOrd < 2) ), std::invalid_argument,
                            ">>> ERROR (Intrepid::CellTools::getReferenceVertex): invalid vertex ordinal for a line. ");
#endif
        return line[vertexOrd];
        break;
        
      // Triangle topologies
      case shards::Triangle<3>::key:
      case shards::Triangle<4>::key:
      case shards::Triangle<6>::key:
      case shards::ShellTriangle<3>::key:
      case shards::ShellTriangle<6>::key:
#ifdef HAVE_INTREPID_DEBUG
        TEST_FOR_EXCEPTION( !( (0 <= vertexOrd) && (vertexOrd < 3) ), std::invalid_argument,
                            ">>> ERROR (Intrepid::CellTools::getReferenceVertex): invalid vertex ordinal for a triangle. ");
#endif
        return triangle[vertexOrd];
        break;
        
      // Quadrilateral topologies  
      case shards::Quadrilateral<4>::key:
      case shards::Quadrilateral<8>::key:
      case shards::Quadrilateral<9>::key:
      case shards::ShellQuadrilateral<4>::key:
      case shards::ShellQuadrilateral<8>::key:
      case shards::ShellQuadrilateral<9>::key:
#ifdef HAVE_INTREPID_DEBUG
        TEST_FOR_EXCEPTION( !( (0 <= vertexOrd) && (vertexOrd < 4) ), std::invalid_argument,
                            ">>> ERROR (Intrepid::CellTools::getReferenceVertex): invalid vertex ordinal for a quadrilateral. ");
#endif
        return quadrilateral[vertexOrd];
        break;
 
      // Tetrahedron topologies
      case shards::Tetrahedron<4>::key:
      case shards::Tetrahedron<8>::key:
      case shards::Tetrahedron<10>::key:
#ifdef HAVE_INTREPID_DEBUG
        TEST_FOR_EXCEPTION( !( (0 <= vertexOrd) && (vertexOrd < 4) ), std::invalid_argument,
                            ">>> ERROR (Intrepid::CellTools::getReferenceVertex): invalid vertex ordinal for a tetrahedron. ");
#endif
        return tetrahedron[vertexOrd];
        break;
        
      // Hexahedron topologies
      case shards::Hexahedron<8>::key:
      case shards::Hexahedron<20>::key:
      case shards::Hexahedron<27>::key:
#ifdef HAVE_INTREPID_DEBUG
        TEST_FOR_EXCEPTION( !( (0 <= vertexOrd) && (vertexOrd < 8) ), std::invalid_argument,
                            ">>> ERROR (Intrepid::CellTools::getReferenceVertex): invalid vertex ordinal for a hexahedron. ");
#endif
        return hexahedron[vertexOrd];
        break;
          
      // Pyramid topologies  
      case shards::Pyramid<5>::key:
      case shards::Pyramid<13>::key:
      case shards::Pyramid<14>::key:
#ifdef HAVE_INTREPID_DEBUG
        TEST_FOR_EXCEPTION( !( (0 <= vertexOrd) && (vertexOrd < 5) ), std::invalid_argument,
                            ">>> ERROR (Intrepid::CellTools::getReferenceVertex): invalid vertex ordinal for a pyramid. ");
#endif
        return pyramid[vertexOrd];
        break;
        
      // Wedge topologies
      case shards::Wedge<6>::key:
      case shards::Wedge<15>::key:
      case shards::Wedge<18>::key:
#ifdef HAVE_INTREPID_DEBUG
        TEST_FOR_EXCEPTION( !( (0 <= vertexOrd) && (vertexOrd < 6) ), std::invalid_argument,
                            ">>> ERROR (Intrepid::CellTools::getReferenceVertex): invalid vertex ordinal for a wedge. ");
#endif
        return wedge[vertexOrd];
        break;
        
      default:
        TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                            ">>> ERROR (Intrepid::CellTools::getReferenceVertex): invalid cell topology.");
    }
    // To disable compiler warning, should never be reached
    return line[0];
  }
    
  
  
  template<class Scalar>
  template<class ArrayOut>
  void CellTools<Scalar>::getReferenceSubcellVertices(ArrayOut&                   subcellVertices,
                                                      const int                   subcellDim,
                                                      const int                   subcellOrd,
                                                      const shards::CellTopology& parentCell){
      
    // Find how many vertices does the specified subcell have.
    int subcVertexCount = parentCell.getVertexCount(subcellDim, subcellOrd);
    
    // Loop over subcell vertices
    for(int subcVertOrd = 0; subcVertOrd < subcVertexCount; subcVertOrd++){
      
      // Get the vertex number relative to the parent reference cell
      int cellVertOrd = parentCell.getNodeMap(subcellDim, subcellOrd, subcVertOrd);
      
      // Loop over vertex Cartesian coordinates
      for(int dim = 0; dim < (int)parentCell.getDimension(); dim++){
        subcellVertices(subcVertOrd, dim) = CellTools::getReferenceVertex(parentCell, cellVertOrd)[dim];
      }
    }
  }  
  
  
  
  template<class Scalar>
  int CellTools<Scalar>::hasReferenceCell(const shards::CellTopology& cell) {
    
    switch(cell.getKey() ) {
      case shards::Line<2>::key:
      case shards::Line<3>::key:
      case shards::ShellLine<2>::key:
      case shards::ShellLine<3>::key:
      case shards::Beam<2>::key:
      case shards::Beam<3>::key:
        
      case shards::Triangle<3>::key:
      case shards::Triangle<4>::key:
      case shards::Triangle<6>::key:
      case shards::ShellTriangle<3>::key:
      case shards::ShellTriangle<6>::key:
        
      case shards::Quadrilateral<4>::key:
      case shards::Quadrilateral<8>::key:
      case shards::Quadrilateral<9>::key:
      case shards::ShellQuadrilateral<4>::key:
      case shards::ShellQuadrilateral<8>::key:
      case shards::ShellQuadrilateral<9>::key:
        
      case shards::Tetrahedron<4>::key:
      case shards::Tetrahedron<8>::key:
      case shards::Tetrahedron<10>::key:
        
      case shards::Hexahedron<8>::key:
      case shards::Hexahedron<20>::key:
      case shards::Hexahedron<27>::key:
        
      case shards::Pyramid<5>::key:
      case shards::Pyramid<13>::key:
      case shards::Pyramid<14>::key:
        
      case shards::Wedge<6>::key:
      case shards::Wedge<15>::key:
      case shards::Wedge<18>::key:
        return 1;
        break;
        
      default:
        return 0;
    }
    return 0;
  }
  
  //============================================================================================//
  //                                                                                            //
  //                     Jacobian, inverse Jacobian and Jacobian determinant                    //
  //                                                                                            //
  //============================================================================================//  
  
  template<class Scalar>
  template<class ArrayScalar>
  void CellTools<Scalar>::setJacobian(ArrayScalar &                jacobian,
                                      const ArrayScalar &          points,
                                      const ArrayScalar &          nodes,
                                      const shards::CellTopology & cellTopo,
                                      const int &                  whichCell) 
  {
    INTREPID_VALIDATE( setJacobianArgs(jacobian, points, nodes, whichCell,  cellTopo) );
    
    int spaceDim  = (int)cellTopo.getDimension();
    int numCells  = nodes.dimension(0);
    int numPoints = points.dimension(0);
    
    // Jacobian is computed using gradients of an appropriate H(grad) basis function: define RCP to the base class
    Teuchos::RCP< Basis< Scalar, FieldContainer<Scalar> > > HGRAD_Basis;
    
    // Choose the H(grad) basis depending on the cell topology. \todo Take into account extended, shell and beam cell topologies
    switch( cellTopo.getKey() ){
      
      // Standard Base topologies (number of nodes = number of vertices)
      case shards::Triangle<3>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TRI_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
      case shards::Quadrilateral<4>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_QUAD_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
      case shards::Tetrahedron<4>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
      case shards::Hexahedron<8>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_HEX_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
      case shards::Wedge<6>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_WEDGE_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      // Standard Extended topologies
      case shards::Triangle<6>::key:                        // curved triangle: use quadratic basis!
      case shards::Quadrilateral<8>::key:
      case shards::Quadrilateral<9>::key:
      case shards::Tetrahedron<10>::key:
      case shards::Hexahedron<20>::key:
      case shards::Hexahedron<27>::key:
      case shards::Wedge<15>::key:
      case shards::Wedge<18>::key:
        TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                            ">>> ERROR (Intrepid::CellTools::setJacobian): Cell topology not supported. ");
        break;
        
      // Base and Extended Line, Beam and Shell topologies  
      case shards::Line<2>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_LINE_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
      case shards::Line<3>::key:
      case shards::Beam<2>::key:
      case shards::Beam<3>::key:
      case shards::ShellLine<2>::key:
      case shards::ShellLine<3>::key:
      case shards::ShellTriangle<3>::key:
      case shards::ShellTriangle<6>::key:
      case shards::ShellQuadrilateral<4>::key:
      case shards::ShellQuadrilateral<8>::key:
      case shards::ShellQuadrilateral<9>::key:
        TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                            ">>> ERROR (Intrepid::CellTools::setJacobian): Cell topology not supported. ");
        break;
      default:
        TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                            ">>> ERROR (Intrepid::CellTools::setJacobian): Cell topology not supported.");        
    }  
    
    // Create an array for the values of basis functions gradients at the reference points
    int basisCardinality = HGRAD_Basis -> getCardinality();
    FieldContainer<Scalar> basisGrads(basisCardinality, numPoints, spaceDim);
    
    // Compute gradients of the basis functions
    HGRAD_Basis -> getValues(basisGrads, points, OPERATOR_GRAD);
    
    // Initialize jacobian
    for(int i = 0; i < jacobian.size(); i++){
      jacobian[i] = 0.0;
    }
    
    // The outer loops select the multi-index of the Jacobian entry: cell, point, row, col
    // If whichcell = -1, all jacobians are computed, otherwise a single cell jacobian is computed
    int cellLoop = (whichCell == -1) ? numCells : 1 ;
    
    for(int cellOrd = 0; cellOrd < cellLoop; cellOrd++) {
      for(int pointOrd = 0; pointOrd < numPoints; pointOrd++) {
        for(int row = 0; row < spaceDim; row++){
          for(int col = 0; col < spaceDim; col++){
            
            // The entry is computed by contracting the basis index. Number of basis functions and vertices must be the same
            for(int bfOrd = 0; bfOrd < basisCardinality; bfOrd++){
              
              if(whichCell == -1) {
                jacobian(cellOrd, pointOrd, row, col) += nodes(cellOrd, bfOrd, row)*basisGrads(bfOrd, pointOrd, col);
              }
              else {
                jacobian(pointOrd, row, col) += nodes(whichCell, bfOrd, row)*basisGrads(bfOrd, pointOrd, col);
              }
            } // bfOrd
            
          } // col
        } // row
      } // cellOrd
    } // pointOrd
  }
  


template<class Scalar>
template<class ArrayScalar>
void CellTools<Scalar>::setJacobianInv(ArrayScalar &        jacobianInv,
                                       const ArrayScalar &  jacobian) 
{
  INTREPID_VALIDATE( setJacobianInvArgs(jacobianInv, jacobian) );

  RealSpaceTools<Scalar>::inverse(jacobianInv, jacobian);
}



template<class Scalar>
template<class ArrayScalar>
void CellTools<Scalar>::setJacobianDet(ArrayScalar &         jacobianDet,
                                        const ArrayScalar &  jacobian)
{
  INTREPID_VALIDATE( setJacobianDetArgs(jacobianDet, jacobian) );

  RealSpaceTools<Scalar>::det(jacobianDet, jacobian);
}

//============================================================================================//
//                                                                                            //
//                      Reference-to-physical frame mapping and its inverse                   //
//                                                                                            //
//============================================================================================//

template<class Scalar>
template<class ArrayScalar>
void CellTools<Scalar>::mapToPhysicalFrame(ArrayScalar &                physPoints,
                                          const ArrayScalar &          refPoints,
                                          const ArrayScalar &          nodes,
                                          const shards::CellTopology & cellTopo,
                                          const int &                  whichCell)
{
  INTREPID_VALIDATE(mapToPhysicalFrameArgs( physPoints, refPoints, nodes, whichCell, cellTopo) );
  
  int spaceDim  = (int)cellTopo.getDimension();
  int numPoints = refPoints.dimension(0);
  int numCells  = nodes.dimension(0);
    
  // Mapping is computed using an appropriate H(grad) basis function: define RCP to the base class
  Teuchos::RCP<Basis<Scalar, FieldContainer<Scalar> > > HGRAD_Basis;
  
  // Choose the H(grad) basis depending on the cell topology. \todo take into account extended topologies
  switch( cellTopo.getKey() ){
    
    // Standard Base topologies (number of nodes = number of vertices)
    case shards::Triangle<3>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TRI_C1_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
    case shards::Quadrilateral<4>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_QUAD_C1_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
    case shards::Tetrahedron<4>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_C1_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
    case shards::Hexahedron<8>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_HEX_C1_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
    case shards::Wedge<6>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_WEDGE_C1_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
      
      // Standard Extended topologies
    case shards::Triangle<6>::key:                        // curved triangle: use quadratic basis!
    case shards::Quadrilateral<8>::key:
    case shards::Quadrilateral<9>::key:
    case shards::Tetrahedron<10>::key:
    case shards::Hexahedron<20>::key:
    case shards::Hexahedron<27>::key:
    case shards::Wedge<15>::key:
    case shards::Wedge<18>::key:
      TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                          ">>> ERROR (Intrepid::CellTools::setJacobian): Cell topology not supported. ");
      break;
      
      // Base and Extended Line, Beam and Shell topologies  
    case shards::Line<2>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_LINE_C1_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
    case shards::Line<3>::key:
    case shards::Beam<2>::key:
    case shards::Beam<3>::key:
    case shards::ShellLine<2>::key:
    case shards::ShellLine<3>::key:
    case shards::ShellTriangle<3>::key:
    case shards::ShellTriangle<6>::key:
    case shards::ShellQuadrilateral<4>::key:
    case shards::ShellQuadrilateral<8>::key:
    case shards::ShellQuadrilateral<9>::key:
      TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                          ">>> ERROR (Intrepid::CellTools::setJacobian): Cell topology not supported. ");
      break;
    default:
      TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                          ">>> ERROR (Intrepid::CellTools::setJacobian): Cell topology not supported.");        
  }// switch  

  // Create an array for the values of basis functions at the reference points
  int basisCardinality = HGRAD_Basis -> getCardinality();
  FieldContainer<Scalar> basisVals(basisCardinality, numPoints);
  
  // Compute basis values
  HGRAD_Basis -> getValues(basisVals, refPoints, OPERATOR_VALUE);
  
  // Initialize physPoints
  for(int i = 0; i < physPoints.size(); i++){
    physPoints[i] = 0.0;
  }
  
  // If whichcell = -1, ref pt. set is mapped to all cells, otherwise, the set is mapped to one cell only
  int cellLoop = (whichCell == -1) ? numCells : 1 ;

  // Compute the map F(refPoints) = sum node_coordinate*basis(refPoints)
  for(int cellOrd = 0; cellOrd < cellLoop; cellOrd++) {
    for(int pointOrd = 0; pointOrd < numPoints; pointOrd++) {
      for(int dim = 0; dim < spaceDim; dim++){
        for(int bfOrd = 0; bfOrd < basisCardinality; bfOrd++){
          
          if(whichCell == -1){
            physPoints(cellOrd, pointOrd, dim) += nodes(cellOrd, bfOrd, dim)*basisVals(bfOrd, pointOrd);
          }
          else{
            physPoints(pointOrd, dim) += nodes(whichCell, bfOrd, dim)*basisVals(bfOrd, pointOrd);
          }
        } // bfOrd
      }// dim
    }// cellOrd
  }//pointOrd
}

  

template<class Scalar>
template<class ArrayScalar>
void CellTools<Scalar>::mapToReferenceFrame(ArrayScalar &                 refPoints,
                                           const ArrayScalar &           physPoints,
                                           const ArrayScalar &           nodes,
                                           const shards::CellTopology &  cellTopo,
                                           const int &                   whichCell)
{
  INTREPID_VALIDATE( mapToReferenceFrameArgs(refPoints, physPoints, nodes, whichCell, cellTopo) );

  int spaceDim  = (int)cellTopo.getDimension();
  int numPoints = physPoints.dimension(0);
    
  // FieldContainer ctor initializes these arrays with zero  
  FieldContainer<Scalar> xOld(numPoints, spaceDim);
  FieldContainer<Scalar> xTem(numPoints, spaceDim);  
  FieldContainer<Scalar> jacobian(numPoints, spaceDim, spaceDim);
  FieldContainer<Scalar> jacobInv(numPoints, spaceDim, spaceDim);
  
  // Newton method to solve the equation F(refPoints) - physPoints = 0:
  // refPoints = xOld - DF^{-1}(xOld)*(F(xOld) - physPoints) = xOld + DF^{-1}(xOld)*(physPoints - F(xOld))
  for(int iter = 0; iter < INTREPID_MAX_NEWTON; ++iter)	{	
    
    // Jacobians at the old iterates and their inverses. 
    setJacobian(jacobian, xOld, nodes, cellTopo, whichCell);
    setJacobianInv(jacobInv, jacobian);
        
    // The Newton step.
    mapToPhysicalFrame( xTem, xOld, nodes, cellTopo, whichCell );      // xTem <- F(xOld)
    RealSpaceTools<Scalar>::subtract( xTem, physPoints, xTem );        // xTem <- physPoints - F(xOld)
    RealSpaceTools<Scalar>::matvec( refPoints, jacobInv, xTem);        // refPoints <- DF^{-1}( physPoints - F(xOld) )
    RealSpaceTools<Scalar>::add( refPoints, xOld );                    // refPoints <- DF^{-1}( physPoints - F(xOld) ) + xOld
    
    // l2 error (Euclidean distance) between old and new iterates: |xOld - xNew|
    FieldContainer<Scalar> error(numPoints); 
    RealSpaceTools<Scalar>::subtract( xTem, xOld, refPoints );
    RealSpaceTools<Scalar>::vectorNorm( error, xTem, NORM_TWO );
    
    // Average l2 error (sum of l2 errors in error divided by numPoints) used to check convergence
    double avError = RealSpaceTools<Scalar>::vectorNorm( error, NORM_ONE )/numPoints; 
    if (avError < INTREPID_TOL) {		          
      break;
    } 
    else if ( iter > INTREPID_MAX_NEWTON - 5) {
      INTREPID_VALIDATE(
      std::cout << " Intrepid::CellTools::mapToReferenceFrame failed to converge to desired tolerance within " 
                << INTREPID_MAX_NEWTON - 5 << " iterations\n" );
      break;
    }
    
    // initialize next Newton step
    xOld = refPoints;	
  } // for(iter)
}



template<class Scalar>
template<class ArrayTypeOut, class ArrayTypeIn>
void CellTools<Scalar>::mapToReferenceSubcell(ArrayTypeOut &                refSubcellPoints,
                                              const ArrayTypeIn &           paramPoints,
                                              const int                     subcellDim,
                                              const int                     subcellOrd,
                                              const shards::CellTopology &  parentCell){
  int cellDim = parentCell.getDimension();
  int numPts  = paramPoints.dimension(0);
    
  // Get the subcell map, i.e., the coefficients of the parametrization function for the subcell
  const FieldContainer<double>& subcellMap = getSubcellParametrization(subcellDim, subcellOrd, parentCell);

  // Apply the parametrization map to every point in parameter space
  if(subcellDim == 2) {
    for(int pt = 0; pt < numPts; pt++){
      double u = paramPoints(pt,0);
      double v = paramPoints(pt,1);
      
      // map_dim(u,v) = c_0(dim) + c_1(dim)*u + c_2(dim)*v because both Quad and Tri ref faces are aaffine!
      for(int  dim = 0; dim < cellDim; dim++){
        refSubcellPoints(pt, dim) = subcellMap(subcellOrd, dim, 0) + \
                                    subcellMap(subcellOrd, dim, 1)*u + \
                                    subcellMap(subcellOrd, dim, 2)*v;
      }
    }
  }
  else if(subcellDim == 1) {    
    for(int pt = 0; pt < numPts; pt++){
      for(int dim = 0; dim < cellDim; dim++) {
        refSubcellPoints(pt, dim) = subcellMap(subcellOrd, dim, 0) + subcellMap(subcellOrd, dim, 1)*paramPoints(pt,0);
      }
    }
  }
  else{
    TEST_FOR_EXCEPTION( !( (subcellDim == 1) || (subcellDim == 2) ), std::invalid_argument, 
                        ">>> ERROR (Intrepid::CellTools::mapToReferenceSubcell): defined only for 1 and 2-subcells");
  }
}


template<class Scalar>
template<class ArrayTypeOut>
void CellTools<Scalar>::getReferenceFaceTangents(ArrayTypeOut &                uTan,
                                                 ArrayTypeOut &                vTan,
                                                 const int &                   faceOrd,
                                                 const shards::CellTopology &  parentCell){
  int spaceDim  = parentCell.getDimension();
  
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( !(spaceDim == 3), std::invalid_argument, 
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceTangents): three-dimensional parent cell required");  
  
  TEST_FOR_EXCEPTION( !( (0 <= faceOrd) && (faceOrd < parentCell.getSubcellCount(2) ) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceTangents): face ordinal out of bounds");  

  TEST_FOR_EXCEPTION( !( (uTan.rank() == 1)  && (vTan.rank() == 1) ), std::invalid_argument,  
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceTangents): rank = 1 required for output arrays"); 
  
  TEST_FOR_EXCEPTION( !( uTan.dimension(0) == spaceDim ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceTangents): dim0 (spatial dim) must match that of parent cell");  

  TEST_FOR_EXCEPTION( !( vTan.dimension(0) == spaceDim ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceTangents): dim0 (spatial dim) must match that of parent cell");  
#endif
  
  // Face parametrizations are computed in setSubcellParametrization and stored in rank-3 array 
  // (subcOrd, coordinate, coefficient): retrieve this array
  const FieldContainer<double>& faceMap = getSubcellParametrization(2, faceOrd, parentCell);
  
  /*  All ref. face maps have affine coordinate functions:  f_dim(u,v) = C_0(dim) + C_1(dim)*u + C_2(dim)*v
   *                           `   => Tangent vectors are:  uTan -> C_1(*);    vTan -> C_2(*)
   */
    // set uTan -> C_1(*)
    uTan(0) = faceMap(faceOrd, 0, 1);
    uTan(1) = faceMap(faceOrd, 1, 1);
    uTan(2) = faceMap(faceOrd, 2, 1);
    
     // set vTan -> C_2(*)
    vTan(0) = faceMap(faceOrd, 0, 2);
    vTan(1) = faceMap(faceOrd, 1, 2);
    vTan(2) = faceMap(faceOrd, 2, 2);
}



template<class Scalar>
template<class ArrayTypeOut>
void CellTools<Scalar>::getReferenceFaceNormal(ArrayTypeOut &                normal,
                                               const int &                   faceOrd,
                                               const shards::CellTopology &  parentCell){
  int spaceDim  = parentCell.getDimension();
  
#ifdef HAVE_INTREPID_DEBUG
  
  TEST_FOR_EXCEPTION( !(spaceDim == 3), std::invalid_argument, 
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceNormal): three-dimensional parent cell required");  
  
  TEST_FOR_EXCEPTION( !( (0 <= faceOrd) && (faceOrd < parentCell.getSubcellCount(2) ) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceNormal): face ordinal out of bounds");  
  
  TEST_FOR_EXCEPTION( !( normal.rank() == 1 ), std::invalid_argument,  
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceNormal): rank = 1 required for output array"); 
    
  TEST_FOR_EXCEPTION( !( normal.dimension(0) == spaceDim ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceNormal): dim0 (spatial dim) must match that of parent cell");  
#endif

  // Reference face normal = vector product of reference face tangents. Allocate temp FC storage:
  FieldContainer<Scalar> uTan(spaceDim);
  FieldContainer<Scalar> vTan(spaceDim);
  getReferenceFaceTangents(uTan, vTan, faceOrd, parentCell);
  
  // Compute the vector product of the reference face tangents:
  RealSpaceTools<Scalar>::vecprod(normal, uTan, vTan);
}



template<class Scalar>
template<class ArrayTypeOut>
void CellTools<Scalar>::getReferenceEdgeTangent(ArrayTypeOut &                edgeTan,
                                                const int &                   edgeOrd,
                                                const shards::CellTopology &  parentCell){
  
  int spaceDim  = parentCell.getDimension();

#ifdef HAVE_INTREPID_DEBUG
  
  TEST_FOR_EXCEPTION( !( (spaceDim == 2) || (spaceDim == 3) ), std::invalid_argument, 
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceTangents): two or three-dimensional parent cell required");
  
  TEST_FOR_EXCEPTION( !( (0 <= edgeOrd) && (edgeOrd < parentCell.getSubcellCount(1) ) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceTangents): edge ordinal out of bounds");  
  
  TEST_FOR_EXCEPTION( !( edgeTan.size() == spaceDim ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceTangents): output array size is required to match space dimension");  
#endif
  
  // Edge parametrizations are computed in setSubcellParametrization and stored in rank-3 array 
  // (subcOrd, coordinate, coefficient)
  const FieldContainer<double>& edgeMap = getSubcellParametrization(1, edgeOrd, parentCell);
  
  // All ref. edge maps have affine coordinate functions: f_dim(u) = C_0(dim) + C_1(dim)*u, 
  //                                     => edge Tangent: -> C_1(*)
  edgeTan(0) = edgeMap(edgeOrd, 0, 1);
  edgeTan(1) = edgeMap(edgeOrd, 1, 1);
  
  // Skip last coordinate for 2D parent cells
  if(spaceDim == 3) {
    edgeTan(2) = edgeMap(edgeOrd, 2, 1);  
  }
}



template<class Scalar>
template<class ArrayTypeOut, class ArrayTypeIn>
void CellTools<Scalar>::getPhysicalFaceTangents(ArrayTypeOut &                uPhysTan,
                                                ArrayTypeOut &                vPhysTan,
                                                const ArrayTypeIn &           uvPoints,
                                                const ArrayTypeIn &           worksetJacobians,
                                                const int &                   worksetFaceOrd,
                                                const shards::CellTopology &  parentCell){
  int worksetSize = worksetJacobians.dimension(0);
  int facePtCount = worksetJacobians.dimension(1); 
  int pCellDim    = parentCell.getDimension();
  
#ifdef HAVE_INTREPID_DEBUG
  std::string errmsg = ">>> ERROR (Intrepid::CellTools::getPhysicalFaceTangents):";

  TEST_FOR_EXCEPTION( !(pCellDim == 3), std::invalid_argument, 
                      ">>> ERROR (Intrepid::CellTools::getPhysicalFaceTangents): three-dimensional parent cell required");  
  
  // (1) uPhysTan and vPhysTan are rank-3 (C,P,D) and D=3 is required
  TEST_FOR_EXCEPTION( !requireRankRange(errmsg, uPhysTan, 3,3), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, uPhysTan, 2, 3,3), std::invalid_argument, errmsg);

  TEST_FOR_EXCEPTION( !requireRankRange(errmsg, vPhysTan, 3,3), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, vPhysTan, 2, 3,3), std::invalid_argument, errmsg);

  TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, uPhysTan,  vPhysTan), std::invalid_argument, errmsg);      

  // (2) uvPoints is rank-2 (P,D) and D=2 is required (points are in 2D parametrization domain)
  TEST_FOR_EXCEPTION( !requireRankRange(errmsg, uvPoints, 2,2), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, uvPoints, 1, 2,2), std::invalid_argument, errmsg);

  // (3) worksetJacobians in rank-4 (C,P,D,D) and D=3 is required
  TEST_FOR_EXCEPTION( !requireRankRange(errmsg, worksetJacobians, 4,4), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 2, 3,3), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 3, 3,3), std::invalid_argument, errmsg);

  // (4) cross-check array dimensions: uPhysTan (C,P,D) vs. worksetJacobians (C,P,D,D) and uvPoints(P,2)
  TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, uPhysTan, 0,1,2,2,  worksetJacobians, 0,1,2,3), std::invalid_argument, errmsg);      
  TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, uPhysTan, 1,  uvPoints, 0), std::invalid_argument, errmsg);      

#endif
    
  // Temp storage for the pair of constant ref. face tangents: rank-1 (D) arrays
  FieldContainer<double> uRefTan(pCellDim);
  FieldContainer<double> vRefTan(pCellDim);
  getReferenceFaceTangents(uRefTan, vRefTan, worksetFaceOrd, parentCell);

  // Loop over workset faces and face points
  for(int pCell = 0; pCell < worksetSize; pCell++){
    for(int pt = 0; pt < facePtCount; pt++){
      
      // Apply parent cell Jacobian to ref. face tangents
      for(int dim = 0; dim < pCellDim; dim++){
        uPhysTan(pCell, pt, dim) = 0.0;
        vPhysTan(pCell, pt, dim) = 0.0;
        
        // Unroll loops: parent cell dimension can only be 3
        uPhysTan(pCell, pt, dim) = \
          worksetJacobians(pCell, pt, dim, 0)*uRefTan(0) + \
          worksetJacobians(pCell, pt, dim, 1)*uRefTan(1) + \
          worksetJacobians(pCell, pt, dim, 2)*uRefTan(2);
        vPhysTan(pCell, pt, dim) = \
          worksetJacobians(pCell, pt, dim, 0)*vRefTan(0) + \
          worksetJacobians(pCell, pt, dim, 1)*vRefTan(1) + \
          worksetJacobians(pCell, pt, dim, 2)*vRefTan(2);
      }// for dim
    }// for pt
  }// for pCell
}



template<class Scalar>
template<class ArrayTypeOut, class ArrayTypeIn>
void CellTools<Scalar>::getPhysicalFaceNormals(ArrayTypeOut &                faceNormals,
                                               const ArrayTypeIn &           uvPoints,
                                               const ArrayTypeIn &           worksetJacobians,
                                               const int &                   worksetFaceOrd,
                                               const shards::CellTopology &  parentCell){
  int worksetSize = worksetJacobians.dimension(0);
  int facePtCount = worksetJacobians.dimension(1); 
  int pCellDim    = parentCell.getDimension();
  
#ifdef HAVE_INTREPID_DEBUG
  std::string errmsg = ">>> ERROR (Intrepid::CellTools::getPhysicalFaceNormals):";
  
  TEST_FOR_EXCEPTION( !(pCellDim == 3), std::invalid_argument, 
                      ">>> ERROR (Intrepid::CellTools::getPhysicalFaceNormals): three-dimensional parent cell required");  
  
  // (1) faceNormals is rank-3 (C,P,D) and D=3 is required
  TEST_FOR_EXCEPTION( !requireRankRange(errmsg, faceNormals, 3,3), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, faceNormals, 2, 3,3), std::invalid_argument, errmsg);
    
  // (2) uvPoints is rank-2 (P,D) and D=2 is required (points are in 2D parametrization domain)
  TEST_FOR_EXCEPTION( !requireRankRange(errmsg, uvPoints, 2,2), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, uvPoints, 1, 2,2), std::invalid_argument, errmsg);
  
  // (3) worksetJacobians in rank-4 (C,P,D,D) and D=3 is required
  TEST_FOR_EXCEPTION( !requireRankRange(errmsg, worksetJacobians, 4,4), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 2, 3,3), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 3, 3,3), std::invalid_argument, errmsg);
  
  // (4) cross-check array dimensions: faceNormals (C,P,D) vs. worksetJacobians (C,P,D,D) and uvPoints(P,2)
  TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, faceNormals, 0,1,2,2,  worksetJacobians, 0,1,2,3), std::invalid_argument, errmsg);      
  TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, faceNormals, 1,  uvPoints, 0), std::invalid_argument, errmsg);      
  
#endif
  
  // Temp storage for physical face tangents: rank-3 (C,P,D) arrays
  FieldContainer<double> uPhysTan(worksetSize, facePtCount, pCellDim);
  FieldContainer<double> vPhysTan(worksetSize, facePtCount, pCellDim);
  getPhysicalFaceTangents(uPhysTan, vPhysTan, uvPoints, worksetJacobians, worksetFaceOrd, parentCell);
  
  // Compute the vector product of the physical face tangents:
  RealSpaceTools<Scalar>::vecprod(faceNormals, uPhysTan, vPhysTan);
  
  
}

//============================================================================================//
//                                                                                            //
//                                        Inclusion tests                                     //
//                                                                                            //
//============================================================================================//


template<class Scalar>
int CellTools<Scalar>::checkPointInclusion(const Scalar*                 point,
                                           const int                     pointDim,
                                           const shards::CellTopology &  cellTopo,
                                           const double &                threshold) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( !(pointDim == (int)cellTopo.getDimension() ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::checkPointInclusion): Point and cell dimensions do not match. ");
#endif
  int testResult = 1;
  
  // Using these values in the tests effectievly inflates the reference element to a larger one
  double minus_one = -1.0 - threshold;
  double plus_one  =  1.0 + threshold;
  double minus_zero = - threshold;
  
  // A cell with extended topology has the same reference cell as a cell with base topology. 
  // => testing for inclusion in a reference Triangle<> and a reference Triangle<6> relies on 
  // on the same set of inequalities. To eliminate unnecessary cases we switch on the base topology
  unsigned key = cellTopo.getBaseTopology() -> key ;
  switch( key ) {
    
    case shards::Line<>::key :
      if( !(minus_one <= point[0] && point[0] <= plus_one))  testResult = 0;
      break;
      
    case shards::Triangle<>::key : {
      Scalar distance = std::max( std::max( -point[0], -point[1] ), point[0] + point[1] - 1.0 );
      if( distance > threshold ) testResult = 0;
      break;
    }
      
    case shards::Quadrilateral<>::key :
      if(!( (minus_one <= point[0] && point[0] <= plus_one) && \
            (minus_one <= point[1] && point[1] <= plus_one) ) ) testResult = 0;   
      break;
      
    case shards::Tetrahedron<>::key : {
      Scalar distance = std::max(  std::max(-point[0],-point[1]), \
                                   std::max(-point[2], point[0] + point[1] + point[2] - 1)  );
      if( distance > threshold ) testResult = 0;
      break;
    }
      
    case shards::Hexahedron<>::key :
      if(!((minus_one <= point[0] && point[0] <= plus_one) && \
           (minus_one <= point[1] && point[1] <= plus_one) && \
           (minus_one <= point[2] && point[2] <= plus_one)))  \
             testResult = 0;
      break;
      
    // The base of the reference prism is the same as the reference triangle => apply triangle test
    // to X and Y coordinates and test whether Z is in [-1,1]
    case shards::Wedge<>::key : {
      Scalar distance = std::max( std::max( -point[0], -point[1] ), point[0] + point[1] - 1 );
      if( distance > threshold  || \
          point[2] < minus_one || point[2] > plus_one) \
            testResult = 0;
      break;
    }

    // The base of the reference pyramid is the same as the reference quad cell => a horizontal plane
    // through a point P(x,y,z) intersects the pyramid at a quadrilateral that equals the base quad 
    // scaled by (1-z) => P(x,y,z) is inside the pyramid <=> (x,y) is in [-1+z,1-z]^2 && 0 <= Z <= 1 
    case shards::Pyramid<>::key : {
      Scalar left  = minus_one + point[2];
      Scalar right = plus_one  - point[2];
      if(!( (left       <= point[0] && point[0] <= right) && \
            (left       <= point[1] && point[1] <= right) && 
            (minus_zero <= point[2] && point[2] <= plus_one) ) )  \
             testResult = 0;  
      break;
    }
      
    default:
      TEST_FOR_EXCEPTION( !( (key == shards::Line<>::key ) ||
                             (key == shards::Triangle<>::key)  ||
                             (key == shards::Quadrilateral<>::key) ||
                             (key == shards::Tetrahedron<>::key)  ||
                             (key == shards::Hexahedron<>::key)  ||
                             (key == shards::Wedge<>::key)  ||
                             (key == shards::Pyramid<>::key) ),
                          std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::checkPointInclusion): Invalid cell topology. ");
  }
  return testResult;
}



template<class Scalar>
template<class ArrayPoint>
int CellTools<Scalar>::checkPointsetInclusion(const ArrayPoint&             points,
                                              const shards::CellTopology &  cellTopo, 
                                              const double &                threshold) {
  
  // create temp output array depending on the rank of the input array 
  int rank = points.rank();
  FieldContainer<int> inRefCell;
  
  switch(rank) {
    case 1: inRefCell.resize(1); break;
    case 2: inRefCell.resize( points.dimension(0) ); break;
    case 3: inRefCell.resize( points.dimension(0), points.dimension(1) ); break;
  }

  // Call the inclusion method which returns inclusion results for all points
  checkPointwiseInclusion(inRefCell, points, cellTopo, threshold);
  
  // Check if any points were outside, break when finding the first one
  int allInside = 1;
  for(int i = 0; i < inRefCell.size(); i++ ){
    if (inRefCell[i] == 0) {
      allInside = 0;
      break;
    }
  }
   return allInside;
}



template<class Scalar>
template<class ArrayInt, class ArrayPoint>
void CellTools<Scalar>::checkPointwiseInclusion(ArrayInt &                    inRefCell,
                                                const ArrayPoint &            points,
                                                const shards::CellTopology &  cellTopo, 
                                                const double &                threshold) {
  // Initializations
  int apRank   = points.rank();
  int dim0     = 1;
  int dim1     = 1;
  int pointDim = 0;
  switch(apRank) {
    case 1:
      pointDim = points.dimension(0);
      break;
    case 2:
      dim1     = points.dimension(0);
      pointDim = points.dimension(1);
      break;
    case 3:
      dim0     = points.dimension(0);
      dim1     = points.dimension(1);
      pointDim = points.dimension(2);
      break;
    default:
      TEST_FOR_EXCEPTION( !( (1 <= apRank) && (apRank <= 3) ), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusion): rank-1, 2 or 3 required for input array. ");      
  }// switch
  
  // If rank of input array is 1,2, or 3, respectively; rank of output array must be 1, 1 and 2, resp.
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( !( (apRank == 1) && (inRefCell.rank() == 1 ) ||
                         (apRank == 2) && (inRefCell.rank() == 1 ) ||
                         (apRank == 3) && (inRefCell.rank() == 2 ) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusion): Output array has invalid rank. ");      
#endif
  
  // This method can handle up to rank-3 input arrays. The spatial dim must be the last dimension. 
  // The method uses [] accessor because array rank is determined at runtime and the appropriate
  // (i,j,..,k) accessor is not known. Use of [] requires the following offsets:
  //    for input array  = i0*dim1*pointDim + i1*dim1  (computed in 2 pieces: inPtr0 and inPtr1, resp)
  //    for output array = i0*dim1                     (computed in one piece: outPtr0)
  int inPtr0  = 0;
  int inPtr1  = 0;
  int outPtr0 = 0;
  Scalar point[3] = {0.0, 0.0, 0.0};
  
  for(int i0 = 0; i0 < dim0; i0++){
    outPtr0 = i0*dim1;
    inPtr0  = outPtr0*pointDim;
    
    for(int i1 = 0; i1 < dim1; i1++) {
      inPtr1 = inPtr0 + i1*pointDim;      
      point[0] = points[inPtr1];
      if(pointDim > 1) {
        point[1] = points[inPtr1 + 1];
        if(pointDim > 2) {
          point[2] = points[inPtr1 + 2];
          if(pointDim > 3) {
            TEST_FOR_EXCEPTION( !( (1 <= pointDim) && (pointDim <= 3)), std::invalid_argument, 
                                ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusion): Input array specifies invalid point dimension ");      
          }
        }
      } //if(pointDim > 1)
      inRefCell[outPtr0 + i1] = checkPointInclusion(point, pointDim, cellTopo, threshold);
    } // for (i1)
  } // for(i2)

}  


template<class Scalar>
template<class ArrayInt, class ArrayPoint, class ArrayScalar>
void CellTools<Scalar>::checkPointwiseInclusion(ArrayInt &                    inCell,
                                                const ArrayPoint &            points,
                                                const ArrayScalar &           nodes,
                                                const int &                   whichCell,
                                                const shards::CellTopology &  cell, 
                                                const double &                threshold)
{
  INTREPID_VALIDATE( checkPointwiseInclusionArgs(inCell, points, nodes, whichCell, cell) );
  
  // For cell topologies with reference cells this test maps the points back to the reference cell
  // and uses the method for reference cells
  unsigned baseKey = cell.getBaseTopology() -> key;
  
  switch(baseKey){
    
    case shards::Line<>::key :
    case shards::Triangle<>::key:
    case shards::Quadrilateral<>::key :
    case shards::Tetrahedron<>::key :
    case shards::Hexahedron<>::key :
    case shards::Wedge<>::key :
    case shards::Pyramid<>::key :
      {
        FieldContainer<Scalar> refPoints(points.dimension(0), points.dimension(1) );
        mapToReferenceFrame(refPoints, points, nodes, cell, whichCell);
        checkPointwiseInclusion(inCell, refPoints, cell, threshold );
        break;
      }
    default: 
      TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                          ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusion): cell topology not supported");
  }// switch
  
}


//============================================================================================//
//                                                                                            //
//                  Validation of input/output arguments for CellTools methods                //
//                                                                                            //
//============================================================================================//

template<class Scalar>
template<class ArrayScalar> 
void CellTools<Scalar>::setJacobianArgs(const ArrayScalar &          jacobian,
                                        const ArrayScalar &          points,
                                        const ArrayScalar &          nodes,
                                        const int &                  whichCell,
                                        const shards::CellTopology & cellTopo){
  
  // Validate nodes array
  TEST_FOR_EXCEPTION( (nodes.rank() != 3), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::setJacobianArgs): rank = 3 required for nodes array");
  
  TEST_FOR_EXCEPTION( (nodes.dimension(0) <= 0), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::setJacobianArgs): dim 0 (number of cells) >= 1 required for nodes array");
  
  TEST_FOR_EXCEPTION( (nodes.dimension(1) != (int)cellTopo.getSubcellCount(0) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::setJacobianArgs): dim 1 (number of nodes) of nodes array does not match cell topology");
  
  TEST_FOR_EXCEPTION( (nodes.dimension(2) != (int)cellTopo.getDimension() ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::setJacobianArgs): dim 2 (spatial dimension) of nodes array  does not match cell dimension");
  
  // Validate points array
  TEST_FOR_EXCEPTION( (points.rank() != 2), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::setJacobianArgs): rank = 2 required for points array");
  
  TEST_FOR_EXCEPTION( (points.dimension(0) <= 0), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::setJacobianArgs): dim 0 (number of points) >= 1 required for points array ");
  
  TEST_FOR_EXCEPTION( (points.dimension(1) != (int)cellTopo.getDimension() ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::setJacobianArgs): dim 1 (spatial dimension) of points array does not match cell dimension");

  
  // validate whichCell. It can be either -1 (default value) or a valid cell ordinal.
  TEST_FOR_EXCEPTION( !( ( (0 <= whichCell ) && (whichCell < nodes.dimension(0) ) || (whichCell == -1) ) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::setJacobianArgs): whichCell = -1 or a valid cell ordinal is required.");
  
  
  // Validate the output array for the Jacobian: if whichCell == -1 all Jacobians are computed
  if(whichCell == -1) {
    TEST_FOR_EXCEPTION( (jacobian.rank() != 4), std::invalid_argument, 
                        ">>> ERROR (Intrepid::CellTools::setJacobianArgs): rank = 4 required for jacobian array");
    
    TEST_FOR_EXCEPTION( (jacobian.dimension(0) != nodes.dimension(0)), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::setJacobianArgs): dim 0 (number of cells) of jacobian array must equal dim 0 of nodes array");
    
    TEST_FOR_EXCEPTION( (jacobian.dimension(1) != points.dimension(0)), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::setJacobianArgs): dim 1 (number of points) of jacobian array must equal dim 0 of points array");
    
    TEST_FOR_EXCEPTION( !(jacobian.dimension(2) == jacobian.dimension(3) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::setJacobianArgs): dim 2 = dim 3 (same spatial dimensions) required for jacobian array. ");
    
    TEST_FOR_EXCEPTION( !( (0 < jacobian.dimension(3) ) && (jacobian.dimension(3) < 4) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::setJacobianArgs): dim 2 and dim 3 (spatial dimensions) must be between 1 and 3. ");
  }     
  // A single cell Jacobian is computed when whichCell != -1 (whichCell value has been already validated, no checks here)
  else {
    TEST_FOR_EXCEPTION( (jacobian.rank() != 3), std::invalid_argument, 
                        ">>> ERROR (Intrepid::CellTools::setJacobianArgs): rank = 3 required for jacobian array");
    
    TEST_FOR_EXCEPTION( (jacobian.dimension(0) != points.dimension(0)), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::setJacobianArgs): dim 0 (number of points) of jacobian array must equal dim 0 of points array");
    
    TEST_FOR_EXCEPTION( !(jacobian.dimension(1) == jacobian.dimension(2) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::setJacobianArgs): dim 1 = dim 2 (same spatial dimensions) required for jacobian array. ");
    
    TEST_FOR_EXCEPTION( !( (0 < jacobian.dimension(1) ) && (jacobian.dimension(1) < 4) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::setJacobianArgs): dim 1 and dim 2 (spatial dimensions) must be between 1 and 3. ");
  }
}



template<class Scalar>
template<class ArrayScalar> 
void CellTools<Scalar>::setJacobianInvArgs(const ArrayScalar &          jacobianInv,
                                           const ArrayScalar &          jacobian)
{
  // Validate input jacobian array: admissible ranks & dimensions are: 
  // - rank-4 with dimensions (C,P,D,D), or rank-3 with dimensions (P,D,D).
  int jacobRank = jacobian.rank();
  TEST_FOR_EXCEPTION( !( (jacobRank == 4) || (jacobRank == 3) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::setJacobianInvArgs): rank = 4 or 3 required for jacobian array. ");
  
  // Verify correctness of spatial dimensions - they are the last two dimensions of the array: rank-2 and rank-1
  TEST_FOR_EXCEPTION( !(jacobian.dimension(jacobRank - 1) == jacobian.dimension(jacobRank - 2) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::setJacobianInvArgs): dim(rank-2) = dim(rank-2) (same spatial dimensions) required for jacobian array. ");
  
  TEST_FOR_EXCEPTION( !( (0 < jacobian.dimension(jacobRank - 1) ) && (jacobian.dimension(jacobRank - 1) < 4) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::setJacobianInvArgs): dim(rank-1) and dim(rank-2) (spatial dimensions) must be between 1 and 3. ");
  
  // Validate output jacobianInv array
  // Global function returns -1 if ranks don't match, the ordinal of first non-matching dimension + 1, or 0 if all dimensions match
  int result = compareArrays(jacobian, jacobianInv);
  
  TEST_FOR_EXCEPTION( result == -1, std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::setJacobianInvArgs): output jacobianInv and input jacobian arrays must have the same rank. ");
  
  if(result > 0) {
    std::ostringstream msg ;
    msg << ">>> ERROR (Intrepid::CellTools::setJacobianInvArgs): non-matching dimension " << result - 1 << " in jacobian and jacobianInv.";
    TEST_FOR_EXCEPTION( true, std::invalid_argument, msg);
  }  
}



template<class Scalar>
template<class ArrayScalar>
void CellTools<Scalar>::setJacobianDetArgs(const ArrayScalar &  jacobianDet,
                                           const ArrayScalar &  jacobian)
{
  // Validate input jacobian array: admissible ranks & dimensions are: 
  // - rank-4 with dimensions (C,P,D,D), or rank-3 with dimensions (P,D,D).
  int jacobRank = jacobian.rank();
  TEST_FOR_EXCEPTION( !( (jacobRank == 4) || (jacobRank == 3) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::setJacobianInvArgs): rank = 4 or 3 required for jacobian array. ");
  
  // Verify correctness of spatial dimensions - they are the last two dimensions of the array: rank-2 and rank-1
  TEST_FOR_EXCEPTION( !(jacobian.dimension(jacobRank - 1) == jacobian.dimension(jacobRank - 2) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::setJacobianInvArgs): dim(rank-2) = dim(rank-2) (same spatial dimensions) required for jacobian array. ");
  
  TEST_FOR_EXCEPTION( !( (0 < jacobian.dimension(jacobRank - 1) ) && (jacobian.dimension(jacobRank - 1) < 4) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::setJacobianInvArgs): dim(rank-1) and dim(rank-2) (spatial dimensions) must be between 1 and 3. ");

  
  // Validate output jacobianDet array: must be rank-2 with dimensions (C,P) if jacobian was rank-4:
  if(jacobRank == 4){
    TEST_FOR_EXCEPTION( !(jacobianDet.rank() == 2), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::setJacobianDetArgs): rank = 2 required for jacobianDet if jacobian is rank-4. ");
    
    TEST_FOR_EXCEPTION( !(jacobianDet.dimension(0) == jacobian.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::setJacobianDetArgs): dim 0 (number of cells) of jacobianDet array must equal dim 0 of jacobian array. ");
    
    TEST_FOR_EXCEPTION( !(jacobianDet.dimension(1) == jacobian.dimension(1) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::setJacobianDetArgs): dim 1 (number of points) of jacobianDet array must equal dim 1 of jacobian array.");  
  }
  
  // must be rank-1 with dimension (P) if jacobian was rank-3
  else {
    TEST_FOR_EXCEPTION( !(jacobianDet.rank() == 1), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::setJacobianDetArgs): rank = 1 required for jacobianDet if jacobian is rank-3. ");
    
    TEST_FOR_EXCEPTION( !(jacobianDet.dimension(0) == jacobian.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::setJacobianDetArgs): dim 0 (number of points) of jacobianDet array must equal dim 0 of jacobian array.");  
  }
}



template<class Scalar>
template<class ArrayScalar>
void CellTools<Scalar>::mapToPhysicalFrameArgs(const ArrayScalar &           physPoints,
                                              const ArrayScalar &           refPoints,
                                              const ArrayScalar &           nodes,
                                              const int&                    whichCell,
                                              const shards::CellTopology &  cellTopo)
{
  // Validate nodes array
  TEST_FOR_EXCEPTION( (nodes.rank() != 3), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrameArgs): rank = 3 required for nodes array");
  
  TEST_FOR_EXCEPTION( (nodes.dimension(0) <= 0), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrameArgs): dim 0 (number of cells) >= 1 required for nodes array");
  
  TEST_FOR_EXCEPTION( (nodes.dimension(1) != (int)cellTopo.getSubcellCount(0) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrameArgs): dim 1 (number of nodes) of nodes array does not match cell topology");
  
  TEST_FOR_EXCEPTION( (nodes.dimension(2) != (int)cellTopo.getDimension() ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrameArgs): dim 2 (spatial dimension) of nodes array  does not match cell dimension");
  
  
  // Validate refPoints array
  TEST_FOR_EXCEPTION( (refPoints.rank() != 2), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrameArgs): rank = 2 required for refPoints array");
  
  TEST_FOR_EXCEPTION( (refPoints.dimension(0) <= 0), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrameArgs): dim 0 (number of points) >= 1 required for refPoints array ");
  
  TEST_FOR_EXCEPTION( (refPoints.dimension(1) != (int)cellTopo.getDimension() ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrameArgs): dim 1 (spatial dimension) of refPoints array does not match cell dimension");

  
  // validate whichCell. It can be either -1 (default value) or a valid cell ordinal.
  TEST_FOR_EXCEPTION( !( ( (0 <= whichCell ) && (whichCell < nodes.dimension(0) ) || (whichCell == -1) ) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::setJacobianArgs): whichCell = -1 or a valid cell ordinal is required.");
  
  
  // Validate output array: whichCell = -1 requires rank-3 array with dimensions (C,P,D)  
  if(whichCell == -1) {
    TEST_FOR_EXCEPTION( (physPoints.rank() != 3), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrameArgs): rank = 3 required for physPoints array");
    
    TEST_FOR_EXCEPTION( (physPoints.dimension(0) != nodes.dimension(0)), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrameArgs): dim 0 (number of cells) of physPoints array must equal dim 0 of nodes array");
    
    TEST_FOR_EXCEPTION( (physPoints.dimension(1) != refPoints.dimension(0)), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrameArgs): dim 1 (number of points) of physPoints array must equal dim 0 of refPoints array"); 
    
    TEST_FOR_EXCEPTION( (physPoints.dimension(2) != (int)cellTopo.getDimension()), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrameArgs): dim 2 (spatial dimension) does not match cell dimension ");  
  }
  // 0 <= whichCell < num cells requires rank-2 array (P,D)
  else{
    TEST_FOR_EXCEPTION( (physPoints.rank() != 2), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrameArgs): rank = 2 required for physPoints array");
    
    TEST_FOR_EXCEPTION( (physPoints.dimension(0) != refPoints.dimension(0)), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrameArgs): dim 0 (number of points) of physPoints array must equal dim 0 of refPoints array"); 
    
    TEST_FOR_EXCEPTION( (physPoints.dimension(1) != (int)cellTopo.getDimension()), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrameArgs): dim 1 (spatial dimension) does not match cell dimension ");      
  }
}



template<class Scalar>
template<class ArrayScalar>
void CellTools<Scalar>::mapToReferenceFrameArgs(const ArrayScalar &          refPoints,
                                               const ArrayScalar &           physPoints,
                                               const ArrayScalar &           nodes,
                                               const int &                   whichCell,
                                               const shards::CellTopology &  cellTopo)
{
  // Validate nodes array
  TEST_FOR_EXCEPTION( (nodes.rank() != 3), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToReferenceFrameArgs): rank = 3 required for nodes array");
  
  TEST_FOR_EXCEPTION( (nodes.dimension(0) <= 0), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToReferenceFrameArgs): dim 0 (number of cells) >= 1 required for nodes array");
  
  TEST_FOR_EXCEPTION( (nodes.dimension(1) != (int)cellTopo.getSubcellCount(0) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToReferenceFrameArgs): dim 1 (number of nodes) of nodes array does not match cell topology");
  
  TEST_FOR_EXCEPTION( (nodes.dimension(2) != (int)cellTopo.getDimension() ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToReferenceFrameArgs): dim 2 (spatial dimension) of nodes array  does not match cell dimension");
  
  
  // Validate whichCell
  TEST_FOR_EXCEPTION( !( (0 <= whichCell ) && (whichCell < nodes.dimension(0) ) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToReferenceFrameArgs): whichCell argument is out of bounds.");
  
  
  // Validate physPoints array
  TEST_FOR_EXCEPTION( (physPoints.rank() != 2), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToReferenceFrameArgs): rank = 2 required for physPoints array");
  
  TEST_FOR_EXCEPTION( (physPoints.dimension(0) <= 0), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToReferenceFrameArgs): dim 0 (number of points) >= 1 required for physPoints array ");
  
  TEST_FOR_EXCEPTION( (physPoints.dimension(1) != (int)cellTopo.getDimension() ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToReferenceFrameArgs): dim 1 (spatial dimension) of physPoints array does not match cell dimension");
  
  

  // Validate refPoints array: refPoints and physPoints must have the same rank and matching dimensions
  int result = compareArrays(physPoints, refPoints);
  
  TEST_FOR_EXCEPTION( result == -1, std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToReferenceFrameArgs): physPoints and refPoints have non-matching ranks. ");
  
  if(result > 0) {
    std::ostringstream msg ;
    msg << ">>> ERROR (Intrepid::CellTools::mapToReferenceFrameArgs): non-matching dimension " << result - 1 << " in physPoints and refPoints.";
    TEST_FOR_EXCEPTION( true, std::invalid_argument, msg);
  }    
}



template<class Scalar>
template<class ArrayInt, class ArrayPoint, class ArrayScalar>
void CellTools<Scalar>::checkPointwiseInclusionArgs(ArrayInt &                    inCell,
                                                    const ArrayPoint &            physPoints,
                                                    const ArrayScalar &           nodes,
                                                    const int &                   whichCell,
                                                    const shards::CellTopology &  cell)
{
  // Validate nodes array
  TEST_FOR_EXCEPTION( (nodes.rank() != 3), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusionArgs): rank = 3 required for nodes array");
  
  TEST_FOR_EXCEPTION( (nodes.dimension(0) <= 0), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusionArgs): dim 0 (number of cells) >= 1 required for nodes array");
  
  TEST_FOR_EXCEPTION( (nodes.dimension(1) != (int)cell.getSubcellCount(0) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusionArgs): dim 1 (number of nodes) of nodes array does not match cell topology");
  
  TEST_FOR_EXCEPTION( (nodes.dimension(2) != (int)cell.getDimension() ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusionArgs): dim 2 (spatial dimension) of nodes array  does not match cell dimension");
  
  
  // Validate whichCell
  TEST_FOR_EXCEPTION( !( (0 <= whichCell ) && (whichCell < nodes.dimension(0) ) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusionArgs): whichCell argument is out of bounds.");
  
  
  // Validate physPoints array
  TEST_FOR_EXCEPTION( (physPoints.rank() != 2), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusionArgs): rank = 2 required for physPoints array");
  
  TEST_FOR_EXCEPTION( (physPoints.dimension(0) <= 0), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusionArgs): dim 0 (number of points) >= 1 required for physPoints array ");
  
  TEST_FOR_EXCEPTION( (physPoints.dimension(1) != (int)cell.getDimension() ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusionArgs): dim 1 (spatial dimension) of physPoints array does not match cell dimension");
  
  // Validate inCell
  TEST_FOR_EXCEPTION( (inCell.rank() != 1), std::invalid_argument, 
                      ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusionArgs): rank = 1 required for inCell array");
  
  TEST_FOR_EXCEPTION( (inCell.dimension(0) != physPoints.dimension(0)), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusionArgs): dim 0 (number of points) of inCell array must equal dim 0 of physPoints array");
}



//============================================================================================//
//                                                                                            //
//                                           Debug                                            //
//                                                                                            //
//============================================================================================//


template<class Scalar>
void CellTools<Scalar>::printSubcellVertices(const int subcellDim,
                                             const int subcellOrd,
                                             const shards::CellTopology & parentCell){
  
  // Get number of vertices for the specified subcell and parent cell dimension
  int subcVertexCount = parentCell.getVertexCount(subcellDim, subcellOrd);
  int cellDim         = parentCell.getDimension();
  
  // Allocate space for the subcell vertex coordinates
  FieldContainer<double> subcellVertices(subcVertexCount, cellDim);
  
  // Retrieve the vertex coordinates
  getReferenceSubcellVertices(subcellVertices,
                              subcellDim,
                              subcellOrd,
                              parentCell);
  
  // Print the vertices
  std::cout 
    << " Subcell " << std::setw(2) << subcellOrd 
    <<  " is " << parentCell.getName(subcellDim, subcellOrd) << " with vertices = {";
  
  // Loop over subcell vertices
  for(int subcVertOrd = 0; subcVertOrd < subcVertexCount; subcVertOrd++){
    std::cout<< "(";
    
    // Loop over vertex Cartesian coordinates
    for(int dim = 0; dim < (int)parentCell.getDimension(); dim++){
      std::cout << subcellVertices(subcVertOrd, dim);
      if(dim < (int)parentCell.getDimension()-1 ) { std::cout << ","; }
    }
    std::cout<< ")";
    if(subcVertOrd < subcVertexCount - 1) { std::cout << ", "; }
  }
  std::cout << "}\n";
}
  

template<class Scalar>
template<class ArrayTypeIn>
void CellTools<Scalar>::printWorksetSubcell(const ArrayTypeIn&            worksetNodes,
                                            const shards::CellTopology&   parentCell,
                                            const int&                    pCellOrd,
                                            const int&                    subcellDim,
                                            const int&                    subcellOrd,
                                            const int&                    fieldWidth){
  
  // Get the ordinals, relative to reference cell, of subcell nodes
  int subcNodeCount = parentCell.getNodeCount(subcellDim, subcellOrd);
  int pCellDim      = parentCell.getDimension();
  std::vector<int> subcNodeOrdinals(subcNodeCount);
  
  for(int i = 0; i < subcNodeCount; i++){
    subcNodeOrdinals[i] = parentCell.getNodeMap(subcellDim, subcellOrd, i);
  }
  
  // Loop over parent cells and print subcell nodes
  
  std::cout 
    << " Subcell " << subcellOrd << " on parent cell " << pCellOrd << " is " 
    << parentCell.getName(subcellDim, subcellOrd) << " with node(s) \n ({";
  
  for(int i = 0; i < subcNodeCount; i++){
    
    // print Cartesian coordinates of the node
    for(int dim = 0; dim < pCellDim; dim++){
      std::cout
      << std::setw(fieldWidth) << std::right << worksetNodes(pCellOrd, subcNodeOrdinals[i], dim); 
      if(dim < pCellDim - 1){ std::cout << ","; }
    }
    std::cout << "}";
    if(i < subcNodeCount - 1){ std::cout <<", {"; }
  }
  std::cout << ")\n\n";
}



} // namespace Intrepid
