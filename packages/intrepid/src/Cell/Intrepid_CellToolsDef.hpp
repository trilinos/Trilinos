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


/** \file   Intrepid_CellToolsDef.hpp
    \brief  Definition file for the Intrepid::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
*/
#ifndef INTREPID_CELLTOOLSDEF_HPP
#define INTREPID_CELLTOOLSDEF_HPP


namespace Intrepid {

  template<class Scalar>
  const FieldContainer<double>& CellTools<Scalar>::getSubcellParametrization(const int subcellDim,
                                                                             const shards::CellTopology& parentCell){
    
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( !(hasReferenceCell(parentCell) ), std::invalid_argument, 
                        ">>> ERROR (Intrepid::CellTools::getSubcellParametrization): the specified cell topology does not have a reference cell.");

    TEST_FOR_EXCEPTION( !( (1 <= subcellDim) && (subcellDim <= 2 ) ), std::invalid_argument,
                           ">>> ERROR (Intrepid::CellTools::getSubcellParametrization): parametrization defined only for 1 and 2-dimensional subcells.");    
#endif
    
    // Coefficients of the coordinate functions defining the parametrization maps are stored in 
    // rank-3 arrays with dimensions (SC, PCD, COEF) where:
    //  - SC    is the subcell count of subcells with the specified dimension in the parent cell
    //  - PCD   is Parent Cell Dimension, which gives the number of coordinate functions in the map:
    //          PCD = 2 for standard 2D cells and non-standard 2D cells: shell line and beam
    //          PCD = 3 for standard 3D cells and non-standard 3D cells: shell Tri and Quad
    //  - COEF  is number of coefficients needed to specify a coordinate function:
    //          COEFF = 2 for edge parametrizations
    //          COEFF = 3 for both Quad and Tri face parametrizations. Because all Quad reference faces
    //          are affine, the coefficient of the bilinear term u*v is zero and is not stored, i.e.,
    //          3 coefficients are sufficient to store Quad face parameterization maps.
    //
    // Arrays are sized and filled only when parametrization of a particular subcell is requested
    // by setSubcellParametrization.
    
    // Edge maps for 2D non-standard cells: ShellLine and Beam
    static FieldContainer<double> lineEdges;        static int lineEdgesSet = 0;
    
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
      case shards::Tetrahedron<11>::key:
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
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( !(hasReferenceCell(parentCell) ), std::invalid_argument, 
                        ">>> ERROR (Intrepid::CellTools::setSubcellParametrization): the specified cell topology does not have a reference cell.");

    TEST_FOR_EXCEPTION( !( (1 <= subcellDim) && (subcellDim <= 2 ) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::setSubcellParametrization): parametrization defined only for 1 and 2-dimensional subcells.");    
#endif
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
      
      // Face parametrizations of 3D cells: (shell Tri and Quad are 3D cells with faces)
      // A 3D cell can have both Tri and Quad faces, but because they are affine images of the
      // parametrization domain, 3 coefficients are enough to store them in both cases.
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
                        ">>> ERROR (Intrepid::CellTools::getReferenceVertex): the specified cell topology does not have a reference cell.");
    
    TEST_FOR_EXCEPTION( !( (0 <= vertexOrd) && (vertexOrd < (int)cell.getVertexCount() ) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::getReferenceVertex): invalid node ordinal for the specified cell topology. ");
#endif
    
    // Simply call getReferenceNode with the base topology of the cell
    return getReferenceNode(cell.getBaseCellTopologyData(), vertexOrd);
  }
    
  
  
  template<class Scalar>
  template<class ArraySubcellVert>
  void CellTools<Scalar>::getReferenceSubcellVertices(ArraySubcellVert &          subcellVertices,
                                                      const int                   subcellDim,
                                                      const int                   subcellOrd,
                                                      const shards::CellTopology& parentCell){
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( !(hasReferenceCell(parentCell) ), std::invalid_argument, 
                        ">>> ERROR (Intrepid::CellTools::getReferenceSubcellVertices): the specified cell topology does not have a reference cell.");

    // subcellDim can equal the cell dimension because the cell itself is a valid subcell! In this case
    // the method will return all cell cellWorkset.
    TEST_FOR_EXCEPTION( !( (0 <= subcellDim) && (subcellDim <= (int)parentCell.getDimension()) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::getReferenceSubcellVertices): subcell dimension out of range.");
    
    TEST_FOR_EXCEPTION( !( (0 <= subcellOrd) && (subcellOrd < (int)parentCell.getSubcellCount(subcellDim) ) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::getReferenceSubcellVertices): subcell ordinal out of range.");
        
    // Verify subcellVertices rank and dimensions
    {
      std::string errmsg = ">>> ERROR (Intrepid::CellTools::getReferenceSubcellVertices):";
      TEST_FOR_EXCEPTION( !( requireRankRange(errmsg, subcellVertices, 2, 2) ), std::invalid_argument, errmsg);
      
      int subcVertexCount = parentCell.getVertexCount(subcellDim, subcellOrd);
      int spaceDim = parentCell.getDimension();
        
      TEST_FOR_EXCEPTION( !( requireDimensionRange(errmsg, subcellVertices, 0,  subcVertexCount, subcVertexCount) ),
                          std::invalid_argument, errmsg);
      
      TEST_FOR_EXCEPTION( !( requireDimensionRange(errmsg, subcellVertices, 1,  spaceDim, spaceDim) ),
                          std::invalid_argument, errmsg);
    }
#endif 
    
    // Simply call getReferenceNodes with the base topology
    getReferenceSubcellNodes(subcellVertices, subcellDim, subcellOrd, parentCell.getBaseCellTopologyData() );
  }  

  
  
  template<class Scalar>
  const double* CellTools<Scalar>::getReferenceNode(const shards::CellTopology& cell,
                                                    const int                   nodeOrd){
    
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( !(hasReferenceCell(cell) ), std::invalid_argument, 
                        ">>> ERROR (Intrepid::CellTools::getReferenceNode): the specified cell topology does not have a reference cell.");

    TEST_FOR_EXCEPTION( !( (0 <= nodeOrd) && (nodeOrd < (int)cell.getNodeCount() ) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::getReferenceNode): invalid node ordinal for the specified cell topology. ");
#endif
    
    // Cartesian coordinates of supported reference cell cellWorkset, padded to three-dimensions.
    // Node order follows cell topology definition in Shards
    static const double line[2][3] ={
      {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0} 
    };
    static const double line_3[3][3] = {
      {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0},     
      // Extension node: edge midpoint
      { 0.0, 0.0, 0.0}
    };
    
    
    // Triangle topologies
    static const double triangle[3][3] = {
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0} 
    };
    static const double triangle_4[4][3] = {
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
      // Extension node: cell center
      { 1/3, 1/3, 0.0}
    };
    static const double triangle_6[6][3] = {
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
      // Extension cellWorkset: 3 edge midpoints
      { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}
    };
    
    
    // Quadrilateral topologies
    static const double quadrilateral[4][3] = {
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}
    };
    static const double quadrilateral_8[8][3] = {
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
      // Extension cellWorkset: 4 edge midpoints
      { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}
    };
    static const double quadrilateral_9[9][3] = {
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
      // Extension cellWorkset: 4 edge midpoints + 1 cell center
      { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}, { 0.0, 0.0, 0.0}
    };
    
    
    // Tetrahedron topologies
    static const double tetrahedron[4][3] = {
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0}
    };
    static const double tetrahedron_8[8][3] = {
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
      // Extension cellWorkset: 4 face centers (do not follow natural face order - see the cell topology!)
      { 1/3, 0.0, 1/3}, { 1/3, 1/3, 1/3}, { 1/3, 1/3, 0.0}, { 0.0, 1/3, 1/3} 
    };
    static const double tetrahedron_10[10][3] = {
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
      // Extension cellWorkset: 6 edge midpoints
      { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}, { 0.0, 0.0, 0.5}, { 0.5, 0.0, 0.5}, { 0.0, 0.5, 0.5}
    };

    static const double tetrahedron_11[10][3] = {
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
      // Extension cellWorkset: 6 edge midpoints
      { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}, { 0.0, 0.0, 0.5}, { 0.5, 0.0, 0.5}, { 0.0, 0.5, 0.5}
    };

    
    // Hexahedron topologies
    static const double hexahedron[8][3] = {
      {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
      {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0}
    };
    static const double hexahedron_20[20][3] = {
      {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
      {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0},
      // Extension cellWorkset: 12 edge midpoints (do not follow natural edge order - see cell topology!)
      { 0.0,-1.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, {-1.0, 0.0,-1.0}, 
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
      { 0.0,-1.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0}, {-1.0, 0.0, 1.0}
    };
    static const double hexahedron_27[27][3] = {
      {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
      {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0},
      // Extension cellWorkset: 12 edge midpoints + 1 cell center + 6 face centers  (do not follow natural subcell order!)
      { 0.0,-1.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, {-1.0, 0.0,-1.0}, 
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
      { 0.0,-1.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0}, {-1.0, 0.0, 1.0},
      { 0.0, 0.0, 0.0},
      { 0.0, 0.0,-1.0}, { 0.0, 0.0, 1.0}, {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, {0.0,-1.0, 0.0}, {0.0, 1.0, 0.0} 
    };
    
    
    // Pyramid topologies
    static const double pyramid[5][3] = {
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0}
    };
    static const double pyramid_13[13][3] = {
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
      // Extension cellWorkset: 8 edge midpoints 
      { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0},
      {-0.5,-0.5, 0.5}, { 0.5,-0.5, 0.5}, { 0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}   
    };
    static const double pyramid_14[14][3] = {
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
      // Extension cellWorkset: 8 edge midpoints + quadrilateral face midpoint 
      { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0},
      {-0.5,-0.5, 0.5}, { 0.5,-0.5, 0.5}, { 0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}, { 0.0, 0.0, 0.0}  
    };
    
    
    // Wedge topologies
    static const double wedge[6][3] = {
      { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0} 
    };
    static const double wedge_15[15][3] = {
      { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0},
      // Extension cellWorkset: 9 edge midpoints (do not follow natural edge order - see cell topology!)
      { 0.5, 0.0,-1.0}, { 0.5, 0.5,-1.0}, { 0.0, 0.5,-1.0}, { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
      { 0.5, 0.0, 1.0}, { 0.5, 0.5, 1.0}, { 0.0, 0.5, 1.0}
    };
    static const double wedge_18[18][3] = {
      { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0},
      // Extension cellWorkset: 9 edge midpoints + 3 quad face centers (do not follow natural subcell order - see cell topology!)
      { 0.5, 0.0,-1.0}, { 0.5, 0.5,-1.0}, { 0.0, 0.5,-1.0}, { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
      { 0.5, 0.0, 1.0}, { 0.5, 0.5, 1.0}, { 0.0, 0.5, 1.0},
      { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}
    };
    
    
    switch(cell.getKey() ) {
      
      // Base line topologies
      case shards::Line<2>::key:
      case shards::ShellLine<2>::key:
      case shards::Beam<2>::key:
        return line[nodeOrd];
        break;
        
      // Extended line topologies
      case shards::Line<3>::key:
      case shards::ShellLine<3>::key:
      case shards::Beam<3>::key:
        return line_3[nodeOrd];
        break;
        
        
      // Base triangle topologies
      case shards::Triangle<3>::key:
      case shards::ShellTriangle<3>::key:
        return triangle[nodeOrd];
        break;
        
      // Extened Triangle topologies
      case shards::Triangle<4>::key:
        return triangle_4[nodeOrd];
        break;
      case shards::Triangle<6>::key:
      case shards::ShellTriangle<6>::key:
        return triangle_6[nodeOrd];
        break;
        
        
      // Base Quadrilateral topologies  
      case shards::Quadrilateral<4>::key:
      case shards::ShellQuadrilateral<4>::key:
        return quadrilateral[nodeOrd];
        break;
        
      // Extended Quadrilateral topologies
      case shards::Quadrilateral<8>::key:
      case shards::ShellQuadrilateral<8>::key:
        return quadrilateral_8[nodeOrd];
        break;
      case shards::Quadrilateral<9>::key:
      case shards::ShellQuadrilateral<9>::key:
        return quadrilateral_9[nodeOrd];
        break;
        
        
      // Base Tetrahedron topology
      case shards::Tetrahedron<4>::key:
        return tetrahedron[nodeOrd];
        break;
        
      // Extended Tetrahedron topologies
      case shards::Tetrahedron<8>::key:
        return tetrahedron_8[nodeOrd];
        break;
      case shards::Tetrahedron<10>::key:
        return tetrahedron_10[nodeOrd];
        break;
      case shards::Tetrahedron<11>::key:
        return tetrahedron_11[nodeOrd];
        break;

        
      // Base Hexahedron topology
      case shards::Hexahedron<8>::key:
        return hexahedron[nodeOrd];
        break;
        
      // Extended Hexahedron topologies
      case shards::Hexahedron<20>::key:
        return hexahedron_20[nodeOrd];
        break;
      case shards::Hexahedron<27>::key:
        return hexahedron_27[nodeOrd];
        break;

        
      // Base Pyramid topology  
      case shards::Pyramid<5>::key:
        return pyramid[nodeOrd];
        break;
        
      // Extended pyramid topologies
      case shards::Pyramid<13>::key:
        return pyramid_13[nodeOrd];
        break;
     case shards::Pyramid<14>::key:
        return pyramid_14[nodeOrd];
        break;
      
        
      // Base Wedge topology
      case shards::Wedge<6>::key:
        return wedge[nodeOrd];
        break;
        
      // Extended Wedge topologies
      case shards::Wedge<15>::key:
        return wedge_15[nodeOrd];
        break;
      case shards::Wedge<18>::key:
        return wedge_18[nodeOrd];
        break;
        
      default:
        TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                            ">>> ERROR (Intrepid::CellTools::getReferenceNode): invalid cell topology.");
    }
    // To disable compiler warning, should never be reached
    return line[0];
  }
  
  
  
  template<class Scalar>
  template<class ArraySubcellNode>
  void CellTools<Scalar>::getReferenceSubcellNodes(ArraySubcellNode &          subcellNodes,
                                                   const int                   subcellDim,
                                                   const int                   subcellOrd,
                                                   const shards::CellTopology& parentCell){
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( !(hasReferenceCell(parentCell) ), std::invalid_argument, 
                        ">>> ERROR (Intrepid::CellTools::getReferenceSubcellNodes): the specified cell topology does not have a reference cell.");
    
    // subcellDim can equal the cell dimension because the cell itself is a valid subcell! In this case
    // the method will return all cell cellWorkset.
    TEST_FOR_EXCEPTION( !( (0 <= subcellDim) && (subcellDim <= (int)parentCell.getDimension()) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::getReferenceSubcellNodes): subcell dimension out of range.");
    
    TEST_FOR_EXCEPTION( !( (0 <= subcellOrd) && (subcellOrd < (int)parentCell.getSubcellCount(subcellDim) ) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::getReferenceSubcellNodes): subcell ordinal out of range.");
    
    // Verify subcellNodes rank and dimensions
    {
      std::string errmsg = ">>> ERROR (Intrepid::CellTools::getReferenceSubcellNodes):";
      TEST_FOR_EXCEPTION( !( requireRankRange(errmsg, subcellNodes, 2, 2) ), std::invalid_argument, errmsg);
      
      int subcNodeCount = parentCell.getNodeCount(subcellDim, subcellOrd);
      int spaceDim = parentCell.getDimension();
      
      TEST_FOR_EXCEPTION( !( requireDimensionRange(errmsg, subcellNodes, 0,  subcNodeCount, subcNodeCount) ),
                          std::invalid_argument, errmsg);
      
      TEST_FOR_EXCEPTION( !( requireDimensionRange(errmsg, subcellNodes, 1,  spaceDim, spaceDim) ),
                          std::invalid_argument, errmsg);
    }
#endif 
    
    // Find how many cellWorkset does the specified subcell have.
    int subcNodeCount = parentCell.getNodeCount(subcellDim, subcellOrd);
    
    // Loop over subcell cellWorkset
    for(int subcNodeOrd = 0; subcNodeOrd < subcNodeCount; subcNodeOrd++){
      
      // Get the node number relative to the parent reference cell
      int cellNodeOrd = parentCell.getNodeMap(subcellDim, subcellOrd, subcNodeOrd);
            
      // Loop over node's Cartesian coordinates
      for(int dim = 0; dim < (int)parentCell.getDimension(); dim++){
        subcellNodes(subcNodeOrd, dim) = CellTools::getReferenceNode(parentCell, cellNodeOrd)[dim];
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
      case shards::Tetrahedron<11>::key:
        
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
  template<class ArrayJac, class ArrayPoint, class ArrayCell>
  void CellTools<Scalar>::setJacobian(ArrayJac &                   jacobian,
                                      const ArrayPoint &           points,
                                      const ArrayCell  &           cellWorkset,
                                      const shards::CellTopology & cellTopo,
                                      const int &                  whichCell) 
  {
    INTREPID_VALIDATE( validateArguments_setJacobian(jacobian, points, cellWorkset, whichCell,  cellTopo) );
    
    int spaceDim  = (int)cellTopo.getDimension();
    int numCells  = cellWorkset.dimension(0);
    //points can be rank-2 (P,D), or rank-3 (C,P,D)
    int numPoints = (points.rank() == 2) ? points.dimension(0) : points.dimension(1);
    
    // Jacobian is computed using gradients of an appropriate H(grad) basis function: define RCP to the base class
    Teuchos::RCP< Basis< Scalar, FieldContainer<Scalar> > > HGRAD_Basis;
    
    // Choose the H(grad) basis depending on the cell topology. \todo define maps for shells and beams
    switch( cellTopo.getKey() ){
      
      // Standard Base topologies (number of cellWorkset = number of vertices)
      case shards::Line<2>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_LINE_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
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
      case shards::Triangle<6>::key:    
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TRI_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
      case shards::Quadrilateral<9>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_QUAD_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Tetrahedron<10>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;

      case shards::Tetrahedron<11>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_COMP12_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Hexahedron<27>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_HEX_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Wedge<18>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_WEDGE_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
        // These extended topologies are not used for mapping purposes
      case shards::Quadrilateral<8>::key:
      case shards::Hexahedron<20>::key:
      case shards::Wedge<15>::key:
        TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                            ">>> ERROR (Intrepid::CellTools::setJacobian): Cell topology not supported. ");
        break;
        
        // Base and Extended Line, Beam and Shell topologies  
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
    
    // Temp (F,P,D) array for the values of basis functions gradients at the reference points
    int basisCardinality = HGRAD_Basis -> getCardinality();
    FieldContainer<Scalar> basisGrads(basisCardinality, numPoints, spaceDim);
    
    // Initialize jacobian
    for(int i = 0; i < jacobian.size(); i++){
      jacobian[i] = 0.0;
    }
        
    // Handle separately rank-2 (P,D) and rank-3 (C,P,D) cases of points arrays.
    switch(points.rank()) {
      
      // refPoints is (P,D): a single or multiple cell jacobians computed for a single set of ref. points
      case 2:
        {
          // getValues requires rank-2 (P,D) input array, but points cannot be passed directly as argument because they are a user type
          FieldContainer<Scalar> tempPoints( points.dimension(0), points.dimension(1) );
          // Copy point set corresponding to this cell oridinal to the temp (P,D) array
          for(int pt = 0; pt < points.dimension(0); pt++){
            for(int dm = 0; dm < points.dimension(1) ; dm++){
              tempPoints(pt, dm) = points(pt, dm);
            }//dm
          }//pt
          HGRAD_Basis -> getValues(basisGrads, tempPoints, OPERATOR_GRAD);
          
          // The outer loops select the multi-index of the Jacobian entry: cell, point, row, col
          // If whichCell = -1, all jacobians are computed, otherwise a single cell jacobian is computed
          int cellLoop = (whichCell == -1) ? numCells : 1 ;
          
          if(whichCell == -1) {
            for(int cellOrd = 0; cellOrd < cellLoop; cellOrd++) {
              for(int pointOrd = 0; pointOrd < numPoints; pointOrd++) {
                for(int row = 0; row < spaceDim; row++){
                  for(int col = 0; col < spaceDim; col++){
                    
                    // The entry is computed by contracting the basis index. Number of basis functions and vertices must be the same.
                    for(int bfOrd = 0; bfOrd < basisCardinality; bfOrd++){
                      jacobian(cellOrd, pointOrd, row, col) += cellWorkset(cellOrd, bfOrd, row)*basisGrads(bfOrd, pointOrd, col);
                    } // bfOrd
                  } // col
                } // row
              } // pointOrd
            } // cellOrd
          }
          else {
            for(int cellOrd = 0; cellOrd < cellLoop; cellOrd++) {
              for(int pointOrd = 0; pointOrd < numPoints; pointOrd++) {
                for(int row = 0; row < spaceDim; row++){
                  for(int col = 0; col < spaceDim; col++){
                  
                    // The entry is computed by contracting the basis index. Number of basis functions and vertices must be the same.
                    for(int bfOrd = 0; bfOrd < basisCardinality; bfOrd++){
                      jacobian(pointOrd, row, col) += cellWorkset(whichCell, bfOrd, row)*basisGrads(bfOrd, pointOrd, col);
                    } // bfOrd
                  } // col
                } // row
              } // pointOrd
            } // cellOrd
          } // if whichcell
        }// case 2
        break;
        
        // points is (C,P,D): multiple jacobians computed at multiple point sets, one jacobian per cell  
      case 3:
        {
          // getValues requires rank-2 (P,D) input array, refPoints cannot be used as argument: need temp (P,D) array
          FieldContainer<Scalar> tempPoints( points.dimension(1), points.dimension(2) );
          
          for(int cellOrd = 0; cellOrd < numCells; cellOrd++) {
            
            // Copy point set corresponding to this cell oridinal to the temp (P,D) array
            for(int pt = 0; pt < points.dimension(1); pt++){
              for(int dm = 0; dm < points.dimension(2) ; dm++){
                tempPoints(pt, dm) = points(cellOrd, pt, dm);
              }//dm
            }//pt
            
            // Compute gradients of basis functions at this set of ref. points
            HGRAD_Basis -> getValues(basisGrads, tempPoints, OPERATOR_GRAD);
            
            // Compute jacobians for the point set corresponding to the current cellordinal
            for(int pointOrd = 0; pointOrd < numPoints; pointOrd++) {
              for(int row = 0; row < spaceDim; row++){
                for(int col = 0; col < spaceDim; col++){
                  
                  // The entry is computed by contracting the basis index. Number of basis functions and vertices must be the same
                  for(int bfOrd = 0; bfOrd < basisCardinality; bfOrd++){
                    jacobian(cellOrd, pointOrd, row, col) += cellWorkset(cellOrd, bfOrd, row)*basisGrads(bfOrd, pointOrd, col);
                  } // bfOrd
                } // col
              } // row
            } // pointOrd
          }//cellOrd
        }// case 3
        
        break;
        
      default:
        TEST_FOR_EXCEPTION( !( (points.rank() == 2) && (points.rank() == 3) ), std::invalid_argument,
                            ">>> ERROR (Intrepid::CellTools::setJacobian): rank 2 or 3 required for points array. ");        
    }//switch
  }
  


template<class Scalar>
template<class ArrayJacInv, class ArrayJac>
void CellTools<Scalar>::setJacobianInv(ArrayJacInv &     jacobianInv,
                                       const ArrayJac &  jacobian) 
{
  INTREPID_VALIDATE( validateArguments_setJacobianInv(jacobianInv, jacobian) );

  RealSpaceTools<Scalar>::inverse(jacobianInv, jacobian);
}



template<class Scalar>
template<class ArrayJacDet, class ArrayJac>
void CellTools<Scalar>::setJacobianDet(ArrayJacDet &     jacobianDet,
                                       const ArrayJac &  jacobian)
{
  INTREPID_VALIDATE( validateArguments_setJacobianDetArgs(jacobianDet, jacobian) );

  RealSpaceTools<Scalar>::det(jacobianDet, jacobian);
}

//============================================================================================//
//                                                                                            //
//                      Reference-to-physical frame mapping and its inverse                   //
//                                                                                            //
//============================================================================================//

template<class Scalar>
template<class ArrayPhysPoint, class ArrayRefPoint, class ArrayCell>
void CellTools<Scalar>::mapToPhysicalFrame(ArrayPhysPoint      &        physPoints,
                                           const ArrayRefPoint &        refPoints,
                                           const ArrayCell     &        cellWorkset,
                                           const shards::CellTopology & cellTopo,
                                           const int &                  whichCell)
{
  INTREPID_VALIDATE(validateArguments_mapToPhysicalFrame( physPoints, refPoints, cellWorkset, cellTopo, whichCell) );
  
  int spaceDim  = (int)cellTopo.getDimension();
  int numCells  = cellWorkset.dimension(0);
  //points can be rank-2 (P,D), or rank-3 (C,P,D)
  int numPoints = (refPoints.rank() == 2) ? refPoints.dimension(0) : refPoints.dimension(1);
    
  // Mapping is computed using an appropriate H(grad) basis function: define RCP to the base class
  Teuchos::RCP<Basis<Scalar, FieldContainer<Scalar> > > HGRAD_Basis;
  
  // Choose the H(grad) basis depending on the cell topology. \todo define maps for shells and beams
  switch( cellTopo.getKey() ){
    
    // Standard Base topologies (number of cellWorkset = number of vertices)
    case shards::Line<2>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_LINE_C1_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
      
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
    case shards::Triangle<6>::key:    
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TRI_C2_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
      
    case shards::Quadrilateral<9>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_QUAD_C2_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
      
    case shards::Tetrahedron<10>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_C2_FEM<Scalar, FieldContainer<Scalar> >() );
      break;

    case shards::Tetrahedron<11>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_COMP12_FEM<Scalar, FieldContainer<Scalar> >() );
      break;

    case shards::Hexahedron<27>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_HEX_C2_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
      
    case shards::Wedge<18>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_WEDGE_C2_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
      
    // These extended topologies are not used for mapping purposes
    case shards::Quadrilateral<8>::key:
    case shards::Hexahedron<20>::key:
    case shards::Wedge<15>::key:
      TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                          ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrame): Cell topology not supported. ");
      break;
      
    // Base and Extended Line, Beam and Shell topologies  
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
                          ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrame): Cell topology not supported. ");
      break;
    default:
      TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                          ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrame): Cell topology not supported.");        
  }// switch  

  // Temp (F,P) array for the values of nodal basis functions at the reference points
  int basisCardinality = HGRAD_Basis -> getCardinality();
  FieldContainer<Scalar> basisVals(basisCardinality, numPoints);
  
  // Initialize physPoints
  for(int i = 0; i < physPoints.size(); i++){
    physPoints[i] = 0.0;
  }
  
  // handle separately rank-2 (P,D) and rank-3 (C,P,D) cases of refPoints
  switch(refPoints.rank()) {
    
    // refPoints is (P,D): single set of ref. points is mapped to one or multiple physical cells
    case 2:
      {
        // getValues requires rank-2 (P,D) input array, but refPoints cannot be passed directly as argument because they are a user type
        FieldContainer<Scalar> tempPoints( refPoints.dimension(0), refPoints.dimension(1) );
        // Copy point set corresponding to this cell oridinal to the temp (P,D) array
        for(int pt = 0; pt < refPoints.dimension(0); pt++){
          for(int dm = 0; dm < refPoints.dimension(1) ; dm++){
            tempPoints(pt, dm) = refPoints(pt, dm);
          }//dm
        }//pt
        HGRAD_Basis -> getValues(basisVals, tempPoints, OPERATOR_VALUE);
        
        // If whichCell = -1, ref pt. set is mapped to all cells, otherwise, the set is mapped to one cell only
        int cellLoop = (whichCell == -1) ? numCells : 1 ;
        
        // Compute the map F(refPoints) = sum node_coordinate*basis(refPoints)
        for(int cellOrd = 0; cellOrd < cellLoop; cellOrd++) {
          for(int pointOrd = 0; pointOrd < numPoints; pointOrd++) {
            for(int dim = 0; dim < spaceDim; dim++){
              for(int bfOrd = 0; bfOrd < basisCardinality; bfOrd++){
                
                if(whichCell == -1){
                  physPoints(cellOrd, pointOrd, dim) += cellWorkset(cellOrd, bfOrd, dim)*basisVals(bfOrd, pointOrd);
                }
                else{
                  physPoints(pointOrd, dim) += cellWorkset(whichCell, bfOrd, dim)*basisVals(bfOrd, pointOrd);
                }
              } // bfOrd
            }// dim
          }// pointOrd
        }//cellOrd
      }// case 2
      break;
      
    // refPoints is (C,P,D): multiple sets of ref. points are mapped to matching number of physical cells.  
    case 3:
      {
        // getValues requires rank-2 (P,D) input array, refPoints cannot be used as argument: need temp (P,D) array
        FieldContainer<Scalar> tempPoints( refPoints.dimension(1), refPoints.dimension(2) );
        
        // Compute the map F(refPoints) = sum node_coordinate*basis(refPoints)
        for(int cellOrd = 0; cellOrd < numCells; cellOrd++) {
          
          // Copy point set corresponding to this cell oridinal to the temp (P,D) array
          for(int pt = 0; pt < refPoints.dimension(1); pt++){
            for(int dm = 0; dm < refPoints.dimension(2) ; dm++){
              tempPoints(pt, dm) = refPoints(cellOrd, pt, dm);
            }//dm
          }//pt
          
          // Compute basis values for this set of ref. points
          HGRAD_Basis -> getValues(basisVals, tempPoints, OPERATOR_VALUE);
          
          for(int pointOrd = 0; pointOrd < numPoints; pointOrd++) {
            for(int dim = 0; dim < spaceDim; dim++){
              for(int bfOrd = 0; bfOrd < basisCardinality; bfOrd++){
                
                physPoints(cellOrd, pointOrd, dim) += cellWorkset(cellOrd, bfOrd, dim)*basisVals(bfOrd, pointOrd);
                
              } // bfOrd
            }// dim
          }// pointOrd
        }//cellOrd        
      }// case 3
      break;
      
    default:
      TEST_FOR_EXCEPTION( !( (refPoints.rank() == 2) && (refPoints.rank() == 3) ), std::invalid_argument,
                             ">>> ERROR (Intrepid::CellTools::mapToPhysicalFrame): rank 2 or 3 required for refPoints array. ");
  }
}



template<class Scalar>
template<class ArrayRefPoint, class ArrayPhysPoint, class ArrayCell>
void CellTools<Scalar>::mapToReferenceFrame(ArrayRefPoint        &        refPoints,
                                            const ArrayPhysPoint &        physPoints,
                                            const ArrayCell      &        cellWorkset,
                                            const shards::CellTopology &  cellTopo,
                                            const int &                   whichCell)
{
  INTREPID_VALIDATE( validateArguments_mapToReferenceFrame(refPoints, physPoints, cellWorkset, cellTopo, whichCell) );
  
  int spaceDim  = (int)cellTopo.getDimension();
  int numPoints;
  int numCells;
  
  // Define initial guesses to be  the Cell centers of the reference cell topology
  FieldContainer<Scalar> cellCenter(spaceDim);
  switch( cellTopo.getKey() ){
    // Standard Base topologies (number of cellWorkset = number of vertices)
    case shards::Line<2>::key:
      cellCenter(0) = 0.0;    break;
      
    case shards::Triangle<3>::key:
    case shards::Triangle<6>::key:    
      cellCenter(0) = 1./3.;    cellCenter(1) = 1./3.;  break;
      
    case shards::Quadrilateral<4>::key:
    case shards::Quadrilateral<9>::key:
      cellCenter(0) = 0.0;      cellCenter(1) = 0.0;    break;
      
    case shards::Tetrahedron<4>::key:
    case shards::Tetrahedron<10>::key:
    case shards::Tetrahedron<11>::key:
      cellCenter(0) = 1./6.;    cellCenter(1) =  1./6.;    cellCenter(2) =  1./6.;  break;
      
    case shards::Hexahedron<8>::key:
    case shards::Hexahedron<27>::key:
      cellCenter(0) = 0.0;      cellCenter(1) =  0.0;       cellCenter(2) =  0.0;   break;
      
    case shards::Wedge<6>::key:
    case shards::Wedge<18>::key:
      cellCenter(0) = 1./3.;    cellCenter(1) =  1./3.;     cellCenter(2) = 0.0;    break;
      
      // These extended topologies are not used for mapping purposes
    case shards::Quadrilateral<8>::key:
    case shards::Hexahedron<20>::key:
    case shards::Wedge<15>::key:
      TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                          ">>> ERROR (Intrepid::CellTools::mapToReferenceFrame): Cell topology not supported. ");
      break;
      
      // Base and Extended Line, Beam and Shell topologies  
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
                          ">>> ERROR (Intrepid::CellTools::mapToReferenceFrame): Cell topology not supported. ");
      break;
    default:
      TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                          ">>> ERROR (Intrepid::CellTools::mapToReferenceFrame): Cell topology not supported.");        
  }// switch key 
  
  // Resize initial guess depending on the rank of the physical points array
  FieldContainer<Scalar> initGuess;
  
  // Default: map (C,P,D) array of physical pt. sets to (C,P,D) array. Requires (C,P,D) initial guess.
  if(whichCell == -1){
    numPoints = physPoints.dimension(1);
    numCells = cellWorkset.dimension(0);
    initGuess.resize(numCells, numPoints, spaceDim);
    // Set initial guess:
    for(int c = 0; c < numCells; c++){
      for(int p = 0; p < numPoints; p++){
        for(int d = 0; d < spaceDim; d++){
          initGuess(c, p, d) = cellCenter(d);
        }// d
      }// p
    }// c
  }
  // Custom: map (P,D) array of physical pts. to (P,D) array. Requires (P,D) initial guess.
  else {
    numPoints = physPoints.dimension(0);
    initGuess.resize(numPoints, spaceDim);
    // Set initial guess:
    for(int p = 0; p < numPoints; p++){
      for(int d = 0; d < spaceDim; d++){
        initGuess(p, d) = cellCenter(d);
      }// d
    }// p
  }
  
  // Call method with initial guess
  mapToReferenceFrameInitGuess(refPoints, initGuess, physPoints, cellWorkset, cellTopo, whichCell);  
}
  
  

template<class Scalar>
template<class ArrayRefPoint, class ArrayInitGuess, class ArrayPhysPoint, class ArrayCell>
void CellTools<Scalar>::mapToReferenceFrameInitGuess(ArrayRefPoint        &        refPoints,
                                                     const ArrayInitGuess &        initGuess,
                                                     const ArrayPhysPoint &        physPoints,
                                                     const ArrayCell      &        cellWorkset,
                                                     const shards::CellTopology &  cellTopo,
                                                     const int &                   whichCell)
{
  INTREPID_VALIDATE( validateArguments_mapToReferenceFrame(refPoints, initGuess, physPoints, cellWorkset, cellTopo, whichCell) );

  int spaceDim  = (int)cellTopo.getDimension();
  int numPoints;
  int numCells=0;
  
  // Temp arrays for Newton iterates and Jacobians. Resize according to rank of ref. point array
  FieldContainer<Scalar> xOld;
  FieldContainer<Scalar> xTem;  
  FieldContainer<Scalar> jacobian;
  FieldContainer<Scalar> jacobInv;
  FieldContainer<Scalar> error; 
  FieldContainer<Scalar> cellCenter(spaceDim);
  
  // Default: map (C,P,D) array of physical pt. sets to (C,P,D) array. Requires (C,P,D) temp arrays and (C,P,D,D) Jacobians.
  if(whichCell == -1){
    numPoints = physPoints.dimension(1);
    numCells = cellWorkset.dimension(0);
    xOld.resize(numCells, numPoints, spaceDim);
    xTem.resize(numCells, numPoints, spaceDim);  
    jacobian.resize(numCells,numPoints, spaceDim, spaceDim);
    jacobInv.resize(numCells,numPoints, spaceDim, spaceDim);
    error.resize(numCells,numPoints); 
    // Set initial guess to xOld
    for(int c = 0; c < numCells; c++){
      for(int p = 0; p < numPoints; p++){
        for(int d = 0; d < spaceDim; d++){
          xOld(c, p, d) = initGuess(c, p, d);
        }// d
      }// p
    }// c
  }
  // Custom: map (P,D) array of physical pts. to (P,D) array. Requires (P,D) temp arrays and (P,D,D) Jacobians.
  else {
    numPoints = physPoints.dimension(0);
    xOld.resize(numPoints, spaceDim);
    xTem.resize(numPoints, spaceDim);  
    jacobian.resize(numPoints, spaceDim, spaceDim);
    jacobInv.resize(numPoints, spaceDim, spaceDim);
    error.resize(numPoints); 
    // Set initial guess to xOld
    for(int p = 0; p < numPoints; p++){
      for(int d = 0; d < spaceDim; d++){
        xOld(p, d) = initGuess(p, d);
      }// d
    }// p
  }
  
  
  // Newton method to solve the equation F(refPoints) - physPoints = 0:
  // refPoints = xOld - DF^{-1}(xOld)*(F(xOld) - physPoints) = xOld + DF^{-1}(xOld)*(physPoints - F(xOld))
  for(int iter = 0; iter < INTREPID_MAX_NEWTON; ++iter) {
    
    // Jacobians at the old iterates and their inverses. 
    setJacobian(jacobian, xOld, cellWorkset, cellTopo, whichCell);
    setJacobianInv(jacobInv, jacobian);
        
    // The Newton step.
    mapToPhysicalFrame( xTem, xOld, cellWorkset, cellTopo, whichCell );      // xTem <- F(xOld)
    RealSpaceTools<Scalar>::subtract( xTem, physPoints, xTem );        // xTem <- physPoints - F(xOld)
    RealSpaceTools<Scalar>::matvec( refPoints, jacobInv, xTem);        // refPoints <- DF^{-1}( physPoints - F(xOld) )
    RealSpaceTools<Scalar>::add( refPoints, xOld );                    // refPoints <- DF^{-1}( physPoints - F(xOld) ) + xOld
    
    // l2 error (Euclidean distance) between old and new iterates: |xOld - xNew|
    RealSpaceTools<Scalar>::subtract( xTem, xOld, refPoints );
    RealSpaceTools<Scalar>::vectorNorm( error, xTem, NORM_TWO );
    
    // Average L2 error for a multiple sets of physical points: error is rank-2 (C,P) array 
    double totalError;
    if(whichCell == -1) {
      FieldContainer<Scalar> cellWiseError(numCells);
      // error(C,P) -> cellWiseError(P)
      RealSpaceTools<Scalar>::vectorNorm( cellWiseError, error, NORM_ONE );
      totalError = RealSpaceTools<Scalar>::vectorNorm( cellWiseError, NORM_ONE );
    }
    //Average L2 error for a single set of physical points: error is rank-1 (P) array
    else{
      totalError = RealSpaceTools<Scalar>::vectorNorm( error, NORM_ONE );   
      totalError = totalError;
    }
    
    // Stopping criterion:
    if (totalError < INTREPID_TOL) {
      break;
    } 
    else if ( iter > INTREPID_MAX_NEWTON) {
      INTREPID_VALIDATE(std::cout << " Intrepid::CellTools::mapToReferenceFrameInitGuess failed to converge to desired tolerance within " 
                      << INTREPID_MAX_NEWTON  << " iterations\n" );
      break;
    }
    
    // initialize next Newton step
    xOld = refPoints;
  } // for(iter)
}



template<class Scalar>
template<class ArraySubcellPoint, class ArrayParamPoint>
void CellTools<Scalar>::mapToReferenceSubcell(ArraySubcellPoint     &       refSubcellPoints,
                                              const ArrayParamPoint &       paramPoints,
                                              const int                     subcellDim,
                                              const int                     subcellOrd,
                                              const shards::CellTopology &  parentCell){
  
  int cellDim = parentCell.getDimension();
  int numPts  = paramPoints.dimension(0);

#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( !(hasReferenceCell(parentCell) ), std::invalid_argument, 
                      ">>> ERROR (Intrepid::CellTools::mapToReferenceSubcell): the specified cell topology does not have a reference cell.");
  
  TEST_FOR_EXCEPTION( !( (1 <= subcellDim) && (subcellDim <= 2 ) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToReferenceSubcell): method defined only for 1 and 2-dimensional subcells.");  
  
  TEST_FOR_EXCEPTION( !( (0 <= subcellOrd) && (subcellOrd < (int)parentCell.getSubcellCount(subcellDim) ) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::mapToReferenceSubcell): subcell ordinal out of range.");
  
  // refSubcellPoints is rank-2 (P,D1), D1 = cell dimension
  std::string errmsg = ">>> ERROR (Intrepid::mapToReferenceSubcell):";
  TEST_FOR_EXCEPTION( !requireRankRange(errmsg, refSubcellPoints, 2,2), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, refSubcellPoints, 1, cellDim, cellDim), std::invalid_argument, errmsg);
                    
  // paramPoints is rank-2 (P,D2) with D2 = subcell dimension
  TEST_FOR_EXCEPTION( !requireRankRange(errmsg, paramPoints, 2,2), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, paramPoints, 1, subcellDim, subcellDim), std::invalid_argument, errmsg);    
  
  // cross check: refSubcellPoints and paramPoints: dimension 0 must match
  TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, refSubcellPoints, 0,  paramPoints, 0), std::invalid_argument, errmsg);      
#endif
  
  
  // Get the subcell map, i.e., the coefficients of the parametrization function for the subcell
  const FieldContainer<double>& subcellMap = getSubcellParametrization(subcellDim, parentCell);

  // Apply the parametrization map to every point in parameter domain
  if(subcellDim == 2) {
    for(int pt = 0; pt < numPts; pt++){
      double u = paramPoints(pt,0);
      double v = paramPoints(pt,1);
      
      // map_dim(u,v) = c_0(dim) + c_1(dim)*u + c_2(dim)*v because both Quad and Tri ref faces are affine!
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
                        ">>> ERROR (Intrepid::CellTools::mapToReferenceSubcell): method defined only for 1 and 2-subcells");
  }
}



template<class Scalar>
template<class ArrayEdgeTangent>
void CellTools<Scalar>::getReferenceEdgeTangent(ArrayEdgeTangent &            refEdgeTangent,
                                                const int &                   edgeOrd,
                                                const shards::CellTopology &  parentCell){
  
  int spaceDim  = parentCell.getDimension();
  
#ifdef HAVE_INTREPID_DEBUG
  
  TEST_FOR_EXCEPTION( !( (spaceDim == 2) || (spaceDim == 3) ), std::invalid_argument, 
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceTangents): two or three-dimensional parent cell required");
  
  TEST_FOR_EXCEPTION( !( (0 <= edgeOrd) && (edgeOrd < (int)parentCell.getSubcellCount(1) ) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceTangents): edge ordinal out of bounds");  
  
  TEST_FOR_EXCEPTION( !( refEdgeTangent.size() == spaceDim ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceTangents): output array size is required to match space dimension");  
#endif
  
  // Edge parametrizations are computed in setSubcellParametrization and stored in rank-3 array 
  // (subcOrd, coordinate, coefficient)
  const FieldContainer<double>& edgeMap = getSubcellParametrization(1, parentCell);
  
  // All ref. edge maps have affine coordinate functions: f_dim(u) = C_0(dim) + C_1(dim)*u, 
  //                                     => edge Tangent: -> C_1(*)
  refEdgeTangent(0) = edgeMap(edgeOrd, 0, 1);
  refEdgeTangent(1) = edgeMap(edgeOrd, 1, 1);
  
  // Skip last coordinate for 2D parent cells
  if(spaceDim == 3) {
    refEdgeTangent(2) = edgeMap(edgeOrd, 2, 1);  
  }
}



template<class Scalar>
template<class ArrayFaceTangentU, class ArrayFaceTangentV>
void CellTools<Scalar>::getReferenceFaceTangents(ArrayFaceTangentU &           uTan,
                                                 ArrayFaceTangentV &           vTan,
                                                 const int &                   faceOrd,
                                                 const shards::CellTopology &  parentCell){
  
#ifdef HAVE_INTREPID_DEBUG
  int spaceDim  = parentCell.getDimension();
  TEST_FOR_EXCEPTION( !(spaceDim == 3), std::invalid_argument, 
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceTangents): three-dimensional parent cell required");  
  
  TEST_FOR_EXCEPTION( !( (0 <= faceOrd) && (faceOrd < (int)parentCell.getSubcellCount(2) ) ), std::invalid_argument,
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
  const FieldContainer<double>& faceMap = getSubcellParametrization(2, parentCell);
  
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
template<class ArraySideNormal>
void CellTools<Scalar>::getReferenceSideNormal(ArraySideNormal &             refSideNormal,
                                               const int &                   sideOrd,
                                               const shards::CellTopology &  parentCell){
  int spaceDim  = parentCell.getDimension();
#ifdef HAVE_INTREPID_DEBUG
  
  TEST_FOR_EXCEPTION( !( (spaceDim == 2) || (spaceDim == 3) ), std::invalid_argument, 
                      ">>> ERROR (Intrepid::CellTools::getReferenceSideNormal): two or three-dimensional parent cell required");
  
  // Check side ordinal: by definition side is subcell whose dimension = spaceDim-1
  TEST_FOR_EXCEPTION( !( (0 <= sideOrd) && (sideOrd < (int)parentCell.getSubcellCount(spaceDim - 1) ) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::getReferenceSideNormal): side ordinal out of bounds");    
#endif  
  
  if(spaceDim == 2){
    
    // 2D parent cells: side = 1D subcell (edge), call the edge tangent method and rotate tangents
    getReferenceEdgeTangent(refSideNormal, sideOrd, parentCell);
    
    // rotate t(t1, t2) to get n(t2, -t1) so that (n,t) is positively oriented: det(n1,n2/t1,t2)>0
    Scalar temp = refSideNormal(0);
    refSideNormal(0) = refSideNormal(1);
    refSideNormal(1) = -temp;
  }
  else{
    // 3D parent cell: side = 2D subcell (face), call the face normal method.
    getReferenceFaceNormal(refSideNormal, sideOrd, parentCell);
  }
}
  


template<class Scalar>
template<class ArrayFaceNormal>
void CellTools<Scalar>::getReferenceFaceNormal(ArrayFaceNormal &             refFaceNormal,
                                               const int &                   faceOrd,
                                               const shards::CellTopology &  parentCell){
  int spaceDim  = parentCell.getDimension();
  
#ifdef HAVE_INTREPID_DEBUG
  
  TEST_FOR_EXCEPTION( !(spaceDim == 3), std::invalid_argument, 
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceNormal): three-dimensional parent cell required");  
  
  TEST_FOR_EXCEPTION( !( (0 <= faceOrd) && (faceOrd < (int)parentCell.getSubcellCount(2) ) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceNormal): face ordinal out of bounds");  
  
  TEST_FOR_EXCEPTION( !( refFaceNormal.rank() == 1 ), std::invalid_argument,  
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceNormal): rank = 1 required for output array"); 
    
  TEST_FOR_EXCEPTION( !( refFaceNormal.dimension(0) == spaceDim ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::getReferenceFaceNormal): dim0 (spatial dim) must match that of parent cell");  
#endif

  // Reference face normal = vector product of reference face tangents. Allocate temp FC storage:
  FieldContainer<Scalar> uTan(spaceDim);
  FieldContainer<Scalar> vTan(spaceDim);
  getReferenceFaceTangents(uTan, vTan, faceOrd, parentCell);
  
  // Compute the vector product of the reference face tangents:
  RealSpaceTools<Scalar>::vecprod(refFaceNormal, uTan, vTan);
}



template<class Scalar>
template<class ArrayEdgeTangent, class ArrayJac>
void CellTools<Scalar>::getPhysicalEdgeTangents(ArrayEdgeTangent &            edgeTangents,
                                                const ArrayJac &              worksetJacobians,
                                                const int &                   worksetEdgeOrd,
                                                const shards::CellTopology &  parentCell){
  int worksetSize = worksetJacobians.dimension(0);
  int edgePtCount = worksetJacobians.dimension(1); 
  int pCellDim    = parentCell.getDimension();
  
#ifdef HAVE_INTREPID_DEBUG
  std::string errmsg = ">>> ERROR (Intrepid::CellTools::getPhysicalEdgeTangents):";
  
  TEST_FOR_EXCEPTION( !( (pCellDim == 3) || (pCellDim == 2) ), std::invalid_argument, 
                      ">>> ERROR (Intrepid::CellTools::getPhysicalEdgeTangents): 2D or 3D parent cell required");  
  
  // (1) edgeTangents is rank-3 (C,P,D) and D=2, or 3 is required
  TEST_FOR_EXCEPTION( !requireRankRange(errmsg, edgeTangents, 3,3), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, edgeTangents, 2, 2,3), std::invalid_argument, errmsg);
 
  // (2) worksetJacobians in rank-4 (C,P,D,D) and D=2, or 3 is required
  TEST_FOR_EXCEPTION( !requireRankRange(errmsg, worksetJacobians, 4,4), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 2, 2,3), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 3, 2,3), std::invalid_argument, errmsg);
  
  // (4) cross-check array dimensions: edgeTangents (C,P,D) vs. worksetJacobians (C,P,D,D)
  TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, edgeTangents, 0,1,2,2,  worksetJacobians, 0,1,2,3), std::invalid_argument, errmsg);      
  
#endif
  
  // Temp storage for constant reference edge tangent: rank-1 (D) arrays
  FieldContainer<double> refEdgeTan(pCellDim);
  getReferenceEdgeTangent(refEdgeTan, worksetEdgeOrd, parentCell);
  
  // Loop over workset faces and edge points
  for(int pCell = 0; pCell < worksetSize; pCell++){
    for(int pt = 0; pt < edgePtCount; pt++){
      
      // Apply parent cell Jacobian to ref. edge tangent
      for(int i = 0; i < pCellDim; i++){
        edgeTangents(pCell, pt, i) = 0.0;
        for(int j = 0; j < pCellDim; j++){
          edgeTangents(pCell, pt, i) +=  worksetJacobians(pCell, pt, i, j)*refEdgeTan(j);
        }// for j
      }// for i
    }// for pt
  }// for pCell
}



template<class Scalar>
template<class ArrayFaceTangentU, class ArrayFaceTangentV, class ArrayJac>
void CellTools<Scalar>::getPhysicalFaceTangents(ArrayFaceTangentU &           faceTanU,
                                                ArrayFaceTangentV &           faceTanV,
                                                const ArrayJac &              worksetJacobians,
                                                const int &                   worksetFaceOrd,
                                                const shards::CellTopology &  parentCell){
  int worksetSize = worksetJacobians.dimension(0);
  int facePtCount = worksetJacobians.dimension(1); 
  int pCellDim    = parentCell.getDimension();
  
#ifdef HAVE_INTREPID_DEBUG
  std::string errmsg = ">>> ERROR (Intrepid::CellTools::getPhysicalFaceTangents):";

  TEST_FOR_EXCEPTION( !(pCellDim == 3), std::invalid_argument, 
                      ">>> ERROR (Intrepid::CellTools::getPhysicalFaceTangents): three-dimensional parent cell required");  
  
  // (1) faceTanU and faceTanV are rank-3 (C,P,D) and D=3 is required
  TEST_FOR_EXCEPTION( !requireRankRange(errmsg, faceTanU, 3,3), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, faceTanU, 2, 3,3), std::invalid_argument, errmsg);

  TEST_FOR_EXCEPTION( !requireRankRange(errmsg, faceTanV, 3,3), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, faceTanV, 2, 3,3), std::invalid_argument, errmsg);

  TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, faceTanU,  faceTanV), std::invalid_argument, errmsg);      

  // (3) worksetJacobians in rank-4 (C,P,D,D) and D=3 is required
  TEST_FOR_EXCEPTION( !requireRankRange(errmsg, worksetJacobians, 4,4), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 2, 3,3), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 3, 3,3), std::invalid_argument, errmsg);

  // (4) cross-check array dimensions: faceTanU (C,P,D) vs. worksetJacobians (C,P,D,D)
  TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, faceTanU, 0,1,2,2,  worksetJacobians, 0,1,2,3), std::invalid_argument, errmsg);      

#endif
    
  // Temp storage for the pair of constant ref. face tangents: rank-1 (D) arrays
  FieldContainer<double> refFaceTanU(pCellDim);
  FieldContainer<double> refFaceTanV(pCellDim);
  getReferenceFaceTangents(refFaceTanU, refFaceTanV, worksetFaceOrd, parentCell);

  // Loop over workset faces and face points
  for(int pCell = 0; pCell < worksetSize; pCell++){
    for(int pt = 0; pt < facePtCount; pt++){
      
      // Apply parent cell Jacobian to ref. face tangents
      for(int dim = 0; dim < pCellDim; dim++){
        faceTanU(pCell, pt, dim) = 0.0;
        faceTanV(pCell, pt, dim) = 0.0;
        
        // Unroll loops: parent cell dimension can only be 3
        faceTanU(pCell, pt, dim) = \
          worksetJacobians(pCell, pt, dim, 0)*refFaceTanU(0) + \
          worksetJacobians(pCell, pt, dim, 1)*refFaceTanU(1) + \
          worksetJacobians(pCell, pt, dim, 2)*refFaceTanU(2);
        faceTanV(pCell, pt, dim) = \
          worksetJacobians(pCell, pt, dim, 0)*refFaceTanV(0) + \
          worksetJacobians(pCell, pt, dim, 1)*refFaceTanV(1) + \
          worksetJacobians(pCell, pt, dim, 2)*refFaceTanV(2);
      }// for dim
    }// for pt
  }// for pCell
}


template<class Scalar>
template<class ArraySideNormal, class ArrayJac>
void CellTools<Scalar>::getPhysicalSideNormals(ArraySideNormal &             sideNormals,
                                               const ArrayJac &              worksetJacobians,
                                               const int &                   worksetSideOrd,
                                               const shards::CellTopology &  parentCell){
  int worksetSize = worksetJacobians.dimension(0);
  int sidePtCount = worksetJacobians.dimension(1);   
  int spaceDim  = parentCell.getDimension();
  
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( !( (spaceDim == 2) || (spaceDim == 3) ), std::invalid_argument, 
                      ">>> ERROR (Intrepid::CellTools::getPhysicalSideNormals): two or three-dimensional parent cell required");
  
  // Check side ordinal: by definition side is subcell whose dimension = spaceDim-1
  TEST_FOR_EXCEPTION( !( (0 <= worksetSideOrd) && (worksetSideOrd < (int)parentCell.getSubcellCount(spaceDim - 1) ) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::getPhysicalSideNormals): side ordinal out of bounds");  
#endif  
  
  if(spaceDim == 2){

    // 2D parent cells: side = 1D subcell (edge), call the edge tangent method and rotate tangents
    getPhysicalEdgeTangents(sideNormals, worksetJacobians, worksetSideOrd, parentCell);
    
    // rotate t(t1, t2) to get n(t2, -t1) so that (n,t) is positively oriented: det(n1,n2/t1,t2)>0
    for(int cell = 0; cell < worksetSize; cell++){
      for(int pt = 0; pt < sidePtCount; pt++){
        Scalar temp = sideNormals(cell, pt, 0);
        sideNormals(cell, pt, 0) = sideNormals(cell, pt, 1);
        sideNormals(cell, pt, 1) = -temp;
      }// for pt
    }// for cell
  }
  else{
    // 3D parent cell: side = 2D subcell (face), call the face normal method.
    getPhysicalFaceNormals(sideNormals, worksetJacobians, worksetSideOrd, parentCell);
  }
}
  
  
  
template<class Scalar>
template<class ArrayFaceNormal, class ArrayJac>
void CellTools<Scalar>::getPhysicalFaceNormals(ArrayFaceNormal &             faceNormals,
                                               const ArrayJac &              worksetJacobians,
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
  
  // (3) worksetJacobians in rank-4 (C,P,D,D) and D=3 is required
  TEST_FOR_EXCEPTION( !requireRankRange(errmsg, worksetJacobians, 4,4), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 2, 3,3), std::invalid_argument, errmsg);
  TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 3, 3,3), std::invalid_argument, errmsg);
  
  // (4) cross-check array dimensions: faceNormals (C,P,D) vs. worksetJacobians (C,P,D,D)
  TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, faceNormals, 0,1,2,2,  worksetJacobians, 0,1,2,3), std::invalid_argument, errmsg);        
#endif
  
  // Temp storage for physical face tangents: rank-3 (C,P,D) arrays
  FieldContainer<double> faceTanU(worksetSize, facePtCount, pCellDim);
  FieldContainer<double> faceTanV(worksetSize, facePtCount, pCellDim);
  getPhysicalFaceTangents(faceTanU, faceTanV, worksetJacobians, worksetFaceOrd, parentCell);
  
  // Compute the vector product of the physical face tangents:
  RealSpaceTools<Scalar>::vecprod(faceNormals, faceTanU, faceTanV);
  
  
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
  unsigned key = cellTopo.getBaseCellTopologyData() -> key ;
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
  
  int rank = points.rank();  
  
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( !( (1 <= points.rank() ) && (points.rank() <= 3) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::checkPointsetInclusion): rank-1, 2 or 3 required for input points array. ");

  // The last dimension of points array at (rank - 1) is the spatial dimension. Must equal the cell dimension.
  TEST_FOR_EXCEPTION( !( points.dimension(rank - 1) == (int)cellTopo.getDimension() ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::checkPointsetInclusion): Point and cell dimensions do not match. ");
#endif
  
  // create temp output array depending on the rank of the input array 
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
template<class ArrayIncl, class ArrayPoint>
void CellTools<Scalar>::checkPointwiseInclusion(ArrayIncl &                   inRefCell,
                                                const ArrayPoint &            points,
                                                const shards::CellTopology &  cellTopo, 
                                                const double &                threshold) {
  int apRank   = points.rank();
  
#ifdef HAVE_INTREPID_DEBUG
  
  // Verify that points and inRefCell have correct ranks and dimensions
  std::string errmsg = ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusion):";
  if(points.rank() == 1) {
    TEST_FOR_EXCEPTION( !(inRefCell.rank() == 1 ), std::invalid_argument, 
                        ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusion): rank-1 input array requires rank-1 output array.");  
    TEST_FOR_EXCEPTION( !(inRefCell.dimension(0) == 1), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusion): rank-1 input array requires dim0 = 1 for output array.");  
  }
  else if(points.rank() == 2){
    TEST_FOR_EXCEPTION( !(inRefCell.rank() == 1 ), std::invalid_argument, 
                        ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusion): rank-2 input array requires rank-1 output array.");  
    // dimension 0 of the arrays must match
    TEST_FOR_EXCEPTION( !requireDimensionMatch( errmsg, inRefCell, 0,  points, 0), std::invalid_argument, errmsg);
  }
  else if (points.rank() == 3) {
    TEST_FOR_EXCEPTION( !(inRefCell.rank() == 2 ), std::invalid_argument, 
                        ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusion): rank-3 input array requires rank-2 output array.");  
    // dimensions 0 and 1 of the arrays must match
    TEST_FOR_EXCEPTION( !requireDimensionMatch( errmsg, inRefCell, 0,1,  points, 0,1), std::invalid_argument, errmsg);
  }
  else{
    TEST_FOR_EXCEPTION( !( (points.rank() == 1) || (points.rank() == 2) || (points.rank() == 3) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusion): rank-1, 2 or 3 required for input points array. ");      
  }    
  
  // The last dimension of points array at (rank - 1) is the spatial dimension. Must equal the cell dimension.
  TEST_FOR_EXCEPTION( !( points.dimension(apRank - 1) == (int)cellTopo.getDimension() ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusion): Point and cell dimensions do not match. ");
  
#endif
  
  // Initializations
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
      TEST_FOR_EXCEPTION( !( (1 <= points.rank() ) && (points.rank() <= 3) ), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::checkPointwiseInclusion): rank-1, 2 or 3 required for input points array. ");      
  }// switch
  
  
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
template<class ArrayIncl, class ArrayPoint, class ArrayCell>
void CellTools<Scalar>::checkPointwiseInclusion(ArrayIncl &                   inCell,
                                                const ArrayPoint &            points,
                                                const ArrayCell &             cellWorkset,
                                                const shards::CellTopology &  cell,
                                                const int &                   whichCell, 
                                                const double &                threshold)
{
  INTREPID_VALIDATE( validateArguments_checkPointwiseInclusion(inCell, points, cellWorkset, whichCell, cell) );
  
  // For cell topologies with reference cells this test maps the points back to the reference cell
  // and uses the method for reference cells
  unsigned baseKey = cell.getBaseCellTopologyData() -> key;
  
  switch(baseKey){
    
    case shards::Line<>::key :
    case shards::Triangle<>::key:
    case shards::Quadrilateral<>::key :
    case shards::Tetrahedron<>::key :
    case shards::Hexahedron<>::key :
    case shards::Wedge<>::key :
    case shards::Pyramid<>::key :
      {
        FieldContainer<Scalar> refPoints;
        
        if(points.rank() == 2){
          refPoints.resize(points.dimension(0), points.dimension(1) );
          mapToReferenceFrame(refPoints, points, cellWorkset, cell, whichCell);
          checkPointwiseInclusion(inCell, refPoints, cell, threshold );
        }
        else if(points.rank() == 3){
          refPoints.resize(points.dimension(0), points.dimension(1), points.dimension(2) );
          mapToReferenceFrame(refPoints, points, cellWorkset, cell, whichCell);
          checkPointwiseInclusion(inCell, refPoints, cell, threshold );          
        }
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
template<class ArrayJac, class ArrayPoint, class ArrayCell>
void CellTools<Scalar>::validateArguments_setJacobian(const ArrayJac    &          jacobian,
                                                      const ArrayPoint  &          points,
                                                      const ArrayCell   &          cellWorkset,
                                                      const int &                  whichCell,
                                                      const shards::CellTopology & cellTopo){
  
  // Validate cellWorkset array
  TEST_FOR_EXCEPTION( (cellWorkset.rank() != 3), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): rank = 3 required for cellWorkset array");
  
  TEST_FOR_EXCEPTION( (cellWorkset.dimension(0) <= 0), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 0 (number of cells) >= 1 required for cellWorkset array");
  
  TEST_FOR_EXCEPTION( (cellWorkset.dimension(1) != (int)cellTopo.getSubcellCount(0) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 1 (number of cell nodes) of cellWorkset array does not match cell topology");
  
  TEST_FOR_EXCEPTION( (cellWorkset.dimension(2) != (int)cellTopo.getDimension() ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 2 (spatial dimension) of cellWorkset array  does not match cell dimension");
    
  // validate whichCell. It can be either -1 (default value) or a valid cell ordinal.
  TEST_FOR_EXCEPTION( !( ( (0 <= whichCell ) && (whichCell < cellWorkset.dimension(0) ) ) || (whichCell == -1) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): whichCell = -1 or a valid cell ordinal is required.");
  
  
  // Validate points array: can be rank-2 (P,D) or rank-3 (C,P,D)
  // If rank-2: admissible jacobians: rank-3 (P,D,D) or rank-4 (C,P,D,D); admissible whichCell: -1 (default) or cell ordinal.
  if(points.rank() == 2) {
    TEST_FOR_EXCEPTION( (points.dimension(0) <= 0), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 0 (number of points) >= 1 required for points array ");
    
    TEST_FOR_EXCEPTION( (points.dimension(1) != (int)cellTopo.getDimension() ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 1 (spatial dimension) of points array does not match cell dimension");
    
    // Validate the output array for the Jacobian: if whichCell == -1 all Jacobians are computed, rank-4 (C,P,D,D) required
    if(whichCell == -1) {
      TEST_FOR_EXCEPTION( (jacobian.rank() != 4), std::invalid_argument, 
                          ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): rank = 4 required for jacobian array");
      
      TEST_FOR_EXCEPTION( (jacobian.dimension(0) != cellWorkset.dimension(0)), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 0 (number of cells) of jacobian array must equal dim 0 of cellWorkset array");
      
      TEST_FOR_EXCEPTION( (jacobian.dimension(1) != points.dimension(0)), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 1 (number of points) of jacobian array must equal dim 0 of points array");

      TEST_FOR_EXCEPTION( (jacobian.dimension(2) != points.dimension(1)), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 2 (spatial dimension) of jacobian array must equal dim 1 of points array");
      
      TEST_FOR_EXCEPTION( !(jacobian.dimension(2) == jacobian.dimension(3) ), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 2 = dim 3 (same spatial dimensions) required for jacobian array. ");
      
      TEST_FOR_EXCEPTION( !( (0 < jacobian.dimension(3) ) && (jacobian.dimension(3) < 4) ), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 2 and dim 3 (spatial dimensions) must be between 1 and 3. ");
    }     
    // A single cell Jacobian is computed when whichCell != -1 (whichCell has been already validated), rank-3 (P,D,D) required
    else {
      TEST_FOR_EXCEPTION( (jacobian.rank() != 3), std::invalid_argument, 
                          ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): rank = 3 required for jacobian array");
      
      TEST_FOR_EXCEPTION( (jacobian.dimension(0) != points.dimension(0)), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 0 (number of points) of jacobian array must equal dim 0 of points array");

      TEST_FOR_EXCEPTION( (jacobian.dimension(1) != points.dimension(1)), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 1 (spatial dimension) of jacobian array must equal dim 1 of points array");
      
      TEST_FOR_EXCEPTION( !(jacobian.dimension(1) == jacobian.dimension(2) ), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 1 = dim 2 (same spatial dimensions) required for jacobian array. ");
      
      TEST_FOR_EXCEPTION( !( (0 < jacobian.dimension(1) ) && (jacobian.dimension(1) < 4) ), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 1 and dim 2 (spatial dimensions) must be between 1 and 3. ");
    }
  }
  // Point array is rank-3 (C,P,D): requires whichCell = -1 and rank-4 (C,P,D,D) jacobians
  else if(points.rank() ==3){
    std::string errmsg  = ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian):";
    TEST_FOR_EXCEPTION( (points.dimension(0) != cellWorkset.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 0 (number of cells) of points array must equal dim 0 of cellWorkset array");

    TEST_FOR_EXCEPTION( (points.dimension(1) <= 0), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 1 (number of points) >= 1 required for points array ");
    
    TEST_FOR_EXCEPTION( (points.dimension(2) != (int)cellTopo.getDimension() ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 2 (spatial dimension) of points array does not match cell dimension");
    
    TEST_FOR_EXCEPTION( (whichCell != -1), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): default value whichCell=-1 required for rank-3 input points");
    
    // rank-4 (C,P,D,D) jacobian required for rank-3 (C,P,D) input points
    TEST_FOR_EXCEPTION( !requireRankRange(errmsg, jacobian,  4, 4), std::invalid_argument,errmsg);
    
    TEST_FOR_EXCEPTION( (jacobian.dimension(0) != points.dimension(0)), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 0 (number of cells) of jacobian array must equal dim 0 of points array");
    
    TEST_FOR_EXCEPTION( (jacobian.dimension(1) != points.dimension(1)), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 1 (number of points) of jacobian array must equal dim 1 of points array");
  
    TEST_FOR_EXCEPTION( (jacobian.dimension(2) != points.dimension(2)), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 2 (spatial dimension) of jacobian array must equal dim 2 of points array");
    
    TEST_FOR_EXCEPTION( !(jacobian.dimension(2) == jacobian.dimension(3) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 2 = dim 3 (same spatial dimensions) required for jacobian array. ");
    
    TEST_FOR_EXCEPTION( !( (0 < jacobian.dimension(3) ) && (jacobian.dimension(3) < 4) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): dim 2 and dim 3 (spatial dimensions) must be between 1 and 3. ");
  }
  else {
    TEST_FOR_EXCEPTION( !( (points.rank() == 2) && (points.rank() ==3) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobian): rank = 2 or 3 required for points array");
  }  
}



template<class Scalar>
template<class ArrayJacInv, class ArrayJac>
void CellTools<Scalar>::validateArguments_setJacobianInv(const ArrayJacInv & jacobianInv,
                                                         const ArrayJac &    jacobian)
{
  // Validate input jacobian array: admissible ranks & dimensions are: 
  // - rank-4 with dimensions (C,P,D,D), or rank-3 with dimensions (P,D,D).
  int jacobRank = jacobian.rank();
  TEST_FOR_EXCEPTION( !( (jacobRank == 4) || (jacobRank == 3) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobianInv): rank = 4 or 3 required for jacobian array. ");
  
  // Verify correctness of spatial dimensions - they are the last two dimensions of the array: rank-2 and rank-1
  TEST_FOR_EXCEPTION( !(jacobian.dimension(jacobRank - 1) == jacobian.dimension(jacobRank - 2) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobianInv): dim(rank-2) = dim(rank-2) (same spatial dimensions) required for jacobian array. ");
  
  TEST_FOR_EXCEPTION( !( (0 < jacobian.dimension(jacobRank - 1) ) && (jacobian.dimension(jacobRank - 1) < 4) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobianInv): dim(rank-1) and dim(rank-2) (spatial dimensions) must be between 1 and 3. ");
  
  // Validate output jacobianInv array: must have the same rank and dimensions as the input array.
  std::string errmsg = ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobianInv):";

  TEST_FOR_EXCEPTION( !(requireRankMatch(errmsg, jacobianInv, jacobian) ), std::invalid_argument, errmsg);
  
  TEST_FOR_EXCEPTION( !(requireDimensionMatch(errmsg, jacobianInv, jacobian) ), std::invalid_argument, errmsg);
}



template<class Scalar>
template<class ArrayJacDet, class ArrayJac>
void CellTools<Scalar>::validateArguments_setJacobianDetArgs(const ArrayJacDet &  jacobianDet,
                                                             const ArrayJac    &  jacobian)
{
  // Validate input jacobian array: admissible ranks & dimensions are: 
  // - rank-4 with dimensions (C,P,D,D), or rank-3 with dimensions (P,D,D).
  int jacobRank = jacobian.rank();
  TEST_FOR_EXCEPTION( !( (jacobRank == 4) || (jacobRank == 3) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobianInv): rank = 4 or 3 required for jacobian array. ");
  
  // Verify correctness of spatial dimensions - they are the last two dimensions of the array: rank-2 and rank-1
  TEST_FOR_EXCEPTION( !(jacobian.dimension(jacobRank - 1) == jacobian.dimension(jacobRank - 2) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobianInv): dim(rank-2) = dim(rank-2) (same spatial dimensions) required for jacobian array. ");
  
  TEST_FOR_EXCEPTION( !( (0 < jacobian.dimension(jacobRank - 1) ) && (jacobian.dimension(jacobRank - 1) < 4) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobianInv): dim(rank-1) and dim(rank-2) (spatial dimensions) must be between 1 and 3. ");

  
  // Validate output jacobianDet array: must be rank-2 with dimensions (C,P) if jacobian was rank-4:
  if(jacobRank == 4){
    TEST_FOR_EXCEPTION( !(jacobianDet.rank() == 2), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobianDetArgs): rank = 2 required for jacobianDet if jacobian is rank-4. ");
    
    TEST_FOR_EXCEPTION( !(jacobianDet.dimension(0) == jacobian.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobianDetArgs): dim 0 (number of cells) of jacobianDet array must equal dim 0 of jacobian array. ");
    
    TEST_FOR_EXCEPTION( !(jacobianDet.dimension(1) == jacobian.dimension(1) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobianDetArgs): dim 1 (number of points) of jacobianDet array must equal dim 1 of jacobian array.");  
  }
  
  // must be rank-1 with dimension (P) if jacobian was rank-3
  else {
    TEST_FOR_EXCEPTION( !(jacobianDet.rank() == 1), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobianDetArgs): rank = 1 required for jacobianDet if jacobian is rank-3. ");
    
    TEST_FOR_EXCEPTION( !(jacobianDet.dimension(0) == jacobian.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_setJacobianDetArgs): dim 0 (number of points) of jacobianDet array must equal dim 0 of jacobian array.");  
  }
}



template<class Scalar>
template<class ArrayPhysPoint, class ArrayRefPoint, class ArrayCell>
void CellTools<Scalar>::validateArguments_mapToPhysicalFrame(const ArrayPhysPoint &        physPoints,
                                                             const ArrayRefPoint  &        refPoints,
                                                             const ArrayCell      &        cellWorkset,
                                                             const shards::CellTopology &  cellTopo,
                                                             const int&                    whichCell)
{
  std::string errmsg = ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame):";
  
  // Validate cellWorkset array
  TEST_FOR_EXCEPTION( (cellWorkset.rank() != 3), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): rank = 3 required for cellWorkset array");
  
  TEST_FOR_EXCEPTION( (cellWorkset.dimension(0) <= 0), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): dim 0 (number of cells) >= 1 required for cellWorkset array");
  
  TEST_FOR_EXCEPTION( (cellWorkset.dimension(1) != (int)cellTopo.getSubcellCount(0) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): dim 1 (number of cell nodes) of cellWorkset array does not match cell topology");
  
  TEST_FOR_EXCEPTION( (cellWorkset.dimension(2) != (int)cellTopo.getDimension() ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): dim 2 (spatial dimension) of cellWorkset array  does not match cell dimension");
  
    
  // validate whichCell. It can be either -1 (default value) or a valid cell ordinal.
  TEST_FOR_EXCEPTION( !( ( (0 <= whichCell ) && (whichCell < cellWorkset.dimension(0) ) ) || (whichCell == -1) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): whichCell = -1 or a valid cell ordinal is required.");
  
  // Validate refPoints array: can be rank-2 (P,D) or rank-3 (C,P,D) array
  // If rank-2: admissible output array is (P,D) or (C,P,D); admissible whichCell: -1 (default) or cell ordinal
  if(refPoints.rank() == 2) {
    TEST_FOR_EXCEPTION( (refPoints.dimension(0) <= 0), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): dim 0 (number of points) >= 1 required for refPoints array ");
    
    TEST_FOR_EXCEPTION( (refPoints.dimension(1) != (int)cellTopo.getDimension() ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): dim 1 (spatial dimension) of refPoints array does not match cell dimension");

    // Validate output array: whichCell = -1 requires rank-3 array with dimensions (C,P,D)  
    if(whichCell == -1) {
      TEST_FOR_EXCEPTION( ( (physPoints.rank() != 3) && (whichCell == -1) ), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): rank = 3 required for physPoints array for the default whichCell value");
      
      TEST_FOR_EXCEPTION( (physPoints.dimension(0) != cellWorkset.dimension(0)), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): dim 0 (number of cells) of physPoints array must equal dim 0 of cellWorkset array");
      
      TEST_FOR_EXCEPTION( (physPoints.dimension(1) != refPoints.dimension(0)), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): dim 1 (number of points) of physPoints array must equal dim 0 of refPoints array"); 
      
      TEST_FOR_EXCEPTION( (physPoints.dimension(2) != (int)cellTopo.getDimension()), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): dim 2 (spatial dimension) does not match cell dimension ");  
    }
    // 0 <= whichCell < num cells requires rank-2 (P,D) arrays for both refPoints and physPoints
    else{
      TEST_FOR_EXCEPTION( (physPoints.rank() != 2), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): rank = 2 required for physPoints array");
      
      TEST_FOR_EXCEPTION( (physPoints.dimension(0) != refPoints.dimension(0)), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): dim 0 (number of points) of physPoints array must equal dim 0 of refPoints array"); 
      
      TEST_FOR_EXCEPTION( (physPoints.dimension(1) != (int)cellTopo.getDimension()), std::invalid_argument,
                          ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): dim 1 (spatial dimension) does not match cell dimension ");      
    }
  }
  // refPoints is (C,P,D): requires physPoints to be (C,P,D) and whichCell=-1  (because all cell mappings are applied)
  else if(refPoints.rank() == 3) {
    
    // 1. validate refPoints dimensions and rank
    TEST_FOR_EXCEPTION( (refPoints.dimension(0) != cellWorkset.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): dim 0 (number of cells) of refPoints and cellWorkset arraya are required to match ");

    TEST_FOR_EXCEPTION( (refPoints.dimension(1) <= 0), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): dim 1 (number of points) >= 1 required for refPoints array ");
    
    TEST_FOR_EXCEPTION( (refPoints.dimension(2) != (int)cellTopo.getDimension() ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): dim 2 (spatial dimension) of refPoints array does not match cell dimension");
    
    // 2. whichCell  must be -1
    TEST_FOR_EXCEPTION( (whichCell != -1), std::invalid_argument, 
                        ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): default value is required for rank-3 refPoints array");

    // 3.  physPoints must match rank and dimensions of refPoints
    TEST_FOR_EXCEPTION( !requireRankMatch(errmsg, refPoints, physPoints), std::invalid_argument, errmsg );
    TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, refPoints, physPoints), std::invalid_argument, errmsg);
  }
  // if rank is not 2 or 3 throw exception
  else {
    TEST_FOR_EXCEPTION( !( (refPoints.rank() == 2) || (refPoints.rank() == 3) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_mapToPhysicalFrame): rank = 2 or 3 required for refPoints array");
  }
}



template<class Scalar>
template<class ArrayRefPoint, class ArrayPhysPoint, class ArrayCell>
void CellTools<Scalar>::validateArguments_mapToReferenceFrame(const ArrayRefPoint  &        refPoints,
                                                              const ArrayPhysPoint &        physPoints,
                                                              const ArrayCell      &        cellWorkset,
                                                              const shards::CellTopology &  cellTopo,
                                                              const int&                    whichCell)
{
  std::string errmsg  = ">>> ERROR (Intrepid::CellTools::validateArguments_mapToReferenceFrame):";
  std::string errmsg1 = ">>> ERROR (Intrepid::CellTools::validateArguments_mapToReferenceFrame):";
  
  // Validate cellWorkset array
  TEST_FOR_EXCEPTION( (cellWorkset.rank() != 3), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_mapToReferenceFrame): rank = 3 required for cellWorkset array");
  
  TEST_FOR_EXCEPTION( (cellWorkset.dimension(0) <= 0), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_mapToReferenceFrame): dim 0 (number of cells) >= 1 required for cellWorkset array");
  
  TEST_FOR_EXCEPTION( (cellWorkset.dimension(1) != (int)cellTopo.getSubcellCount(0) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_mapToReferenceFrame): dim 1 (number of cell nodes) of cellWorkset array does not match cell topology");
  
  TEST_FOR_EXCEPTION( (cellWorkset.dimension(2) != (int)cellTopo.getDimension() ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_mapToReferenceFrame): dim 2 (spatial dimension) of cellWorkset array  does not match cell dimension");
    
  // Validate whichCell. It can be either -1 (default value) or a valid cell ordinal.
  TEST_FOR_EXCEPTION( !( ( (0 <= whichCell ) && (whichCell < cellWorkset.dimension(0) ) ) || (whichCell == -1) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_mapToReferenceFrame): whichCell = -1 or a valid cell ordinal is required.");
  
  // Admissible ranks and dimensions of refPoints and physPoints depend on whichCell value:
  // default is to map multiple sets of points to multiple sets of points. (C,P,D) arrays required
  int validRank;
  if(whichCell == -1) {
    validRank = 3;
    errmsg1 += " default value of whichCell requires rank-3 arrays:";
  }
  // whichCell is valid cell ordinal => we map single set of pts to a single set of pts. (P,D) arrays required
  else{
    errmsg1 += " rank-2 arrays required when whichCell is valid cell ordinal";
    validRank = 2;
  }
  TEST_FOR_EXCEPTION( !requireRankRange(errmsg1, refPoints,  validRank,validRank), std::invalid_argument, errmsg1);
  TEST_FOR_EXCEPTION( !requireRankMatch(errmsg1, physPoints, refPoints),           std::invalid_argument, errmsg1);
  TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg1, refPoints, physPoints),      std::invalid_argument, errmsg1);
}



template<class Scalar>
template<class ArrayRefPoint, class ArrayInitGuess, class ArrayPhysPoint, class ArrayCell>
void CellTools<Scalar>::validateArguments_mapToReferenceFrame(const ArrayRefPoint  &        refPoints,
                                                              const ArrayInitGuess &        initGuess,
                                                              const ArrayPhysPoint &        physPoints,
                                                              const ArrayCell      &        cellWorkset,
                                                              const shards::CellTopology &  cellTopo,
                                                              const int&                    whichCell)
{
  // Call the method that validates arguments with the default initial guess selection
  validateArguments_mapToReferenceFrame(refPoints, physPoints, cellWorkset, cellTopo, whichCell);
  
  // Then check initGuess: its rank and dimensions must match those of physPoints.
  std::string errmsg = ">>> ERROR (Intrepid::CellTools::validateArguments_mapToReferenceFrame):";
  TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, initGuess, physPoints), std::invalid_argument, errmsg);  
}


template<class Scalar>
template<class ArrayIncl, class ArrayPoint, class ArrayCell>
void CellTools<Scalar>::validateArguments_checkPointwiseInclusion(ArrayIncl &                   inCell,
                                                                  const ArrayPoint &            physPoints,
                                                                  const ArrayCell &             cellWorkset,
                                                                  const int &                   whichCell,
                                                                  const shards::CellTopology &  cell)
{
  // Validate cellWorkset array
  TEST_FOR_EXCEPTION( (cellWorkset.rank() != 3), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): rank = 3 required for cellWorkset array");
  
  TEST_FOR_EXCEPTION( (cellWorkset.dimension(0) <= 0), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): dim 0 (number of cells) >= 1 required for cellWorkset array");
  
  TEST_FOR_EXCEPTION( (cellWorkset.dimension(1) != (int)cell.getSubcellCount(0) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): dim 1 (number of cell nodes) of cellWorkset array does not match cell topology");
  
  TEST_FOR_EXCEPTION( (cellWorkset.dimension(2) != (int)cell.getDimension() ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): dim 2 (spatial dimension) of cellWorkset array  does not match cell dimension");
  
  
  // Validate whichCell It can be either -1 (default value) or a valid cell ordinal.
  TEST_FOR_EXCEPTION( !( ( (0 <= whichCell ) && (whichCell < cellWorkset.dimension(0) ) ) || (whichCell == -1) ), std::invalid_argument,
                      ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): whichCell = -1 or a valid cell ordinal is required.");
  
  
  // Validate points array: can be rank-2 (P,D) or rank-3 (C,P,D)
  // If rank-2: admissible inCell is rank-1 (P); admissible whichCell is valid cell ordinal but not -1.
  if(physPoints.rank() == 2) {
    
    TEST_FOR_EXCEPTION( (whichCell == -1), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): whichCell = a valid cell ordinal is required with rank-2 input array.");

    TEST_FOR_EXCEPTION( (physPoints.dimension(0) <= 0), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): dim 0 (number of points) >= 1 required for physPoints array ");
    
    TEST_FOR_EXCEPTION( (physPoints.dimension(1) != (int)cell.getDimension() ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): dim 1 (spatial dimension) of physPoints array does not match cell dimension");
    
    // Validate inCell
    TEST_FOR_EXCEPTION( (inCell.rank() != 1), std::invalid_argument, 
                        ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): rank = 1 required for inCell array");
    
    TEST_FOR_EXCEPTION( (inCell.dimension(0) != physPoints.dimension(0)), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): dim 0 (number of points) of inCell array must equal dim 0 of physPoints array");
  }
  // If rank-3: admissible inCell is rank-2 (C,P); admissible whichCell = -1.
  else if (physPoints.rank() == 3){
    
    TEST_FOR_EXCEPTION( !(whichCell == -1), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): whichCell = -1 is required with rank-3 input array.");
    
    TEST_FOR_EXCEPTION( (physPoints.dimension(0) != cellWorkset.dimension(0) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): dim 0 (number of cells)  of physPoints array must equal dim 0 of cellWorkset array ");

    TEST_FOR_EXCEPTION( (physPoints.dimension(1) <= 0), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): dim 1 (number of points) >= 1 required for physPoints array ");
    
    TEST_FOR_EXCEPTION( (physPoints.dimension(2) != (int)cell.getDimension() ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): dim 2 (spatial dimension) of physPoints array does not match cell dimension");
    
    // Validate inCell
    TEST_FOR_EXCEPTION( (inCell.rank() != 2), std::invalid_argument, 
                        ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): rank = 2 required for inCell array");
    
    TEST_FOR_EXCEPTION( (inCell.dimension(0) != physPoints.dimension(0)), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): dim 0 (number of cells) of inCell array must equal dim 0 of physPoints array");    

    TEST_FOR_EXCEPTION( (inCell.dimension(1) != physPoints.dimension(1)), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): dim 1 (number of points) of inCell array must equal dim 1 of physPoints array");    
  }
  else {
    TEST_FOR_EXCEPTION( !( (physPoints.rank() == 2) && (physPoints.rank() ==3) ), std::invalid_argument,
                        ">>> ERROR (Intrepid::CellTools::validateArguments_checkPointwiseInclusion): rank = 2 or 3 required for points array");
  }
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
template<class ArrayCell>
void CellTools<Scalar>::printWorksetSubcell(const ArrayCell &             cellWorkset,
                                            const shards::CellTopology &  parentCell,
                                            const int&                    pCellOrd,
                                            const int&                    subcellDim,
                                            const int&                    subcellOrd,
                                            const int&                    fieldWidth){
  
  // Get the ordinals, relative to reference cell, of subcell cellWorkset
  int subcNodeCount = parentCell.getNodeCount(subcellDim, subcellOrd);
  int pCellDim      = parentCell.getDimension();
  std::vector<int> subcNodeOrdinals(subcNodeCount);
  
  for(int i = 0; i < subcNodeCount; i++){
    subcNodeOrdinals[i] = parentCell.getNodeMap(subcellDim, subcellOrd, i);
  }
  
  // Loop over parent cells and print subcell cellWorkset
  
  std::cout 
    << " Subcell " << subcellOrd << " on parent cell " << pCellOrd << " is " 
    << parentCell.getName(subcellDim, subcellOrd) << " with node(s) \n ({";
  
  for(int i = 0; i < subcNodeCount; i++){
    
    // print Cartesian coordinates of the node
    for(int dim = 0; dim < pCellDim; dim++){
      std::cout
      << std::setw(fieldWidth) << std::right << cellWorkset(pCellOrd, subcNodeOrdinals[i], dim); 
      if(dim < pCellDim - 1){ std::cout << ","; }
    }
    std::cout << "}";
    if(i < subcNodeCount - 1){ std::cout <<", {"; }
  }
  std::cout << ")\n\n";
}



} // namespace Intrepid
#endif
