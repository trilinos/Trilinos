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


/** \file   Intrepid_MultiCellDef.hpp
    \brief  Definition file for the Intrepid::MultiCell class.
    \author Created by P. Bochev, D. Ridzal, and D. Day.
*/

#include "Intrepid_CellTemplates.hpp"

namespace Intrepid {

template<class Scalar>
const char* MultiCell<Scalar>::cell_type_names[] = {
"Node",
"Edge",
"Tri",
"Quad",
"Tet",
"Hex",
"Pyramid",
"Prism",
"Polygon",
"Polyhedron",
"CellSet",
"Max_Canonical_Type",
"Polygon1",
"Polygon2",
"Polygon3",
"Polygon4",
"Polygon5",
"Polyhedron1",
"Polyhedron2",
"Polyhedron3",
"Polyhedron4",
"Polyhedron5",
"Max_Custom_Type"
};

template<class Scalar>
MultiCell<Scalar>::MultiCell(const int      num_cells_,
                             const int      ambient_dimension_,
                             const CellType my_cell_type_,
                             const Scalar*  node_coords_) :
// initializations
num_cells(num_cells_),
ambient_dimension(ambient_dimension_),
my_cell_type(my_cell_type_) 
{
  // By default, the Mcell does not compute its pullback information (if pullback is applicable)
  pullback_status = UNDEFINED;
  
  // This constructor does not accept orientation data
  orients_status = UNDEFINED;                             
  
  int nnodes_per_cell = this->getNumNodes(my_cell_type);
  // allocate a chunk of memory first, creating a zero vector of Point vectors
  Point<Scalar> tmp_point(ambient_dimension);
  std::vector< Point<Scalar> > tmp_point_vec(nnodes_per_cell, tmp_point);
  node_coords.assign(num_cells, tmp_point_vec);
  // fill in the coordinate info
  for (int i=0; i<num_cells; i++) {
    for (int j=0; j<nnodes_per_cell; j++) {
      // do some pointer arithmetic
      Point<Scalar> workpoint(node_coords_+(i*nnodes_per_cell+j)*ambient_dimension, ambient_dimension);
      node_coords[i][j] = workpoint;
    }
  }
  
  int num_edges = this->getNumSubcells(my_cell_type, 1);
  // allocate a zero vector of 'short' vectors
  std::vector<short> tmp_edge_vec(num_edges);
  edge_orients.assign(num_cells, tmp_edge_vec);
  
  int num_faces = this->getNumSubcells(my_cell_type, 2);
  // allocate a zero vector of 'short' vectors
  std::vector<short> tmp_face_vec(num_faces);
  face_orients.assign(num_cells, tmp_face_vec);
}

template<class Scalar>
MultiCell<Scalar>::MultiCell(const int      num_cells_,
                             const int      ambient_dimension_,
                             const CellType my_cell_type_,
                             const Scalar*  node_coords_,
                             const short*   edge_orients_,
                             const short*   face_orients_) :
// intializations
num_cells(num_cells_),
ambient_dimension(ambient_dimension_),
my_cell_type(my_cell_type_) 
{
  // By default, the Mcell does not compute its pullback information (if pullback is applicable)
  pullback_status = UNDEFINED;
  
  // This constructor accepts orientation data
  orients_status = DEFINED;                             
  
  int nnodes_per_cell = this->getNumNodes(my_cell_type);
  // allocate a chunk of memory first, creating a zero vector of Point vectors
  Point<Scalar> tmp_point(ambient_dimension);
  std::vector< Point<Scalar> > tmp_point_vec(nnodes_per_cell, tmp_point);
  node_coords.assign(num_cells, tmp_point_vec);
  // fill in the coordinate info
  for (int i=0; i<num_cells; i++) {
    for (int j=0; j<nnodes_per_cell; j++) {
      // do some pointer arithmetic
      Point<Scalar> workpoint(node_coords_+(i*nnodes_per_cell+j)*ambient_dimension, ambient_dimension);
      node_coords[i][j] = workpoint;
    }
  }
  
  // edge orientations
  int num_edges = this->getNumSubcells(my_cell_type, 1);
  // allocate a zero vector of 'short' vectors
  std::vector<short> tmp_edge_vec(num_edges);
  edge_orients.assign(num_cells, tmp_edge_vec);
  // fill in the orientation info
  for (int i=0; i<num_cells; i++) {
    edge_orients[i].assign(edge_orients_+i*num_edges, edge_orients_+(i+1)*num_edges);
  }

  // face orientations
  int num_faces = this->getNumSubcells(my_cell_type, 2);
  // allocate a zero vector of 'short' vectors
  std::vector<short> tmp_face_vec(num_faces);
  face_orients.assign(num_cells, tmp_face_vec);
  // fill in the orientation info
  for (int i=0; i<num_cells; i++) {
    face_orients[i].assign(face_orients_+i*num_faces, face_orients_+(i+1)*num_faces);
  }
}

template<class Scalar>
void MultiCell<Scalar>::setConnMapCustom(const CellType        cell_type_,
                                         const ConnMapTemplate conn_map_template_[]) {
  assert((cell_type_ > MAXCANONICAL) && (cell_type_ < MAXTYPE)); // can only write to a custom template
  // define connectivity map of the new cell type
  conn_map_custom[cell_type_-MAXCANONICAL-1][0] = conn_map_template_[0];
  conn_map_custom[cell_type_-MAXCANONICAL-1][1] = conn_map_template_[1];
  conn_map_custom[cell_type_-MAXCANONICAL-1][2] = conn_map_template_[2];
}


template<class Scalar>
void MultiCell<Scalar>::printMyInfo(std::ostream & out) const {
  out.setf(std::ios_base::scientific, std::ios_base::floatfield);
  out.precision(6);
  
  out  << "\n>>>>>>>>>> MultiCell info: \n\n"; 
  out << std::left;
  out << std::setw(30) << "\t Generating cell type:" << this -> getMyName() << "\n";
  out << std::setw(30) << "\t Pullback status:" << this -> getPullbackInfo() <<"\n";
  out << std::setw(30) << "\t Orientation status:" << StatusNames[orients_status] <<"\n";
  out << std::setw(30) << "\t Number of cells:" << num_cells << "\n\n";
  out << "Nodes:\n";
  int nnodes_per_cell = this->getMyNumNodes();
  for (int i=0; i < num_cells; i++) {
    for (int j=0; j < nnodes_per_cell; j++) {
      out << node_coords[i][j] << "\n";
    }
    out << std::setw(16) << "END CELL" << "\n";
  }
  out << "Edge template:\n";
  int nedges_per_cell = this->getNumSubcells(my_cell_type, 1);
  for (int i=0; i < nedges_per_cell; i++) {
    std::vector<int> tmpnodes;
    this->getSubcellNodes(my_cell_type, 1, i, tmpnodes);
    out << "  " << std::setw(3) << i << " -> {" <<  tmpnodes[0] << ", " << tmpnodes[1] << "}" << "\n";
  }
  out << "Edge orientations:\n";
  for (int i=0; i < num_cells; i++) {
    for (int j=0; j < nedges_per_cell; j++) {
      out << std::setw(5) << edge_orients[i][j];
    }
    out << "\n";
  }
  out << "Face template:\n";
  int nfaces_per_cell = this->getNumSubcells(my_cell_type, 2);
  for (int i=0; i < nfaces_per_cell; i++) {
    std::vector<int> tmpnodes;
    CellType tmptype = this->getSubcellType(my_cell_type, 2, i);
    this->getSubcellNodes(my_cell_type, 2, i, tmpnodes);
    out << "    " << i << " -> " << std::setw(12) << getName(tmptype) << ": ";
    for (int j=0; j < this->getNumNodes(tmptype); j++) {
      out << std::setw(4) << tmpnodes[j];
    }
    out << "\n";
  }
  out << "Face orientations:\n";
  for (int i=0; i < num_cells; i++) {
    for (int j=0; j < nfaces_per_cell; j++) {
      out << std::setw(5) << face_orients[i][j];
    }
    out << "\n";
  }
  out << "<<<<<<<<<< END MULTICELL INFO \n";

}

template<class Scalar>
void MultiCell<Scalar>::setPullback() {
  // First check if pullback was already set
  if(pullback_status == DEFINED) return; 
  // make sure there's enough room to hold all pullbacks to avoid resizing of the vector
  pullback.resize(num_cells);
  for(int cell_id = 0; cell_id < num_cells; cell_id++){
    this -> getPullback(cell_id,pullback[cell_id]);
  }
  pullback_status = DEFINED;
}

template<class Scalar>
void MultiCell<Scalar>::getPullback(const int cell_id_, PullbackTemplate<Scalar>& pullback_) const{
  
  // First check if the pullback has been computed ans set in this instance. If so,
  // copy from the private member (cant's send back non-const reference to private data)
  if(pullback_status == DEFINED) {
    pullback_ = pullback[cell_id_];
    return;
  }
  
  switch(my_cell_type) {
    case EDGE: {
      if(ambient_dimension != 1){
        std::cerr<<"MultiCell::getPullback error: \n"
                   "\t cell dimension does not match ambient dimension\n";
        exit(1);
      }
      Scalar v0 = this->getPoint(cell_id_,0)[0];
      Scalar v1 = this->getPoint(cell_id_,1)[0];
      
      pullback_.cell_type     = my_cell_type;
      pullback_.cell_id       = cell_id_;
      pullback_.pullback_type = AFFINE;
      pullback_.F[0][0] = (v1 - v0)/2.0;
      pullback_.F[0][1] = (v1 + v0)/2.0;
      break;
    }
    case TRI:
      if(ambient_dimension != 2){
        std::cerr<<"MultiCell::getPullback error: \n"
        "\t cell dimension does not match ambient dimension\n";
        exit(1);
      }
      pullback_.cell_type     = my_cell_type;        
      pullback_.cell_id       = cell_id_;
      pullback_.pullback_type = AFFINE;
      for(int dim=0; dim < ambient_dimension; dim++){
        Scalar v0 = this->getPoint(cell_id_,0)[dim];
        Scalar v1 = this->getPoint(cell_id_,1)[dim];
        Scalar v2 = this->getPoint(cell_id_,2)[dim];
          
        pullback_.F[dim][0] = (-v0 + v1);                     // x coefficient
        pullback_.F[dim][1] = (-v0 + v2);                     // y coefficient
        pullback_.F[dim][2] =   v0;                           // free term
      }
      break;
    case QUAD:
      if(ambient_dimension != 2){
        std::cerr<<"MultiCell::getPullback error: \n"
        "\t cell dimension does not match ambient dimension\n";
        exit(1);
      }
      pullback_.cell_type     = my_cell_type;                
      pullback_.cell_id       = cell_id_;
      pullback_.pullback_type = NON_AFFINE;
      for(int dim=0; dim < ambient_dimension; dim++){
        Scalar v0 = this->getPoint(cell_id_,0)[dim];
        Scalar v1 = this->getPoint(cell_id_,1)[dim];
        Scalar v2 = this->getPoint(cell_id_,2)[dim];
        Scalar v3 = this->getPoint(cell_id_,3)[dim];
          
        pullback_.F[dim][0] = ( v0 - v1 + v2 - v3)/4.0;            // xy coefficient
        pullback_.F[dim][1] = (-v0 + v1 + v2 - v3)/4.0;            // x  coefficient
        pullback_.F[dim][2] = (-v0 - v1 + v2 + v3)/4.0;            // y  coefficient
        pullback_.F[dim][3] = ( v0 + v1 + v2 + v3)/4.0;            // free term
      }
      break;
    case TET:
      if(ambient_dimension != 3){
        std::cerr<<"MultiCell::getPullback error: \n"
        "\t cell dimension does not match ambient dimension\n";
        exit(1);
      }
      pullback_.cell_type     = my_cell_type;        
      pullback_.cell_id       = cell_id_;
      pullback_.pullback_type = AFFINE;
      for(int dim=0; dim < ambient_dimension; dim++){
        Scalar v0 = this->getPoint(cell_id_,0)[dim];
        Scalar v1 = this->getPoint(cell_id_,1)[dim];
        Scalar v2 = this->getPoint(cell_id_,2)[dim];
        Scalar v3 = this->getPoint(cell_id_,3)[dim];
        
        pullback_.F[dim][0] = (-v0 + v1);                          // x coefficient
        pullback_.F[dim][1] = (-v0 + v2);                          // y coefficient
        pullback_.F[dim][2] = (-v0 + v3);                          // z coefficient
        pullback_.F[dim][3] = (-v0 + v3);                          // free term
      }
      break;
    case HEX:
      if(ambient_dimension != 3){
        std::cerr<<"MultiCell::getPullback error: \n"
        "\t cell dimension does not match ambient dimension\n";
        exit(1);
      }
      pullback_.cell_type     = my_cell_type;        
      pullback_.cell_id       = cell_id_;
      pullback_.pullback_type = NON_AFFINE;
      for(int dim=0; dim < ambient_dimension; dim++){
        Scalar v0 = this->getPoint(cell_id_,0)[dim];
        Scalar v1 = this->getPoint(cell_id_,1)[dim];
        Scalar v2 = this->getPoint(cell_id_,2)[dim];
        Scalar v3 = this->getPoint(cell_id_,3)[dim];
        
        Scalar v4 = this->getPoint(cell_id_,4)[dim];
        Scalar v5 = this->getPoint(cell_id_,5)[dim];
        Scalar v6 = this->getPoint(cell_id_,6)[dim];
        Scalar v7 = this->getPoint(cell_id_,7)[dim];

        pullback_.F[dim][0] = (-v0 + v1 - v2 + v3 + v4 - v5 + v6 - v7)/8.0;   // xyz
        pullback_.F[dim][1] = ( v0 - v1 + v2 - v3 + v4 - v5 + v6 - v7)/8.0;   // xy
        pullback_.F[dim][2] = ( v0 - v1 - v2 + v3 - v4 + v5 + v6 - v7)/8.0;   // xz
        pullback_.F[dim][3] = ( v0 + v1 - v2 - v3 - v4 - v5 + v6 + v7)/8.0;   // yz
        pullback_.F[dim][4] = (-v0 + v1 + v2 - v3 - v4 + v5 + v6 - v7)/8.0;   // x
        pullback_.F[dim][5] = (-v0 - v1 + v2 + v3 - v4 - v5 + v6 + v7)/8.0;   // y
        pullback_.F[dim][6] = (-v0 - v1 - v2 - v3 + v4 + v5 + v6 + v7)/8.0;   // z
        pullback_.F[dim][7] = ( v0 + v1 + v2 + v3 + v4 + v5 + v6 + v7)/8.0;   // free term
      }
      break;
    case PRISM:
    case PYRAMID:
      if(ambient_dimension != 3){
        std::cerr<<"MultiCell::getPullback error: \n"
        "\t cell dimension does not match ambient dimension\n";
        exit(1);
      }
      std::cout << " Pullback method not implemented...\n";
      exit(1);
      break;
    default:
      std::cerr<<"MultiCell::getPullback error: pullback not defined for this cell type \n";
      exit(1);
      break;
  }
}



template<class Scalar>
LinearMap<Scalar> MultiCell<Scalar>::Jacobian(
  const int                cell_id_,
  const Point<Scalar>&     argument_) const{
  //
  if(pullback_status == DEFINED){                 // First check if pullback was precomputed
    return this -> Jacobian(cell_id_, argument_, pullback[cell_id_]);
  }
  else{                                           // If not, compute the pullback on the fly
    PullbackTemplate<Scalar> temp_pullback;
    this -> getPullback(cell_id_,temp_pullback);
    return this -> Jacobian(cell_id_,argument_,temp_pullback);
  }
  
}



template<class Scalar>
LinearMap<Scalar> MultiCell<Scalar>::Jacobian(
  const int                        cell_id_,
  const Point<Scalar>&             argument_,
  const PullbackTemplate<Scalar>&  pullback_) const{
  //
  // for now, TET and TRI cells have only AFFINE pullbacks implemented
  // First make sure that the pullback corresponds to the cell with cell_id, and
  // that the PointType is REFERENCE if the pullback is not affine! For an
  // affine pullback the Jacobian is constant and the Point argument_ is never
  // used. However, for non-affine pullbacks, the Jacobian is actually evaluated
  // at the argument_ which must be inside the reference cell. 
  //
  if(pullback_.cell_id != cell_id_){
    std::cerr << "Multicell::Jacobian error: \n"
    "\t Cell id does not match the cell id of the pullback\n";
    exit(1);
  }
  if(pullback_.pullback_type == NON_AFFINE && \
     argument_.getType()     == AMBIENT){
    std::cerr << "Multicell::Jacobian error: \n"
    "\t Admissible Point type is REFERENCE, while the argument is AMBIENT type/n";
    exit(1);
  }
  //
  // Temp storage for the linear map coefficients by row for the LinearMap ctor.
  //
  Scalar DF[9];                                   
  switch(my_cell_type) {                         
    case EDGE:  // the EDGE pullback is always affine
      DF[0] = pullback_.F[0][0];                    
      break;
    case TRI:
    case TET:
      for(int row = 0; row < ambient_dimension; row++){
        for(int col = 0; col < ambient_dimension; col++){
          DF[col + row*ambient_dimension] = pullback_.F[row][col]; 
        }
      }
      break;
      // For QUAD and HEX rows contain grad of row-th coordinate function evaluated at argument_
      //
    case QUAD:                                     
      for(int row = 0; row < ambient_dimension; row++){
        DF[0 + row*ambient_dimension] = pullback_.F[row][0]*argument_[1] + pullback_.F[row][1];
        DF[1 + row*ambient_dimension] = pullback_.F[row][0]*argument_[0] + pullback_.F[row][2];
      }
      break;
    case HEX:
      for(int row = 0; row < ambient_dimension; row++){
        DF[0 + row*ambient_dimension] = \
          pullback_.F[row][0]*argument_[1]*argument_[2]+\
          pullback_.F[row][1]*argument_[1]+pullback_.F[row][2]*argument_[2]+pullback_.F[row][4];
        //
        DF[1 + row*ambient_dimension] = \
          pullback_.F[row][0]*argument_[0]*argument_[2]+ \
          pullback_.F[row][1]*argument_[0]+pullback_.F[row][3]*argument_[2]+pullback_.F[row][5];
        //
        DF[2 + row*ambient_dimension] = \
          pullback_.F[row][0]*argument_[0]*argument_[1]+ \
          pullback_.F[row][2]*argument_[0]+pullback_.F[row][3]*argument_[1]+pullback_.F[row][6];
      }
      break;
    default:
      std::cerr << "Multicell::Jacobian error: not implemented\n";
      exit(1);
  }
  LinearMap<Scalar> DF_temp(DF,ambient_dimension);
  return DF_temp;
}



template<class Scalar>
Point<Scalar> MultiCell<Scalar>::Pullback(
  const int            cell_id_, 
  const Point<Scalar>& preimage_) const{
  //
  if(pullback_status == DEFINED){                 // First check if pullback was precomputed
    return this -> Pullback(cell_id_, preimage_, pullback[cell_id_]);
  }
  else{                                           // If not, compute the pullback on the fly
    PullbackTemplate<Scalar> temp_pullback;
    this -> getPullback(cell_id_,temp_pullback);
    return this -> Pullback(cell_id_,preimage_,temp_pullback);
  }
}


template<class Scalar>
Point<Scalar> MultiCell<Scalar>::Pullback(
  const int                       cell_id_, 
  const Point<Scalar>&            preimage_,
  const PullbackTemplate<Scalar>& pullback_) const {
  //
  // First make sure that the pullback corresponds to the cell with cell_id, and
  // that the argument PointType = REFERENCE
  //
  if(pullback_.cell_id != cell_id_){
    std::cerr << " Multicell::Pullback error: \n"
    "\t Cell id does not match the cell id of the pullback\n";
    exit(1);
  }
  if(preimage_.getType() != REFERENCE){
    std::cerr <<  " Multicell::Pullback error: \n"
    "\t Admissible Point type is REFERENCE but argument has AMBIENT type.\n";
    exit(1);
  }
  Scalar x[3];                                    
  switch(my_cell_type){
    case EDGE:
      if(!(-1.0 <= preimage_[0] && preimage_[0] <= 1.0)) {
        std::cerr<< " Multicell::Pullback error: \n"
        " \t The specified 1D Point argument is not in the reference EDGE cell\n";
        exit(1);
      }            
      x[0] = pullback_.F[0][0]*preimage_[0] + pullback_.F[0][1];
      break;
    case TRI:
      if(!((preimage_[0] >= 0) && \
           (preimage_[1] >= 0) && \
           (preimage_[0] + preimage_[1] <= 1))) {
        std::cerr<< " Multicell::Pullback error: \n"
        " \t The specified 2D Point argument is not in the reference TRI cell\n";
        exit(1);
      }      
      for(int dim = 0;dim < ambient_dimension; dim++){
        x[dim] = \
        pullback_.F[dim][0]*preimage_[0] + \
        pullback_.F[dim][1]*preimage_[1] + \
        pullback_.F[dim][2];
      }
      break;
    case QUAD:
      if(!((-1.0 <= preimage_[0] && preimage_[0] <= 1.0) && \
           (-1.0 <= preimage_[1] && preimage_[1] <= 1.0))) {
        std::cerr<< " Multicell::Pullback error: \n"
        " \t The specified 2D Point argument is not in the reference QUAD cell\n";
        exit(1);
      }
      for(int dim = 0; dim < ambient_dimension; dim++){
        x[dim] = \
        pullback_.F[dim][0]*preimage_[0]*preimage_[1] +\
        pullback_.F[dim][1]*preimage_[0] +\
        pullback_.F[dim][2]*preimage_[1] +\
        pullback_.F[dim][3];
      }
      break;
    case TET:
      if(!((preimage_[0] >= 0) && \
           (preimage_[1] >= 0) && \
           (preimage_[2] >= 0) && \
           (preimage_[0] + preimage_[1] + preimage_[2] <= 1))) {
        std::cerr<< " Multicell::Pullback error: \n"
        " \t The specified 3D Point argument is not in the reference TET cell\n";
        exit(1);
      }      
      for(int dim = 0; dim < ambient_dimension; dim++){
        x[dim] = \
        pullback_.F[dim][0]*preimage_[0] + \
        pullback_.F[dim][1]*preimage_[1] + \
        pullback_.F[dim][2]*preimage_[2] + \
        pullback_.F[dim][3];
      }
      break;
    case HEX:
      if(!((-1.0 <= preimage_[0] && preimage_[0] <= 1.0) && \
           (-1.0 <= preimage_[1] && preimage_[1] <= 1.0) && \
           (-1.0 <= preimage_[2] && preimage_[2] <= 1.0))) {
        std::cerr<< " Multicell::Pullback error: \n"
        " \t The specified 3D Point argument is not in the reference HEX cell\n";
        exit(1);
      }
      for(int dim = 0; dim < ambient_dimension; dim++){
        x[dim] = \
        pullback_.F[dim][0]*preimage_[0]*preimage_[1]*preimage_[2]+\
        pullback_.F[dim][1]*preimage_[0]*preimage_[1] + \
        pullback_.F[dim][2]*preimage_[0]*preimage_[2] + \
        pullback_.F[dim][3]*preimage_[1]*preimage_[2] + \
        pullback_.F[dim][4]*preimage_[0] + \
        pullback_.F[dim][5]*preimage_[1] + \
        pullback_.F[dim][6]*preimage_[2] + \
        pullback_.F[dim][7]; 
      }
      break;
    case PRISM:
    case PYRAMID:
      std::cout<<"Not Implemented \n";
      exit(1);
      break;
    default:
      std::cerr<<"MultiCell::Pullback error: unexpected cell type\n";
      exit(1);
  }
  // The return point is in AMBIENT space, set its type accordingly:
  Point<Scalar> temp_point(x,ambient_dimension,AMBIENT);
  return temp_point;
}


template<class Scalar>
Point<Scalar> MultiCell<Scalar>::InversePullback(
  const int            cell_id_, 
  const Point<Scalar>& image_) const{
  //
  if(pullback_status == DEFINED){                 
    // First check if pullback was precomputed
    return this -> InversePullback(cell_id_, image_, pullback[cell_id_]);
  }
  else{                                           
    // If not, compute the pullback on the fly
    PullbackTemplate<Scalar> temp_pullback;
    this -> getPullback(cell_id_,temp_pullback);
    return this -> InversePullback(cell_id_,image_,temp_pullback);
  }
}


template<class Scalar>
Point<Scalar> MultiCell<Scalar>::InversePullback(
  const int                       cell_id_, 
  const Point<Scalar>&            image_,
  const PullbackTemplate<Scalar>& pullback_) const {
  //
  // First make sure that the pullback corresponds to the cell with cell_id
  // and that the point type is AMBIENT:
  //
  if(pullback_.cell_id != cell_id_){
    std::cerr << " Multicell::InversePullback error: \n"
    "\t Cell id does not match the cell id of the pullback.\n";
    exit(1);
  }
  if(image_.getType() != AMBIENT){
    std::cerr <<  " Multicell::InversePullback error: \n"
    "\t Admissible Point type is AMBIENT, argument has REFERENCE type.\n";
    exit(1);
  }
  //
  // Temp storage for the preimage point with the right dimension, do not set type yet
  //
  Point<Scalar> preimage(ambient_dimension);            
  LinearMap<Scalar> A(ambient_dimension);
  Scalar temp[3];                                         
  switch(my_cell_type){
    case EDGE:  // Always affine!
      temp[0] = (image_[0] - pullback_.F[0][1])/pullback_.F[0][0];
      preimage.setData(temp,ambient_dimension);
      preimage.setType(REFERENCE);
      if( preimage.inRefCell(EDGE) == FAILURE) {
        std::cerr<< " Local_0Form::InversePullback error: \n"
        " \t The 1D Point argument is not in the specified EDGE cell\n";
        exit(1);
      }
        return preimage;
      break;
    case TRI:
    case TET:
      //
      // For TRI and TET the affine pullback has the form x = A*x_ref + b and
      // the InversePullback is: preimage = A^{-1}*(image - b). The coefficients 
      // of b are in F[*][ambient_dimension] and A can be obtained by using the 
      // Jacobian method of MultiCell class with any Point argument - it will 
      // not be used because the Jacobian is a constant. Note that for a 
      // non-affine pullback, Jacobian requires point whose PointType = REFERENCE. 
      //
      switch(pullback_.pullback_type){                      
        case AFFINE:                                           
          for(int dim = 0; dim < ambient_dimension; dim++){  
            temp[dim] = pullback_.F[dim][ambient_dimension];             
          }                                                  
          preimage.setData(temp,ambient_dimension); 
          // Compute the preimage:
          preimage = ((this -> Jacobian(cell_id_,image_,pullback_)).\
                        getInverse())*(image_ - preimage);
          // The PointType of the preimage is REFERENCE:
          preimage.setType(REFERENCE);
          // If preimage is not in the reference cell, the specified image was not in the specified cell!
          if( preimage.inRefCell(TRI) == FAILURE) {
            std::cerr<< " Local_0Form::InversePullback error: \n"
            " \t The 2D Point argument is not in the specified TRI cell\n";
            exit(1);
          }
          return preimage;
          break;
        case NON_AFFINE:
          std::cerr << "MultiCell::InversePullback error: method not implemented\n";
          exit(1);
          break;
        default:
          std::cerr << "MultiCell::InversePullback error: unexpected pullback type\n";
          exit(1);
          break;          
      }
      break;
    case QUAD:
    case HEX:
    case PRISM:
    case PYRAMID:
      std::cerr << "MultiCell::InversePullback error: method not implemented\n";
      exit(1);
      break;
    default:
      std::cerr << "MultiCell::InversePullback error: unexpected cell type \n";
      exit(1);
  }
}
template<class Scalar>
std::ostream& operator << (std::ostream& os,
                           const MultiCell<Scalar>& base) {
  base.printMyInfo(os);
  return os;
}

} // namespace Intrepid
