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

/** \file   Intrepid_CellDef.hpp
\brief  Definition file for the Intrepid::Cell class.
\author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {
  
template<class Scalar>
Cell<Scalar>::Cell(const ECell    generatingCellType,
                   const Scalar*  vertices) :
  MultiCell<Scalar>::MultiCell(1, generatingCellType, vertices) {}
    

template<class Scalar>
const Teuchos::Array<short> & Cell<Scalar>::getCellEdgeSigns() const {
  return this -> MultiCell<Scalar>::getCellEdgeSigns(0);
}
    
 
template<class Scalar>
const Teuchos::Array<short> & Cell<Scalar>::getCellFaceSigns() const {
  return this -> MultiCell<Scalar>::getCellFaceSigns(0);
}


template<class Scalar>
const Teuchos::Array<short> & Cell<Scalar>::getCellEdgeTags() const {
  return this -> MultiCell<Scalar>::getCellEdgeTags(0);
}


template<class Scalar>
const Teuchos::Array<short> & Cell<Scalar>::getCellFaceTags() const {
  return this -> MultiCell<Scalar>::getCellFaceTags(0);
}


template<class Scalar>
const Point<Scalar> & Cell<Scalar>::getCellVertex(const int vertexID) const {
  return this->MultiCell<Scalar>::getCellVertex(0, vertexID);
}
    
    
template<class Scalar>
const Teuchos::Array<Point<Scalar> > & Cell<Scalar>::getCellVertices() const {
  return this->MultiCell<Scalar>::getCellVertices(0);
}
    

template<class Scalar>
void Cell<Scalar>::setChart() {
  this -> MultiCell<Scalar>::setChart(0);
}


template<class Scalar>
void Cell<Scalar>::setChart(const ShapePoints<Scalar>& shapePoints) {
  this -> MultiCell<Scalar>::setChart(0, shapePoints);
}


template<class Scalar>
void Cell<Scalar>::setAtlas() {
  this -> MultiCell<Scalar>::setAtlas();
}


template<class Scalar> 
void Cell<Scalar>::setAtlas(const ShapePoints<Scalar>& shapePoints) {
  Teuchos::Array< ShapePoints<Scalar> > dummyArray;
  dummyArray.resize(1);
  dummyArray[0] = shapePoints;
  this -> MultiCell<Scalar>::setAtlas(dummyArray);
}

    
    
template<class Scalar>
Matrix<Scalar> Cell<Scalar>::jacobian(const Point<Scalar>& refPoint,
                                      const double threshold) const {
  return this->MultiCell<Scalar>::jacobian(0, refPoint, threshold);
}


template<class Scalar>
Point<Scalar> Cell<Scalar>::mapToPhysicalCell(const Point<Scalar>& refPoint,
                                              const double threshold) const {
  return this->MultiCell<Scalar>::mapToPhysicalCell(0, refPoint, threshold);
}


template<class Scalar>
Point<Scalar> Cell<Scalar>::mapToReferenceCell(const Point<Scalar>& physPoint) const {
  return this->MultiCell<Scalar>::mapToReferenceCell(0, physPoint);
}
    

template<class Scalar>
bool Cell<Scalar>::inPhysicalCell(const Point<Scalar>& physPoint,
                                      const double threshold) const {
  return this->MultiCell<Scalar>::inPhysicalCell(0, physPoint, threshold);
}


template<class Scalar>
const EMapping Cell<Scalar>::getCellMappingType() const {
  return this -> MultiCell<Scalar>::getCellMappingType(0);
}

  
} // namespace Intrepid
