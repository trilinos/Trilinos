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
Cell<Scalar>::Cell(const ECell    generatingCellType,
                   const Scalar*  vertices,
                   const short*   subcellSigns,
                   const int      subcellDim) :
  MultiCell<Scalar>::MultiCell(1, generatingCellType, vertices, subcellSigns, subcellDim) {}
    
    
template<class Scalar>
Cell<Scalar>::Cell(const ECell    generatingCellType,
                   const Scalar*  vertices,
                   const short*   edgeSigns,
                   const short*   faceSigns) :
  MultiCell<Scalar>::MultiCell(1, generatingCellType, vertices, edgeSigns, faceSigns) {}
    

template<class Scalar>
const Teuchos::Array<short> & Cell<Scalar>::getMySubcellSigns(const int subcellDim) const {
  return this->MultiCell<Scalar>::getMySubcellSigns(0, subcellDim);
}
    
    
template<class Scalar>
const Point<Scalar> & Cell<Scalar>::getVertex(const int vertexID) const {
  return this->MultiCell<Scalar>::getVertex(0, vertexID);
}
    
    
template<class Scalar>
const Teuchos::Array<Point<Scalar> > & Cell<Scalar>::getCell() const {
  return this->MultiCell<Scalar>::getCell(0);
}
    

template<class Scalar>
const ChartTemplate<Scalar> & Cell<Scalar>::getChart() const {
  return this->MultiCell<Scalar>::getChart(0);
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
EFailCode Cell<Scalar>::insidePhysicalCell(const Point<Scalar>& physPoint,
                                           const double threshold) const {
  return this->MultiCell<Scalar>::insidePhysicalCell(0, physPoint, threshold);
}
  
} // namespace Intrepid
