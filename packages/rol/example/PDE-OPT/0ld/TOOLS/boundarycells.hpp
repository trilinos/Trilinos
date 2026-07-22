// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde.hpp
    \brief Implements the local PDE interface for the Poisson control problem.
*/

#ifndef PDEOPT_TOOLS_BOUNDARY_CELLS_HPP
#define PDEOPT_TOOLS_BOUNDARY_CELLS_HPP

#include "Intrepid_FieldContainer.hpp"
#include "ROL_Ptr.hpp"
#include <vector>

enum EBoundaryType {
  DIRICHLET = 0,
  NEUMANN,
  ROBIN,
  USERDEFINED,
  LAST
};

template <class Real>
class BoundaryCells {
private:
  const int boundaryID_;
  const EBoundaryType boundaryType_;
  const std::vector<int> localCellIndex_;
  const ROL::Ptr<Intrepid::FieldContainer<Real> > cellNodes_;
  const int localCellSideID_;

public:
  BoundaryCells(const int &boundaryID,
                const EBoundaryType &boundaryType,
                const std::vector<int> &localCellIndex,
                const ROL::Ptr<Intrepid::FieldContainer<Real> > &cellNodes,
                const int &localCellSideID)
    : boundaryID_(boundaryID),
      boundaryType_(boundaryType),
      localCellIndex_(localCellIndex),
      cellNodes_(cellNodes),
      localCellSideID_(localCellSideID) {}
  
  const int getBoundaryID(void) const {
    return boundaryID_; 
  }
   
  const EBoundaryType getBoundaryType(void) const {
    return boundaryType_; 
  }
  
  const std::vector<int> getBoundaryCellLocalIndex(void) const {
    return localCellIndex_; 
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real> > getCellNodes(void) const {
    return cellNodes_; 
  }

  const int getLocalCellSideID(void) const {
    return localCellSideID_;
  }
  
};
#endif
