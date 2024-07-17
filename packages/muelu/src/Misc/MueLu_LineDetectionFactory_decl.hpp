// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_LINEDETECTIONFACTORY_DECL_HPP
#define MUELU_LINEDETECTIONFACTORY_DECL_HPP

// same as in SemiCoarsenPFactory (TODO rework this)
#define VERTICAL 1
#define HORIZONTAL 2
#define GRID_SUPPLIED -1

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_LineDetectionFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

namespace MueLu {

/*!
  @class LineDetectionFactory class.
  @brief Factory for building line detection information
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class LineDetectionFactory : public SingleLevelFactoryBase {
#undef MUELU_LINEDETECTIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  using coordinate_type       = typename Teuchos::ScalarTraits<SC>::coordinateType;
  using CoordinateMultiVector = typename Xpetra::MultiVector<coordinate_type, LO, GO, NO>;

  //! @name Constructors/Destructors.
  //@{

  LineDetectionFactory()
    : Zorientation_(VERTICAL) {}

  //! Destructor.
  virtual ~LineDetectionFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*!
    @brief Build method.

    Builds line detection information and stores it in currentLevel
    */
  void Build(Level& currentLevel) const;

  //@}

 private:
  void sort_coordinates(LO numCoords, LO* OrigLoc,
                        coordinate_type* xvals,
                        coordinate_type* yvals,
                        coordinate_type* zvals,
                        coordinate_type* xtemp,
                        coordinate_type* ytemp,
                        coordinate_type* ztemp,
                        bool flipXY = false) const;

  LO ML_compute_line_info(LO LayerId[], LO VertLineId[],
                          LO Ndof, LO DofsPerNode,
                          LO MeshNumbering, LO NumNodesPerVertLine,
                          coordinate_type* xvals, coordinate_type* yvals, coordinate_type* zvals,
                          const Teuchos::Comm<int>& comm) const;

  void ML_az_dsort2(coordinate_type dlist[], LO N, LO list2[]) const;

  //! internally stores line detection mode
  //! can be either vertical, horizontal or coordinates
  //! for the first run. On the coarser levels we automatically
  //! switch to vertical mode
  mutable LO Zorientation_;

};  // class LineDetectionFactory

}  // namespace MueLu

#define MUELU_LINEDETECTIONFACTORY_SHORT
#endif  // MUELU_LINEDETECTIONFACTORY_DECL_HPP
