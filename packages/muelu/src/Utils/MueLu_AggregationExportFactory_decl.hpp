// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * MueLu_AggregationExportFactory_decl.hpp
 *
 *  Created on: Feb 10, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AGGREGATIONEXPORTFACTORY_DECL_HPP_
#define MUELU_AGGREGATIONEXPORTFACTORY_DECL_HPP_

#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_VisualizationHelpers.hpp"
#include "MueLu_AggregationExportFactory_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Graph_fwd.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_AmalgamationFactory_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"

namespace MueLu {

class Level;

/*!
  @class AggregationExportFactory class.
  @brief Factory to export aggregation info or visualize aggregates using VTK

  Note, that some routines only work for 1 dof per node.

  @ingroup MueLuVisualizationClasses

  ## Input/output of AggregationExportFactory ##

  ### User parameters of AggregationExportFactory ###
  Parameter | type | default | master.xml | validated | requested | description
  ----------|------|---------|:----------:|:---------:|:---------:|------------
  | aggregation: output filename           | string  |   |  | * |   | filename for VTK-style visualization output |
  | aggregation: output file: time step    | int     | 0 |  | * |   | time step (overwrites '%TIMESTEP' in output file name) |
  | aggregation: output file: iter         | int     | 0 |  | * |   | nonlinear iteration (overwrites '%ITER' in output file name) |
  | aggregation: output file: agg style    | string  | Point Cloud |   | * |  | style of aggregation visualization for VTK output. Can be either "Point Cloud", "Jacks", or "Convex Hulls" |
  | aggregation: output file: fine graph edges | bool | false  |   | * |  | Draw fine node connections in VTK output (only works for 1 dofs per node!) |
  | aggregation: output file: build colormap | bool | false  |   | * |  | Output a random color map for paraView in a separate xml file. |
  | Output filename | string |   |    | * |  | Output file name for aggregation data export (outdated, do not use) |
  | Output file: time step | int | 0  |   | * |  | time step variable for output filename (outdated, do not use) |
  | Output file: iter      | int | 0  |   | * |  | nonlinear iteration variable for output filename (outdated, do not use) |
  | A | Factory | Teuchos::null  |   | * | * | Factory for A |
  | Coordinates | Factory | Teuchos::null  |   | * | * | Factory for Coordinates (only necessary for vtk output) |
  | Graph | Factory | Teuchos::null  |   | * | * | Factory for Graph of A (only necessary for vtk output) |
  | Aggregates | Factory | Teuchos::null  |   | * | * | Factory for Aggregates |
  | UnAmalgamationInfo | Factory | Teuchos::null  |   | * | * | Factory for UnAmalgamationInfo |
  | DofsPerNode | Factory | Teuchos::null  |   | * | * | Factory for DofsPerNode |

  The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
  The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see AggregationExportFactory::GetValidParameters).<br>
  The * in the @c requested column states that the data is requested as input with all dependencies (see AggregationExportFactory::DeclareInput).
*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class AggregationExportFactory : public TwoLevelFactoryBase, public VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_AGGREGATIONEXPORTFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  AggregationExportFactory()
    : doFineGraphEdges_(false)
    , doCoarseGraphEdges_(false)
    , numNodes_(0)
    , numAggs_(0)
    , dims_(0)
    , myRank_(-1)
    , aggsOffset_(0) {}

  //! Destructor.
  virtual ~AggregationExportFactory() {}
  //@}

  RCP<const ParameterList> GetValidParameterList() const;

  //! Input
  //@{

  void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

  //@}

  //@{
  //! @name Build methods.

  //! Build an object with this factory.
  void Build(Level& fineLevel, Level& coarseLevel) const;

  using coordinate_type       = typename Teuchos::ScalarTraits<SC>::coordinateType;
  using CoordinateMultiVector = typename Xpetra::MultiVector<coordinate_type, LO, GO, NO>;

  //@}

 private:
  // Break different viz styles into separate functions for organization:
  void doJacksPlus_(std::vector<int>& vertices, std::vector<int>& geomSizes) const;
  void doConvexHulls(std::vector<int>& vertices, std::vector<int>& geomSizes) const;
  void doGraphEdges_(std::ofstream& fout, Teuchos::RCP<Matrix>& A, Teuchos::RCP<GraphBase>& G, bool fine, int dofs) const;  // add geometry to display node connections from a matrix. Connections in graph but not matrix have different color.

  // write VTK data
  void writeFile_(std::ofstream& fout, std::string styleName, std::vector<int>& vertices, std::vector<int>& geomSizes) const;
  void buildColormap_() const;
  void writePVTU_(std::ofstream& pvtu, std::string baseFname, int numProcs) const;

  static const int CONTRAST_1_ = -1;
  static const int CONTRAST_2_ = -2;
  static const int CONTRAST_3_ = -3;

  // Data that the different styles need to have available when building geometry
  mutable Teuchos::RCP<CoordinateMultiVector> coords_;        // fine local coordinates
  mutable Teuchos::RCP<CoordinateMultiVector> coordsCoarse_;  // coarse local coordinates
  mutable Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId_;
  mutable Teuchos::ArrayRCP<LocalOrdinal> aggSizes_;
  mutable std::vector<bool> isRoot_;
  mutable bool doFineGraphEdges_;
  mutable bool doCoarseGraphEdges_;
  mutable int numNodes_;
  mutable int numAggs_;
  mutable int dims_;
  mutable int myRank_;
  mutable Teuchos::RCP<const Map> nodeMap_;        // map used in A and Coordinates to map local ordinals to global ordinals. Need the whole map especially if it's not contiguous.
  mutable Teuchos::RCP<const Map> nodeMapCoarse_;  // Map for Ac
  mutable int aggsOffset_;                         // in a global list of aggregates, the offset of local aggregate indices
};                                                 // class AggregationExportFactory
}  // namespace MueLu

#define MUELU_AGGREGATIONEXPORTFACTORY_SHORT

#endif /* MUELU_AGGREGATIONEXPORTFACTORY_DECL_HPP_ */
