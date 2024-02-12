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

#ifndef MUELU_COARSENINGVISUALIZATIONFACTORY_DECL_HPP_
#define MUELU_COARSENINGVISUALIZATIONFACTORY_DECL_HPP_

#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_VisualizationHelpers.hpp"
#include "MueLu_CoarseningVisualizationFactory_fwd.hpp"

namespace MueLu {

class Level;

/*!
  @class CoarseningVisualizationFactory class.
  @brief Factory to visualize coarsening information using prolongation operators

  @ingroup MueLuVisualizationClasses

  ## Input/output of CoarseningVisualizationFactory ##

  ### User parameters of CoarseningVisualizationFactory ###
  Parameter | type | default | master.xml | validated | requested | description
  ----------|------|---------|:----------:|:---------:|:---------:|------------
  | visualization: output filename           | string  | viz%LEVELID  |  | * |   | filename for VTK-style visualization output |
  | visualization: output file: time step    | int     | 0 |  | * |   | time step (overwrites '%TIMESTEP' in output file name) |
  | visualization: output file: iter         | int     | 0 |  | * |   | nonlinear iteration (overwrites '%ITER' in output file name) |
  | visualization: style                     | string  | Point Cloud |   | * |  | style of aggregation visualization for VTK output. Can be either "Point Cloud", "Jacks", or "Convex Hulls" |
  | visualization: fine graph edges | bool | false  |   | * |  | Draw fine node connections in VTK output (only works for 1 dofs per node!) |
  | visualization: build colormap | bool | false  |   | * |  | Output a random color map for paraView in a separate xml file. |
  | P | Factory | Teuchos::null  |   | * | * | Prolongator factory. The user has to declare either P or Ptent but not both at the same time. |
  | Ptent | Factory | Teuchos::null  |   | * | * | Tentative prolongator factory. The user has to declare either P or Ptent but not both at the same time. |
  | Coordinates | Factory | Teuchos::null  |   | * | * | Factory for Coordinates on fine level |
  | Graph | Factory | Teuchos::null  |   | * | * | Factory for Graph of A |

  The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
  The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see CoarseningVisualizationFactory::GetValidParameters).<br>
  The * in the @c requested column states that the data is requested as input with all dependencies (see CoarseningVisualizationFactory::DeclareInput).
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class CoarseningVisualizationFactory : public TwoLevelFactoryBase, public VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_COARSENINGVISUALIZATIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  CoarseningVisualizationFactory() {}

  //! Destructor.
  virtual ~CoarseningVisualizationFactory() {}
  //@}

  RCP<const ParameterList> GetValidParameterList() const;

  //! Input
  //@{

  void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

  //@}

  //@{
  //! @name Build methods.

  //! Build an object with this factory.
  void Build(Level &fineLevel, Level &coarseLevel) const;

  //@}

};  // class CoarseningVisualizationFactory
}  // namespace MueLu

#define MUELU_COARSENINGVISUALIZATIONFACTORY_SHORT

#endif /* MUELU_COARSENINGVISUALIZATIONFACTORY_DECL_HPP_ */
