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
#ifndef MUELU_GENERALGEOMETRICPFACTORY_DECL_HPP
#define MUELU_GENERALGEOMETRICPFACTORY_DECL_HPP

#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_GeneralGeometricPFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"

namespace MueLu {


/*!
  @class GenearlGeometricPFactory
  @ingroup MueLuTransferClasses
  @brief Prolongator factory performing geometric coarsening.

  The geometric algorithm assumes the underlying mesh is reasonably structured.
  Any rate of coarsening can be applied, and the rate is automatically decrease along
  an edge if the number of element is not divisible by the coarsening rate.
  The coarsening rate is allowed to be different in all direction which means that
  semi-coarsening can be achieved within this algorithm in 1 or 2 directions.
  The main difficulty is to obtain the number of elements/nodes in each directions
  to identify coarse nodes and fine nodes.

  ## Input/output of GeneralGeometricPFactory ##

  ### User parameters of SemiCoarsenPFactory ###
  | Parameter   | type    | default   | master.xml | validated | requested | description                                                                      |
  | ------------|---------|-----------|:----------:|:---------:|:---------:|----------------------------------------------------------------------------------|
  | Coarsen     | string  | -1        |            |           |           | A string that specify the coarsening rate, if it is a single character, it will  |
  |             |         |           |            |           |           | indicate a unique coarsening rate in each direction, if it is longer, it will be | 
  |             |         |           |            |           |           | processed as a vector with 3 entries, one for each spatial direction             |
  | A           | Factory | null      |            | *         | *         | Generating factory of the matrix A used during the prolongator smoothing process |
  | Nullspace   | Factory | null      |            | *         | *         | Generating factory of the nullspace. The GeneralGeometricPFactory provides       |
  |             |         |           |            |           |           | a coarse version of the given Nullspace.                                         |
  | Coordinates | Factory | NoFactory |            | *         | *         | Generating factory for coorindates. The coordinates are expected to be provided  |
  |             |         |           |            |           |           | on the finest level using the NoFactory mechanism. The coordinates are used to   |
  |             |         |           |            |           |           | compute the coarsening stencil and coarse coordinates are generated for the next |
  |             |         |           |            |           |           | level.                                                                           |



  The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
  The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see GeneralGeometricCoarsenPFactory::GetValidParameters).<br>
  The * in the @c requested column states that the data is requested as input with all dependencies (see GeneralGeometricCoarsenPFactory::DeclareInput).

  ### Variables provided by GeneralGeometricPFactory ###

  After GeneralGeometricPFactory::Build the following data is available (if requested)

  | Parameter         | generated by             | description
  |-------------------|--------------------------|------------------------------------------------------------------------------------------------------------------|
  | P                 | GeneralGeometricPFactory | Prolongator                                                                                                      |
  | Nullspace         | GeneralGeometricPFactory | Coarse nullspace (the fine level nullspace information is coarsened using P to generate a coarse version         |
  |                   |                          | of the nullspace. No scaling is applied.                                                                         |
  | coarseCoordinates | NoFactory                | Coarse coordinates that will be used on the next fine level to compute the coarsening stencils                   |

*/
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class GeneralGeometricPFactory : public PFactory {
#undef MUELU_GENERALGEOMETRICPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    GeneralGeometricPFactory() { }

    //! Destructor.
    virtual ~GeneralGeometricPFactory() { }
    //@}

    RCP<const ParameterList> GetValidParameterList() const;

    //! Input
    //@{

    void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

    //@}

    //! @name Build methods.
    //@{

    void Build (Level& fineLevel, Level& coarseLevel) const;
    void BuildP(Level& fineLevel, Level& coarseLevel) const;

    //@}

  private:
    LO MakeGeneralGeometricP(LO const ndm, LO const fNodePerDir[], LO const cNodePerDir[], LO const coarse_rate[], LO const end_rate[],
                             const RCP<Xpetra::MultiVector<double,LO,GO,NO> >& fCoords, LO const nnzP, LO const dofsPerNode, RCP<const Map>& coarseMap,
                             RCP<Matrix> & Amat, RCP<Matrix>& P, RCP<Xpetra::MultiVector<double,LO,GO,NO> >& cCoords) const;

    void ComputeStencil(const LO ndm, const LO ordering, const SC coord[9][3], SC stencil[8]) const;

  }; //class GeneralGeometricPFactory

} //namespace MueLu

#define MUELU_GENERALGEOMETRICPFACTORY_SHORT
#endif // MUELU_GENERALGEOMETRICPFACTORY_DECL_HPP
