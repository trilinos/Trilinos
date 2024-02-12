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
#ifndef MUELU_TENTATIVEPFACTORY_DECL_HPP
#define MUELU_TENTATIVEPFACTORY_DECL_HPP

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialQRDenseSolver.hpp>

#include <Xpetra_CrsGraphFactory_fwd.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_Import_fwd.hpp>
#include <Xpetra_ImportFactory_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TentativePFactory_fwd.hpp"

#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_AmalgamationFactory_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"
#include "MueLu_CoarseMapFactory_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

/*!
  @class TentativePFactory class.
  @brief Factory for building tentative prolongator.

  Factory for creating tentative prolongator.   Nullspace vectors are split across aggregates so that they
  have local support.  The vectors with local support are factored via LAPACK QR.  The Q becomes the
  tentative prolongator, and the R becomes the coarse nullspace.

  Note that the TentativePFactory also creates the coarse null space vectors, that is, it serves as generating
  factory for the Nullspace vectors on the coarse levels. To feed in the near null space vectors on the finest
  level one has to use the NullspaceFactory.

  @ingroup MueLuTransferClasses

  ## Input/output of TentativePFactory ##

  ### User parameters of TentativePFactory ###
  Parameter | type | default | master.xml | validated | requested | description
  ----------|------|---------|:----------:|:---------:|:---------:|------------
   A                  | Factory | null |   | * | * | Generating factory of the matrix A
   Aggregates         | Factory | null |   | * | * | Generating factory of the aggregates (of type "Aggregates" produced, e.g., by the UncoupledAggregationFactory)
   Nullspace          | Factory | null |   | * | * | Generating factory of the fine nullspace vectors (of type "MultiVector"). In the default case the same instance of the TentativePFactory is also the generating factory for the null space vectors (on the next coarser levels). Therefore, it is recommended to declare the TentativePFactory to be the generating factory of the "Nullspace" variable globally using the FactoryManager object! For defining the near null space vectors on the finest level one should use the NullspaceFactory.
   UnAmalgamationInfo | Factory | null |   | * | * | Generating factory of UnAmalgamationInfo. This data (of type "AmalgamationInfo") is produced by the AmalgamationFactory class. The default option should be fine for most standard cases though.
   CoarseMap          | Factory | null |   | * | * | Generating factory of the coarse map. The map generates the coarse domain map of the prolongation operator. The default choice should be fine as long as you do not wish some offset in the domain map (e.g. for block operators in multiphysics problems). The data is generated by the CoarseMapFactory.

  The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
  The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see TentativePFactory::GetValidParameters).<br>
  The * in the @c requested column states that the data is requested as input with all dependencies (see TentativePFactory::DeclareInput).

  ### Variables provided by TentativePFactory ###

  After TentativePFactory::Build the following data is available (if requested)

  Parameter | generated by | description
  ----------|--------------|------------
  | P       | TentativePFactory   | Non-smoothed "tentative" prolongation operator (with piece-wise constant transfer operator basis functions)
  | Nullspace | TentativePFactory | Coarse near null space vectors. Please also check the documentation of the NullspaceFactory for the special dependency tree of the "Nullspace" variable throughout all multigrid levels.
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class TentativePFactory : public PFactory {
#undef MUELU_TENTATIVEPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  TentativePFactory() {}

  //! Destructor.
  virtual ~TentativePFactory() {}
  //@}

  RCP<const ParameterList> GetValidParameterList() const;

  //! Input
  //@{

  void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{

  void Build(Level& fineLevel, Level& coarseLevel) const;
  void BuildP(Level& fineLevel, Level& coarseLevel) const;

  //@}

 private:
  void BuildPuncoupled(RCP<Matrix> A, RCP<Aggregates> aggregates, RCP<AmalgamationInfo> amalgInfo, RCP<MultiVector> fineNullspace,
                       RCP<const Map> coarseMap, RCP<Matrix>& Ptentative, RCP<MultiVector>& coarseNullspace, const int levelID) const;
  void BuildPcoupled(RCP<Matrix> A, RCP<Aggregates> aggregates, RCP<AmalgamationInfo> amalgInfo, RCP<MultiVector> fineNullspace,
                     RCP<const Map> coarseMap, RCP<Matrix>& Ptentative, RCP<MultiVector>& coarseNullspace) const;
  void BuildPuncoupledBlockCrs(RCP<Matrix> A, RCP<Aggregates> aggregates, RCP<AmalgamationInfo> amalgInfo, RCP<MultiVector> fineNullspace,
                               RCP<const Map> coarseMap, RCP<Matrix>& Ptentative, RCP<MultiVector>& coarseNullspace, const int levelID) const;

  mutable bool bTransferCoordinates_ = false;

};  // class TentativePFactory

}  // namespace MueLu

// TODO: noQR_

#define MUELU_TENTATIVEPFACTORY_SHORT
#endif  // MUELU_TENTATIVEPFACTORY_DECL_HPP
