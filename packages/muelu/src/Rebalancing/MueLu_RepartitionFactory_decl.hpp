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
#ifndef MUELU_REPARTITIONFACTORY_DECL_HPP
#define MUELU_REPARTITIONFACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#ifdef HAVE_MPI

// Some classes are only used in the definition (_def.hpp) of this class
// but forward declarations are needed here to enable the UseShortNames mechanism.
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Import_fwd.hpp>
#include <Xpetra_ImportFactory_fwd.hpp>
#include <Xpetra_Export_fwd.hpp>
#include <Xpetra_ExportFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_RepartitionFactory_fwd.hpp"
#include "MueLu_CloneRepartitionInterface_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"

namespace MueLu {

/*!
  @class RepartitionFactory class.
  @brief Factory for building permutation matrix that can be be used to shuffle data (matrices, vectors) among processes

  This factory acts on both the number of partitions and a vector (usually created by Zoltan) that indicates to which partitions
  the current level's system matrix's DOFS belong.

  We always call the Interface routines in ZoltanInterface, Zoltan2Interface or IsorropiaInterface, but they have special short-cuts
  if no rebalancing is necessary (then return "Partition = null") or all data is moved to one processor.

  ## Input/output of RepartitionInterface ##

  ### User parameters of RepartitionInterface ###
  Parameter | type | default | master.xml | validated | requested | description
  ----------|------|---------|:----------:|:---------:|:---------:|------------
  | A                                      | Factory | null  |   | * | * | Generating factory of the matrix A used during the prolongator smoothing process |
  | number of partitions                   | Factory | null  |   | * | * | Factory calculating the "number of partitions" for rebalancing (usually an instance of RepartitionHeuristicFactory)
  | Partition                              | Factory | null  |   | * | * | Factory generating the "Partition" variable containing the distribution of the DOFs over the reduced number of processors (given by "number of partitions")
  | repartition: print partition distribution | bool | false | * | * |   | Partition distribution printout.
  | repartition: remap parts                  | bool | false | * | * |   | Postprocessing for partitioning to reduce data migration.
  | repartition: remap num values             | int  |  4    | * | * |   | Number of maximum components from each processor used to construct partial bipartite graph.

  The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
  The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see RepartitionFactory::GetValidParameters).<br>
  The * in the @c requested column states that the data is requested as input with all dependencies (see RepartitionFactory::DeclareInput).

  ### Variables provided by RepartitionInterface ###

  After RepartitionFactory::Build the following data is available (if requested)

  Parameter | generated by | description
  ----------|--------------|------------
  | Importer | RepartitionFactory   | Importer used for rebalancing

  Importer contains the Xpetra::Import object to rebalance the data (matrix or vectors). It is Teuchos::null if no rebalancing is necessary.

  @note For blocked systems you need one RepartionFactory for each sub block together with a separate Instance of an Interface class.

*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class RepartitionFactory : public SingleLevelFactoryBase {
#undef MUELU_REPARTITIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  RepartitionFactory() {}

  //! Destructor.
  virtual ~RepartitionFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! @name Input
  //@{

  /*! @brief Determines the data that RepartitionFactory needs, and the factories that generate that data.

      If this class requires some data, but the generating factory is not specified in DeclareInput, then this class
      will fall back to the settings in FactoryManager.
  */
  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  //! Build an object with this factory.
  void Build(Level& currentLevel) const;

  //@}

  //! @name Helper methods.
  //@{

  /*!
    @brief Determine which process should own each partition.

    Partitions are assigned to processes in order to minimize data movement.  The basic idea is that a good choice for partition
    owner is to choose the pid that already has the greatest number of nonzeros for a particular partition. If willAcceptPartition=false
    is set on a rank, no partition will be placed there.

  */
  void DeterminePartitionPlacement(const Matrix& A, GOVector& decomposition, GO numPartitions, bool willAcceptPartition = true, bool allSubdomainsAcceptPartitions = true) const;

};  // class RepartitionFactory

}  // namespace MueLu

#define MUELU_REPARTITIONFACTORY_SHORT

#endif  // ifdef HAVE_MPI
#endif  // MUELU_REPARTITIONFACTORY_DECL_HPP
