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
#ifndef MUELU_COALESCEDROPFACTORY_KOKKOS_DECL_HPP
#define MUELU_COALESCEDROPFACTORY_KOKKOS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#include "Xpetra_Matrix_fwd.hpp"

#include "MueLu_CoalesceDropFactory_kokkos_fwd.hpp"

#include "MueLu_AmalgamationInfo_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_LWGraph_kokkos_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

/*!
  @class CoalesceDropFactory_kokkos
  @brief Factory for creating a graph based on a given matrix.

  Factory for creating graphs from matrices with entries selectively dropped.
  This factory combines the functionality of CoalesceDropFactory and FilteredAFactory from the non-Kokkos
  code path.

  For an in-depth discussion, see https://github.com/trilinos/Trilinos/issues/1676.

  ## Code paths ##

  Both the classic dropping strategy as well as a coordinate-based distance
  laplacian method is implemented. For performance reasons there are four
  distinctive code paths for the classical method:

  - one DOF per node without dropping (i.e. "aggregation: drop tol" = 0.0)
  - one DOF per node with dropping (i.e. "aggregation: drop tol" > 0.0)
  - DOFs per node > 1 withouth dropping
  - DOFs per node > 1 with dropping

  Additionally there is a code path for the distance-laplacian mode.

  ## Input/output of CoalesceDropFactory_kokkos ##

  ### User parameters of CoalesceDropFactory_kokkos ###
  Parameter | type | default | master.xml | validated | requested | description
  ----------|------|---------|:----------:|:---------:|:---------:|------------
   A              | Factory | null |   | * | * | Generating factory of the operator A
   UnAmalgamationInfo        | Factory | null |   | * | * | Generating factory of type AmalgamationFactory which generates the variable 'UnAmalgamationInfo'. Do not change the default unless you know what you are doing.
   Coordinates               | Factory | null |   | * | (*) | Generating factory for variable 'Coordinates'. The coordinates are only needed if "distance laplacian" is chosen for the parameter "aggregation: drop scheme"
   "aggregation: drop scheme" | std::string | "classical" | * | * |   | Coalescing algorithm. You can choose either "classical" (=default) or "distance laplacian"
   "aggregation: drop tol" | double | 0.0 | * | * |   | Threshold parameter for dropping small entries
   "aggregation: Dirichlet threshold" | double | 0.0 | * | * |   | Threshold for determining whether entries are zero during Dirichlet row detection
   "lightweight wrap" | bool | true |   |   *  |   | hidden switch between fast implementation based on MueLu::LWGraph and a failsafe slower implementation based on Xpetra::Graph (for comparison). The user should not change the default value (=true)

  The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
  The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see CoalesceDropFactory_kokkos::GetValidParameters).<br>
  The * in the @c requested column states that the data is requested as input with all dependencies (see CoalesceDropFactory_kokkos::DeclareInput).

  ### Variables provided by UncoupledAggregationFactory ###

  After CoalesceDropFactory_kokkos::Build the following data is available (if requested)

  Parameter | generated by | description
  ----------|--------------|------------
  Graph   | CoalesceDropFactory_kokkos   | Graph of matrix A
  DofsPerNode | CoalesceDropFactory_kokkos | number of DOFs per node. Note, that we assume a constant number of DOFs per node for all nodes associated with the operator A.

  ## Amalgamation process ##

  The CoalesceDropFactory_kokkos is internally using the AmalgamationFactory
  for amalgamating the dof-based maps to node-based maps. The
  AmalgamationFactory creates the "UnAmalgamationInfo" container which
  basically stores all the necessary information for translating dof based
  data to node based data and vice versa. The container is used, since this
  way the amalgamation is only done once and later reused by other factories.

  Of course, often one does not need the information from the
  "UnAmalgamationInfo" container since the same information could be
  extracted of the "Graph" or the map from the "Coordinates" vector.
  However, there are also some situations (e.g. when doing rebalancing based
  on HyperGraph partitioning without coordinate information) where one has
  not access to a "Graph" or "Coordinates" variable.
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class CoalesceDropFactory_kokkos;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
class CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosDeviceWrapperNode<DeviceType> > : public SingleLevelFactoryBase {
 public:
  using local_ordinal_type  = LocalOrdinal;
  using global_ordinal_type = GlobalOrdinal;
  using execution_space     = typename DeviceType::execution_space;
  using range_type          = Kokkos::RangePolicy<local_ordinal_type, execution_space>;
  using node_type           = Tpetra::KokkosCompat::KokkosDeviceWrapperNode<DeviceType>;

 private:
  // For compatibility
  using Node = node_type;
#undef MUELU_COALESCEDROPFACTORY_KOKKOS_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  CoalesceDropFactory_kokkos() {}

  //! Destructor
  virtual ~CoalesceDropFactory_kokkos() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  void Build(Level& currentLevel) const;
};

}  // namespace MueLu

#define MUELU_COALESCEDROPFACTORY_KOKKOS_SHORT
#endif  // MUELU_COALESCEDROPFACTORY_KOKKOS_DECL_HPP
