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
#ifndef MUELU_SEMICOARSENPFACTORY_KOKKOS_DECL_HPP
#define MUELU_SEMICOARSENPFACTORY_KOKKOS_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_SemiCoarsenPFactory_kokkos_fwd.hpp"

namespace MueLu {

/*!
  @class SemiCoarsenPFactory_kokkos
  @ingroup MueLuTransferClasses
  @brief Prolongator factory performing semi-coarsening.

  The semi-coarsening is performed along user-provided "vertical lines" (in z-direction).
  The line detection algorithm can be found in the LineDetectionFactory.
  Usually, the SemiCoarsenPFactory_kokkos is used together with the TogglePFactory and a
  second TentativePFactory_kokkos which allows to dynamically switch from semi-coarsening to
  aggregation-based coarsening (or any other compatible coarsening algorithm).

  ## Input/output of SemiCoarsenPFactory_kokkos ##

  ### User parameters of SemiCoarsenPFactory_kokkos ###
  Parameter | type | default | master.xml | validated | requested | description
  ----------|------|---------|:----------:|:---------:|:---------:|------------
  | A                         | Factory | null      | | * | * | Generating factory of the matrix A used during the prolongator smoothing process |
  | Nullspace                 | Factory | null      | | * | * | Generating factory of the nullspace. The SemiCoarsenPFactory_kokkos provides a coarse version of the given Nullspace. |
  | Coordinates               | Factory | NoFactory | | * | * | Generating factory for coorindates. The coordinates are expected to be provided on the finest level using the NoFactory mechanism. The coordinates are used to determine the number of z-layers if not otherwise provided by the user. |
  | LineDetection_VertLineIds | Factory | null      | | * | * | Generating factory for LineDetection information. Usually provided by the LineDetectionFactory. Array with vertical line ids for all nodes on current processor. |
  | LineDetection_Layers      | Factory | null      | | * | * | Generating factory for LineDetection information. Usually provided by the LineDetectionFactory. Array with layer id for all nodes on current processor. |
  | CoarseNumZLayers          | Factory | null      | | * | * | Generating factory for LineDetection information. Usually provided by the LineDetectionFactory. Number of remaining z-layers after semi-coarsening. |
  | semicoarsen: coarsen rate | int     | null      | | * | * | Coarsening rate along vertical lines (2 corresponds to classical semicoarsening. Values > 2 for more aggressive coarsening). |

  The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
  The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see SemiCoarsenPFactory_kokkos::GetValidParameters).<br>
  The * in the @c requested column states that the data is requested as input with all dependencies (see SemiCoarsenPFactory_kokkos::DeclareInput).

  ### Variables provided by SemiCoarsenPFactory_kokkos ###
  After SemiCoarsenPFactory_kokkos::Build the following data is available (if requested)

  Parameter | generated by | description
  ----------|--------------|------------
  | P          | SemiCoarsenPFactory_kokkos | Prolongator
  | Nullspace  | SemiCoarsenPFactory_kokkos | Coarse nullspace (the fine level nullspace information is coarsened using P to generate a coarse version of the nullspace. No scaling is applied.
  | NumZLayers | NoFactory                  | Number of z layers after coarsening. Necessary input for LineDetectionFactory. Useful input for TogglePFactory.
*/

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class SemiCoarsenPFactory_kokkos;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal,
          class DeviceType>
class SemiCoarsenPFactory_kokkos<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>> : public PFactory {
public:
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef typename DeviceType::execution_space execution_space;
  typedef Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> node_type;

private:
  // For compatibility
  typedef node_type Node;

#undef MUELU_SEMICOARSENPFACTORY_KOKKOS_SHORT
#include "MueLu_UseShortNames.hpp"

public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  SemiCoarsenPFactory_kokkos() : bTransferCoordinates_(false) {}

  //! Destructor.
  virtual ~SemiCoarsenPFactory_kokkos() = default;
  //@}

  RCP<const ParameterList> GetValidParameterList() const;

  //! Input
  //@{

  void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{

  void Build(Level &fineLevel, Level &coarseLevel) const;
  void BuildP(Level &fineLevel, Level &coarseLevel) const;

  //@}

  void BuildSemiCoarsenP(Level &coarseLevel, const LO NFRows, const LO NFNodes,
                         const LO DofsPerNode, const LO NFLayers,
                         const LO NCLayers, const ArrayRCP<LO> LayerId,
                         const ArrayRCP<LO> VertLineId, const RCP<Matrix> &Amat,
                         const RCP<MultiVector> fineNullspace, RCP<Matrix> &P,
                         RCP<MultiVector> &coarseNullspace) const;

  // boolean which is true if coordinate information is available to be
  // transferred to coarse coordinate information
  mutable bool bTransferCoordinates_;
}; // class SemiCoarsenPFactory_kokkos

} // namespace MueLu

#define MUELU_SEMICOARSENPFACTORY_KOKKOS_SHORT
#endif // MUELU_SEMICOARSENPFACTORY_KOKKOS_DECL_HPP
