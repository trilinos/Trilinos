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

#ifndef MUELU_REGIONRFACTORY_KOKKOS_DECL_HPP
#define MUELU_REGIONRFACTORY_KOKKOS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Types.hpp"

#ifdef HAVE_MUELU_KOKKOS_REFACTOR
#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include <Kokkos_Core.hpp>

#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_RegionRFactory_kokkos_fwd.hpp"

namespace MueLu {

  /*!
    @class RegionRFactory_kokkos class
    @brief Factory that builds a restriction operator for region multigrid
  */

  template <class Scalar = DefaultScalar,
            class LocalOrdinal = DefaultLocalOrdinal,
            class GlobalOrdinal = DefaultGlobalOrdinal,
            class Node = DefaultNode>
  class RegionRFactory_kokkos : public TwoLevelFactoryBase {
#undef MUELU_REGIONRFACTORY_KOKKOS_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    using real_type = typename Teuchos::ScalarTraits<SC>::coordinateType;
    using realvaluedmultivector_type = typename Xpetra::MultiVector<real_type, LO, GO, NO>;
    using execution_space = typename Node::execution_space;
    using memory_space    = typename Node::memory_space;
    using device_type     = Kokkos::Device<execution_space, memory_space>;
    using intTupleView    = typename Kokkos::View<int[3], device_type>;
    using LOTupleView     = typename Kokkos::View<LO[3],  device_type>;

    //! @name Constructors/Destructors.
    //@{

    //! Default Constructor
    RegionRFactory_kokkos() = default;

    //!Destructor
    virtual ~RegionRFactory_kokkos() = default;
    //@}

    //! Input
    //@{
    RCP<const ParameterList> GetValidParameterList() const;

    void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

    //@}

    //! @name Build methods.
    //@{

    void Build(Level& fineLevel, Level& coarseLevel) const;

    void Build3D(const int numDimensions,
                 Array<LO>& lFineNodesPerDim,
                 const RCP<Matrix>& A,
                 const RCP<realvaluedmultivector_type>& fineCoordinates,
                 RCP<Matrix>& R,
                 RCP<realvaluedmultivector_type>& coarseCoordinates,
                 Array<LO>& lCoarseNodesPerDim) const;

    //@}

  }; // class RegionRFactory_kokkos

} // namespace MueLu

#define MUELU_REGIONRFACTORY_KOKKOS_SHORT
#endif //ifdef HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_REGIONRFACTORY_KOKKOS_DECL_HPP
