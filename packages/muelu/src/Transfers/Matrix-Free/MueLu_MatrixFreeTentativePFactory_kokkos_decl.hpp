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
#ifndef MUELU_MATRIXFREETENTATIVEPFACTORY_KOKKOS_DECL_HPP
#define MUELU_MATRIXFREETENTATIVEPFACTORY_KOKKOS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_MatrixFreeTentativePFactory_kokkos_fwd.hpp"

#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#include "Teuchos_ScalarTraits.hpp"

#include "MueLu_Aggregates_kokkos_fwd.hpp"
#include "MueLu_AmalgamationFactory_kokkos_fwd.hpp"
#include "MueLu_AmalgamationInfo_kokkos_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_Utilities_kokkos_fwd.hpp"

namespace MueLu {

  /*!
    @class MatrixFreeTentativePFactory class.
    @brief Factory for building the matrix-free tentative restrictor.

    Factory for creating the matrix-free tentative restrictor.   Nullspace vectors are split across aggregates so that they
    have local support.  The vectors with local support are factored via LAPACK QR.  The Q becomes the
    tentative prolongator, and the R becomes the coarse nullspace.

    Note that the MatrixFreeTentativePFactory also creates the coarse null space vectors, that is, it serves as generating
    factory for the Nullspace vectors on the coarse levels. To feed in the near null space vectors on the finest
    level one has to use the NullspaceFactory.

    @ingroup MueLuTransferClasses

    ## Input/output of MatrixFreeTentativePFactory ##

    ### User parameters of MatrixFreeTentativePFactory ###
    Parameter | type | default | master.xml | validated | requested | description
    ----------|------|---------|:----------:|:---------:|:---------:|------------
     Aggregates         | Factory | null |   | * | * | Generating factory of the aggregates (of type "Aggregates" produced, e.g., by the UncoupledAggregationFactory)
     Nullspace          | Factory | null |   | * | * | Generating factory of the fine nullspace vectors (of type "MultiVector"). In the default case the same instance of the MatrixFreeTentativePFactory is also the generating factory for the null space vectors (on the next coarser levels). Therefore, it is recommended to declare the MatrixFreeTentativePFactory to be the generating factory of the "Nullspace" variable globally using the FactoryManager object! For defining the near null space vectors on the finest level one should use the NullspaceFactory.
     UnAmalgamationInfo | Factory | null |   | * | * | Generating factory of UnAmalgamationInfo. This data (of type "AmalgamationInfo") is produced by the AmalgamationFactory class. The default option should be fine for most standard cases though.
     CoarseMap          | Factory | null |   | * | * | Generating factory of the coarse map. The map generates the coarse domain map of the prolongation operator. The default choice should be fine as long as you do not wish some offset in the domain map (e.g. for block operators in multiphysics problems). The data is generated by the CoarseMapFactory.

    The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
    The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see MatrixFreeTentativePFactory::GetValidParameters).<br>
    The * in the @c requested column states that the data is requested as input with all dependencies (see MatrixFreeTentativePFactory::DeclareInput).

    ### Variables provided by MatrixFreeTentativePFactory ###

    After MatrixFreeTentativePFactory::Build the following data is available (if requested)

    Parameter | generated by | description
    ----------|--------------|------------
    | R       | MatrixFreeTentativePFactory   | Non-smoothed "tentative" prolongation operator (with piece-wise constant transfer operator basis functions)
    | Nullspace | MatrixFreeTentativePFactory | Coarse near null space vectors. Please also check the documentation of the NullspaceFactory for the special dependency tree of the "Nullspace" variable throughout all multigrid levels.
  */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class MatrixFreeTentativePFactory_kokkos;

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  class MatrixFreeTentativePFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > : public PFactory {
  public:
    typedef LocalOrdinal                                             local_ordinal_type;
    typedef GlobalOrdinal                                            global_ordinal_type;
    typedef typename DeviceType::execution_space                     execution_space;
    typedef Kokkos::RangePolicy<local_ordinal_type, execution_space> range_type;
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>      node_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::coordinateType   real_type;

  private:
    // For compatibility
    typedef node_type                                           Node;
#undef MUELU_MATRIXFREETENTATIVEPFACTORY_KOKKOS_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    MatrixFreeTentativePFactory_kokkos() { }

    //! Destructor.
    virtual ~MatrixFreeTentativePFactory_kokkos() { }
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
  };

} //namespace MueLu

#define MUELU_MATRIXFREETENTATIVEPFACTORY_KOKKOS_SHORT
#endif // MUELU_MATRIXFREETENTATIVEPFACTORY_KOKKOS_DECL_HPP
