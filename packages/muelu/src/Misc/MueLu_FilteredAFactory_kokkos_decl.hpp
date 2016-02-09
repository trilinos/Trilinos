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
#ifndef MUELU_FILTEREDAFACTORY_KOKKOS_DECL_HPP
#define MUELU_FILTEREDAFACTORY_KOKKOS_DECL_HPP

#include <string>

#include "MueLu_ConfigDefs.hpp"
#ifdef HAVE_MUELU_KOKKOS_REFACTOR

#include "MueLu_FilteredAFactory_kokkos_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

namespace MueLu {

  /*!
    @class FilteredAFactory_kokkos class.
    @brief Factory for building filtered matrices using filtered graphs.
  */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class FilteredAFactory_kokkos;

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  class FilteredAFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > : public SingleLevelFactoryBase {
  public:
    typedef LocalOrdinal                                        local_ordinal_type;
    typedef GlobalOrdinal                                       global_ordinal_type;
    typedef typename DeviceType::execution_space                execution_space;
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> node_type;

  private:
    // For compatibility
    typedef node_type                                           Node;
#undef MUELU_FILTEREDAFACTORY_KOKKOS_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{

    FilteredAFactory_kokkos() { }

    //! Destructor.
    virtual ~FilteredAFactory_kokkos() { }

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

      Builds filtered matrix and returns it in <tt>currentLevel</tt>.
      */
    void Build(Level& currentLevel) const;

    //@}
  private:
    void BuildReuse(const Matrix& A, const LWGraph_kokkos& G, const bool lumping, Matrix& filteredA) const;
    void BuildNew  (const Matrix& A, const LWGraph_kokkos& G, const bool lumping, Matrix& filteredA) const;

  }; //class FilteredAFactory_kokkos

} //namespace MueLu

#define MUELU_FILTEREDAFACTORY_KOKKOS_SHORT
#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_FILTEREDAFACTORY_KOKKOS_DECL_HPP
