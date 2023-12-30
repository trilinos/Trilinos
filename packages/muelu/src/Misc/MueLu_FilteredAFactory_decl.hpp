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
#ifndef MUELU_FILTEREDAFACTORY_DECL_HPP
#define MUELU_FILTEREDAFACTORY_DECL_HPP

#include <string>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FilteredAFactory_fwd.hpp"

#include "MueLu_GraphBase.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
namespace MueLu {

/*!
  @class FilteredAFactory class.
  @brief Factory for building filtered matrices using filtered graphs.
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class FilteredAFactory : public SingleLevelFactoryBase {
#undef MUELU_FILTEREDAFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  FilteredAFactory() {}

  //! Destructor.
  virtual ~FilteredAFactory() {}

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
  void BuildReuse(const Matrix& A, const GraphBase& G, const bool lumping, double dirichletThresh, Matrix& filteredA) const;
  void BuildNew(const Matrix& A, const GraphBase& G, const bool lumping, double dirichletThresh, Matrix& filteredA) const;
  void BuildNewUsingRootStencil(const Matrix& A, const GraphBase& G, double dirichletThresh, Level& currentLevel, Matrix& filteredA, bool use_spread_lumping, double DdomAllowGrowthRate, double DdomCap) const;
  void ExperimentalLumping(const Matrix& A, Matrix& filteredA, double rho, double rho2) const;

};  // class FilteredAFactory

}  // namespace MueLu

#define MUELU_FILTEREDAFACTORY_SHORT
#endif  // MUELU_FILTEREDAFACTORY_DECL_HPP
