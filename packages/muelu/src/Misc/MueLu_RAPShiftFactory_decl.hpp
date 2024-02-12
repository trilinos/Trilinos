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
#ifndef MUELU_RAPSHIFTFACTORY_DECL_HPP
#define MUELU_RAPSHIFTFACTORY_DECL_HPP

#include <string>

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_RAPShiftFactory_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"

namespace MueLu {
/*!
  @class RAPShiftFactory
  @brief Factory for building coarse grid matrices, when the matrix
         is of the form K+a*M. Useful when you want to change the shift
         variable ("a") at every level. Each level must store the stiffness
         matrix K and mass matrix M separately.
*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class RAPShiftFactory : public TwoLevelFactoryBase {
#undef MUELU_RAPSHIFTFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  RAPShiftFactory();

  virtual ~RAPShiftFactory() {}

  //@}

  //! @name Input
  //@{

  RCP<const ParameterList> GetValidParameterList() const;

  void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{
  void Build(Level &fineLevel, Level &coarseLevel) const;
  //@}

  //! @name Handling of user-defined transfer factories
  //@{

  //! Indicate that the restriction operator action should be implicitly defined by the transpose of the prolongator.
  void SetImplicitTranspose(bool const &implicit) {
    implicitTranspose_ = implicit;
  }

  void SetShifts(std::vector<Scalar> &shifts) {
    shifts_.clear();
    shifts_ = shifts;
  }

  //@}

  //@{
  /*! @brief Add transfer factory in the end of list of transfer factories in RepartitionAcFactory.

  Transfer factories are derived from TwoLevelFactoryBase and project some data from the fine level to
  the next coarser level.
  */
  void AddTransferFactory(const RCP<const FactoryBase> &factory);

  // TODO add a function to remove a specific transfer factory?

  //! Returns number of transfer factories.
  size_t NumTransferFactories() const { return transferFacts_.size(); }

  //@}

 private:
  //! If true, the action of the restriction operator action is implicitly defined by the transpose of the prolongator.
  bool implicitTranspose_;

  //! list of user-defined transfer Factories
  std::vector<RCP<const FactoryBase> > transferFacts_;

  // vector of shifting terms
  std::vector<Scalar> shifts_;

};  // class RAPShiftFactory

}  // namespace MueLu

#define MUELU_RAPSHIFTFACTORY_SHORT
#endif  // MUELU_RAPSHIFTFACTORY_DECL_HPP
