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
#ifndef MUELU_MHDRAPFACTORY_DECL_HPP
#define MUELU_MHDRAPFACTORY_DECL_HPP

#include <string>

#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
//#include "MueLu_MHDRAPFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"

namespace MueLu {
/*!
  @class MHDRAPFactory
  @brief Factory for building coarse matrices.
*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class MHDRAPFactory : public TwoLevelFactoryBase {
#undef MUELU_MHDRAPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  MHDRAPFactory();

  virtual ~MHDRAPFactory() {}

  //@}

  //! @name Input
  //@{

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

  //@}

  //! @name internal print methods.
  //    static std::string PerfUtils::PrintMatrixInfo(const Matrix & Ac, const std::string & msgTag);

  static std::string PrintLoadBalancingInfo(const Matrix &Ac, const std::string &msgTag);

 private:
  //! If true, the action of the restriction operator action is implicitly defined by the transpose of the prolongator.
  bool implicitTranspose_;

};  // class MHDRAPFactory

}  // namespace MueLu

#define MUELU_MHDRAPFACTORY_SHORT
#endif  // MUELU_MHDRAPFACTORY_DECL_HPP
