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
#ifndef MUELU_HIERARCHYFACTORY_HPP
#define MUELU_HIERARCHYFACTORY_HPP

#include "Teuchos_RCP.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

#include "MueLu_Hierarchy_fwd.hpp"

namespace MueLu {

//!
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class HierarchyFactory : public BaseClass {
#undef MUELU_HIERARCHYFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //@{ Constructors/Destructors.

  //! Destructor.
  virtual ~HierarchyFactory() {}

  //@}

  //! Create an empty Hierarchy object
  // Note: This function is not very useful at the moment as MueLu only have on Hierarchy class.
  //       In the future, we might have an abstract Hierarchy class and several derived Hierarchy classes.
  //       Using this function will then be the recommended way to generate a Hierarchy.
  //
  // This method is called Create() instead of Build(), because it return an non-initialized
  // object (ie: MG setup is not done).
  // Build() function in MueLu returns initialized objects.
  virtual RCP<Hierarchy> CreateHierarchy() const = 0;

  //! Create a labeled empty Hierarchy object
  virtual RCP<Hierarchy> CreateHierarchy(const std::string& label) const = 0;

  //! Setup Hierarchy object
  virtual void SetupHierarchy(Hierarchy& H) const = 0;

};  // class HierarchyFactoryBase

}  // namespace MueLu

#define MUELU_HIERARCHYFACTORY_SHORT
#endif  // ifndef MUELU_HIERARCHYFACTORY_HPP
