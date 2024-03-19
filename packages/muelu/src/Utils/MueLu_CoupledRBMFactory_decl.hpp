// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2013 Sandia Corporation
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
#ifndef MUELU_COUPLEDRBMFACTORY_DECL_HPP
#define MUELU_COUPLEDRBMFACTORY_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_CoupledRBMFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"

namespace MueLu {

/*!
  @class CoupledRBMFactory
  @ingroup MueLuTransferClasses
  @brief Nullspace Factory for coupled acoustic-elastic problems.

   Combines standard nullspace with rigid body modes.
   Assumes that acoustic pressure DOFs are padded with 2 extra DOFs
   (so that there are 3 DOFs at each mesh grid point)
*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class CoupledRBMFactory : public SingleLevelFactoryBase {
#undef MUELU_COUPLEDRBMFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{
  //! Constructor
  CoupledRBMFactory(const int numPDEs)
    : nspName_("Nullspace")
    , numPDEs_(numPDEs) {}
  //! Constructor
  CoupledRBMFactory(const std::string& nspName = "Nullspace")
    : nspName_(nspName)
    , numPDEs_(3) {}

  //! Destructor.
  virtual ~CoupledRBMFactory();

  //@}

  //! @name Input
  //@{

  /*! @brief Specifies the data that this class needs, and the factories that generate that data.

      If the Build method of this class requires some data, but the generating factory is not specified in DeclareInput, then this class
      will fall back to the settings in FactoryManager.
  */
  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  //! Build an object with this factory.
  void Build(Level& currentLevel) const;

  void BuildRBM(RCP<Matrix>& A, RCP<MultiVector>& Coords, RCP<MultiVector>& nullspace) const;

  //@}
  void setNumPDEs(int numPDEs) {
    numPDEs_ = numPDEs;
  }

  void setLastAcousticDOF(int lastDOF) {
    lastAcousticDOF_ = lastDOF;
  }

 private:
  //! name of nullspace vector on finest level
  std::string nspName_;

  int numPDEs_;

  int lastAcousticDOF_;

};  // class CoupledRBMFactory

}  // namespace MueLu

#define MUELU_COUPLEDRBMFACTORY_SHORT
#endif  // MUELU_COUPLEDRBMFACTORY_DECL_HPP
