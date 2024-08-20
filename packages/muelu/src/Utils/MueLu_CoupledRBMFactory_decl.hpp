// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
