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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_NULLSPACEFACTORY_DECL_HPP
#define MUELU_NULLSPACEFACTORY_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_NullspaceFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"

namespace MueLu {

  /*!
     @class NullspaceFactory class.
     @brief Factory for generating nullspace

      The NullspaceFactory is meant to generate a default approximation for
      the fine level nullspace (Level 0 only). For all other levels it is used
      only to act as "generating factory" for the "Nullspace", which is actually
      handled by the TentativePFactory.

      There are two types of constructors:
      \code{.cpp}
      NullspaceFactory(RCP<const FactoryBase> AFact = Teuchos::null, RCP<const FactoryBase> nullspaceFact = Teuchos::null);
      \endcode
      This constructor uses AFact (or information from FactoyManager) for accessing
      the fine level variable "A" and generates a constant null space approximation
      for A on the finest level (with nPDE constant vectors).
      The nullspaceFact_ defines which factory generates the intermedium and coarse-
      level nullspaces. It must not be Teuchos::null but a TentativePFactory.
      \note If there is a "Nullspace" variable stored on the finest level (given by
      the user) this is is preferred to generating the nullspace from A.

      \code{.cpp}
      NullspaceFactory(std::string nspName, RCP<const FactoryBase> nullspaceFact = Teuchos::null);
      \endcode
      This constructor uses the variable with the variable nspName on the finest level
      as null space for the finest multigrid level.
      The nullspaceFact_ defines which factory generates the intermedium and coarse-
      level nullspaces. It must not be Teuchos::null but a TentativePFactory.

     @ingroup MueLuTransferClasses
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class NullspaceFactory : public SingleLevelFactoryBase {
#undef MUELU_NULLSPACEFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    NullspaceFactory(RCP<const FactoryBase> AFact = Teuchos::null, RCP<const FactoryBase> nullspaceFact = Teuchos::null);

    //! alternative Constructor
    NullspaceFactory(std::string nspName, RCP<const FactoryBase> nullspaceFact = Teuchos::null);

    //! Destructor
    virtual ~NullspaceFactory();

    //@}

    //! @name Input
    //@{
    /*! @brief Specifies the data that this class needs, and the factories that generate that data.

        If the Build method of this class requires some data, but the generating factory is not specified in DeclareInput, then this class
        will fall back to the settings in FactoryManager.
    */

    void DeclareInput(Level &currentLevel) const;

    //@}

    //! @name Build methods.
    //@{

    //! Build an object with this factory.
    void Build(Level &currentLevel) const;

    //@}

  private:

    //! name of nullspace vector on finest level
    std::string nspName_;

    //! Factory for generating matrix A.
    RCP<const FactoryBase> AFact_;

    //! Other nullspace factory (used for level !=0)
    RCP<const FactoryBase> nullspaceFact_;

  }; //class NullspaceFactory

} //namespace MueLu

#define MUELU_NULLSPACEFACTORY_SHORT
#endif // MUELU_NULLSPACEFACTORY_DECL_HPP
