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
/*
 * MueLu_SchurComplementFactory_decl.hpp
 *
 *  Created on: Jun 18, 2012
 *      Author: wiesner
 */

#ifndef MUELU_SCHURCOMPLEMENTFACTORY_DECL_HPP_
#define MUELU_SCHURCOMPLEMENTFACTORY_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_MapExtractor_fwd.hpp>
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_StridedMap_fwd.hpp>


#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_FactoryManagerBase_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_SmootherBase_fwd.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Utilities_fwd.hpp"


namespace MueLu {

  /*!
    @class SchurComplementFactory class.
    @brief Factory for building the Schur Complement for a 2x2 block matrix.

    For a blocked matrix
        A = [A_00  A_01; A_10  A_11]
    it computes the Schur complement
        S = A_11 - 1/\omega A_10 \hat{A_00}^{-1} A_01,
    where \omega is some scaling factor and \hat{A_00} an approximation of A_00 (just the diagonal of A_00)
  */
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class SchurComplementFactory : public SingleLevelFactoryBase {
#undef MUELU_SCHURCOMPLEMENTFACTORY_SHORT
    #include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    SchurComplementFactory() { }

    //! Destructor.
    virtual ~SchurComplementFactory() { }
    //@}

    //! Input
    //@{

    void DeclareInput(Level& currentLevel) const;

    RCP<const ParameterList> GetValidParameterList() const;

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level& currentLevel) const;

    //@}


  private:

  }; // class SchurComplementFactory

} // namespace MueLu

#define MUELU_SCHURCOMPLEMENTFACTORY_SHORT
#endif /* MUELU_SCHURCOMPLEMENTFACTORY_DECL_HPP_ */
