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
#ifndef MUELU_DIRECTSOLVER_DECL_HPP
#define MUELU_DIRECTSOLVER_DECL_HPP

#include <Teuchos_ParameterList.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_DirectSolver_fwd.hpp"

#include "MueLu_SmootherPrototype.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Amesos2Smoother_fwd.hpp"
#include "MueLu_AmesosSmoother_fwd.hpp"

// Note: DirectSolver is a SmootherPrototype that cannot be turned into a smoother using Setup().
//       When this prototype is cloned using Copy(), the clone is an Amesos or an Amesos2 smoother.
//       The clone can be used as a smoother after calling Setup().

namespace MueLu {

  /*!
    @class DirectSolver
    @brief Class that encapsulates direct solvers. Autoselection of AmesosSmoother or Amesos2Smoother according to the compile time configuration of Trilinos
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class DirectSolver : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {
#undef MUELU_DIRECTSOLVER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors / destructors
    //@{

    //! @brief Constructor
    //! Note: only parameters shared by Amesos and Amesos2 should be used for type and paramList (example: type= "Klu", "Superlu", paramList = <empty>) .
    DirectSolver(std::string const & type = "", Teuchos::ParameterList const & paramList = Teuchos::ParameterList(), RCP<FactoryBase> AFact = Teuchos::null);

    //! Destructor
    virtual ~DirectSolver() { }

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    //! @name Setup and Apply methods.
    //@{

    //! DirectSolver cannot be turned into a smoother using Setup(). Setup() always returns a RuntimeError exception.
    void Setup(Level &currentLevel);

    //! DirectSolver cannot be applied. Apply() always returns a RuntimeError exception.
    void Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero=false) const;

    //@}

    //! When this prototype is cloned using Copy(), the clone is an Amesos or an Amesos2 smoother.
    RCP<SmootherPrototype> Copy() const;

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const;

    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

    //@}

  private:

    //
    // Parameters
    //

    //! amesos1/2-specific key phrase that denote smoother type
    std::string type_;

    //! parameter list that is used by Amesos internally
    Teuchos::ParameterList paramList_;

    //! A Factory
    RCP<FactoryBase> AFact_;

    //
    // Underlying Smoother
    //

    //! Smoother
    RCP<SmootherPrototype> s_;

  }; // class DirectSolver

} // namespace MueLu

#define MUELU_DIRECTSOLVER_SHORT
#endif // MUELU_DIRECTSOLVER_DECL_HPP
