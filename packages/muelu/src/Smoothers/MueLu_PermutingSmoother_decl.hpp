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
/*
 * MueLu_PermutingSmoother_decl.hpp
 *
 *  Created on: Nov 28, 2012
 *      Author: wiesner
 */

#ifndef MUELU_PERMUTINGSMOOTHER_DECL_HPP
#define MUELU_PERMUTINGSMOOTHER_DECL_HPP

#include <Teuchos_ParameterList.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_TrilinosSmoother_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"
#include "MueLu_PermutationFactory_fwd.hpp"

namespace MueLu {

  /*!
    @class PermutingSmoother
    @brief This class first calculates row- and column permutation operators and applies a smoother to the permuted linear system.

  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class PermutingSmoother : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {
#undef MUELU_PERMUTINGSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors / destructors
    //@{

    //! @brief Constructor
    //!    @param[in] mapName   Name of map object in level class, which rows/cols can be permuted
    //!    @param[in] mapFact   generating factory of map with name mapName
    //!    @param[in] type      string that contains type of smoother (e.g. "RELAXATION" or "ILU")
    //!    @param[in] paramList parameter list with parameters for smoother (default: empty)
    //!    @param[in] overlap   LocalOrdinal with overlap inforation (default: 0)
    //!    @param[in] permFact  factory, generating permutation and scaling matrices (default: Teuchos::null -> use internal PermutationFactory instance)
    PermutingSmoother(std::string const & mapName, const RCP<const FactoryBase> & mapFact, std::string const & type = "", Teuchos::ParameterList const & paramList = Teuchos::ParameterList(), LO const & overlap = 0, RCP<FactoryBase> permFact = Teuchos::null);

    //! Destructor
    virtual ~PermutingSmoother();
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    //! @name Setup and Apply methods.
    //@{

    //! @brief Set up the direct solver.
    void Setup(Level &currentLevel);

    /*! @brief Apply the direct solver.
    Solves the linear system <tt>AX=B</tt> using the constructed solver.
    @param X initial guess
    @param B right-hand side
    @param InitialGuessIsZero This option has no effect.
    */
    void Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero=false) const;
    //@}

    RCP<SmootherPrototype> Copy() const;

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

    //@}

  private:
    //! ifpack1/2-specific key phrase that denote smoother type
    std::string type_;

    //! parameter list that is used by Ifpack/Ifpack2 internally
    Teuchos::ParameterList paramList_;

    //! overlap when using the smoother in additive Schwarz mode
    LO overlap_;

    //! Permutation Factory
    RCP<FactoryBase> permFact_;

    //! permQT matrix object
    RCP<Matrix> permQT_;

    //! permP matrix object
    RCP<Matrix> permP_;

    //! scaling matrix object
    Teuchos::RCP<Matrix> diagScalingOp_;

    //
    // Underlying Smoother
    //

    //! Smoother
    RCP<SmootherPrototype> s_; // TrilinosSmoother object

  }; // class PermutingSmoother

} // namespace MueLu

#define MUELU_PERMUTINGSMOOTHER_SHORT
#endif /* MUELU_PERMUTINGSMOOTHER_DECL_HPP */
