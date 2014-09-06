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
#ifndef MUELU_SCHWARZSMOOTHER_DECL_HPP
#define MUELU_SCHWARZSMOOTHER_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2) && defined(HAVE_MUELU_IFPACK2)

#include <Teuchos_ParameterList.hpp>

#include "MueLu_SchwarzSmoother_fwd.hpp"

#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Factory_decl.hpp"
#include "Ifpack2_Factory_def.hpp"

#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace Amesos2 { template<class OP, class MV> class Solver; }

namespace MueLu {

  /*!
    @class SchwarzSmoother
    @brief Class that uses Amesos2 direct solvers and Ifpack2 preconditioners in an additive schwarz setting.
  */

  template <class Scalar = SmootherPrototype<>::scalar_type,
            class LocalOrdinal = typename SmootherPrototype<Scalar>::local_ordinal_type,
            class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
            class Node = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class SchwarzSmoother : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node>
  {
#undef MUELU_SCHWARZSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors / destructors
    //@{

    /*! @brief Constructor
      Creates a MueLu interface to the direct solvers in Amesos2 and preconditioners in Ifpack2.
    */
    SchwarzSmoother(const std::string& type = "", const Teuchos::ParameterList& paramList = Teuchos::ParameterList(), const LocalOrdinal& overlapLevel = 0);

    //! Destructor
    virtual ~SchwarzSmoother();
    //@}

    //! Input
    //@{

    void DeclareInput(Level& currentLevel) const;

    //@}

    //! @name Setup and Apply methods.
    //@{

    /*! @brief Set up the local subdomain solver or preconditioner.
      This creates the underlying Amesos2 solver object or Ifpack2 preconditioning object according
      to the parameter list options passed into the SchwarzSmoother constructor.
    */
    void Setup(Level& currentLevel);

    /*! @brief Apply Additive Schwarz with Amesos2/Ifpack2 as the local subdomain solver.
    @param X initial guess
    @param B right-hand side
    @param InitialGuessIsZero This option has no effect.
    */
    void Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero = false) const;
    //@}

    RCP<SmootherPrototype> Copy() const;

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    void print(Teuchos::FancyOStream& out, const VerbLevel verbLevel = Default) const;

    //@}

  private:
    typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Tpetra_CrsMatrix;
    typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> Tpetra_MultiVector;

    //! local smoother type
    std::string type_;

    //! parameter list
    Teuchos::ParameterList paramList_;

    //! level of overlap
    LocalOrdinal overlapLevel_;

    //! pointer to maps
    Teuchos::RCP< Tpetra::Map<LO,GO,NO> > localRowMap_;
    Teuchos::RCP< Tpetra::Export<LO,GO,NO> > TpetraExporter_;
    Teuchos::RCP< Tpetra::Import<LO,GO,NO> > TpetraImporter_;
    Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > UniqueMap_;
    Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > OverlapMap_;

    //! pointer to matrix
    RCP<Matrix> A_;

    //! pointer to Amesos2 solver object
    RCP<Amesos2::Solver<Tpetra_CrsMatrix, Tpetra_MultiVector> > prec_;

    //! pointer to Ifpack2 preconditiioner
    RCP< Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > ifpack2prec_;

  }; // class SchwarzSmoother

} // namespace MueLu

#define MUELU_SCHWARZSMOOTHER_SHORT
#endif // HAVE_MUELU_TPETRA && HAVE_MUELU_AMESOS2
#endif // MUELU_SCHWARZSMOOTHER_DECL_HPP

// TODO: PARAMETER LIST NOT TAKE INTO ACCOUNT !!!
