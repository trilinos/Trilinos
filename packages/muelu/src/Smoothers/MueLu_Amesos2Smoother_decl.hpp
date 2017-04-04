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
#ifndef MUELU_AMESOS2SMOOTHER_DECL_HPP
#define MUELU_AMESOS2SMOOTHER_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2)

#include <Teuchos_ParameterList.hpp>

#include "MueLu_Amesos2Smoother_fwd.hpp"

#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace Amesos2 { template<class OP, class MV> class Solver; }

namespace MueLu {

  /*!
    @class Amesos2Smoother
    @ingroup MueLuSmootherClasses
    @brief Class that encapsulates Amesos2 direct solvers.

    This class creates an Amesos2 preconditioner factory.  The factory is capable of generating direct solvers
    based on the type and ParameterList passed into the constructor.  See the constructor for more information.
  */

  template <class Scalar = SmootherPrototype<>::scalar_type,
            class LocalOrdinal = typename SmootherPrototype<Scalar>::local_ordinal_type,
            class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
            class Node = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class Amesos2Smoother : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
#undef MUELU_AMESOS2SMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors / destructors
    //@{

    /*! @brief Constructor
      Creates a MueLu interface to the direct solvers in the Amesos2 package.
      If you are using type=="", then either SuperLU or KLU2 are used by default.
    */
    Amesos2Smoother(const std::string& type = "", const Teuchos::ParameterList& paramList = Teuchos::ParameterList());

    //! Destructor
    virtual ~Amesos2Smoother();
    //@}

    //! Input
    //@{

    void DeclareInput(Level& currentLevel) const;

    //@}

    //! @name Setup and Apply methods.
    //@{

    /*! @brief Set up the direct solver.
      This creates the underlying Amesos2 solver object according to the parameter list options passed into the
      Amesos2Smoother constructor.  This includes doing a numeric factorization of the matrix.
    */
    void Setup(Level& currentLevel);

    /*! @brief Apply the direct solver.
    Solves the linear system <tt>AX=B</tt> using the constructed solver.
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
   
    //! Get a rough estimate of cost per iteration
    size_t getNodeSmootherComplexity() const;

    //@}

  private:
    typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Tpetra_CrsMatrix;
    typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> Tpetra_MultiVector;

    //! amesos2-specific key phrase that denote smoother type
    std::string type_;

    //! pointer to Amesos2 solver object
    RCP<Amesos2::Solver<Tpetra_CrsMatrix, Tpetra_MultiVector> > prec_;

    bool useTransformation_;
    RCP<MultiVector> X_, B_;

  }; // class Amesos2Smoother

#ifdef HAVE_MUELU_EPETRA

# if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
  // Stub specialization for missing Epetra template args
  template<>
  class Amesos2Smoother<double,int,int,Xpetra::EpetraNode> : public SmootherPrototype<double,int,int,Xpetra::EpetraNode> {
    typedef double              Scalar;
    typedef int                 LocalOrdinal;
    typedef int                 GlobalOrdinal;
    typedef Xpetra::EpetraNode  Node;
#undef MUELU_AMESOS2SMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    Amesos2Smoother(const std::string& type = "", const Teuchos::ParameterList& paramList = Teuchos::ParameterList()) {
      MUELU_TPETRA_ETI_EXCEPTION("Amesos2Smoother<double,int,int,EpetraNode>","Amesos2Smoother<double,int,int,EpetraNode>","int");
    }
    virtual ~Amesos2Smoother() {};
    void DeclareInput(Level& currentLevel) const {};
    void Setup(Level& currentLevel) {};
    void Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero = false) const {};

    RCP<SmootherPrototype> Copy() const { return Teuchos::null; };

    std::string description() const { return std::string(""); };
    void print(Teuchos::FancyOStream& out, const VerbLevel verbLevel = Default) const {};
   
    //! Get a rough estimate of cost per iteration
    size_t getNodeSmootherComplexity() const {size_t cplx=0; return cplx;};
  };
# endif
#endif // HAVE_MUELU_EPETRA

} // namespace MueLu

#define MUELU_AMESOS2SMOOTHER_SHORT
#endif // HAVE_MUELU_TPETRA && HAVE_MUELU_AMESOS2
#endif // MUELU_AMESOS2SMOOTHER_DECL_HPP
