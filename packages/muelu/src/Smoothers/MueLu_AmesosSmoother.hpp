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
#ifndef MUELU_AMESOSSMOOTHER_HPP
#define MUELU_AMESOSSMOOTHER_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_AMESOS)

#include <Teuchos_ParameterList.hpp>
#include <Xpetra_Matrix_fwd.hpp>

class Epetra_LinearProblem;
class Amesos_BaseSolver;

#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Exceptions.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

  /*!
    @class AmesosSmoother
    @brief Class that encapsulates Amesos direct solvers.

    This class creates an Amesos preconditioner factory.  The factory is capable of generating direct solvers
    based on the type and ParameterList passed into the constructor.  See the constructor for more information.
  */

  class AmesosSmoother : public SmootherPrototype<double, int, int>
  {
    typedef double Scalar;
    typedef int    LocalOrdinal;
    typedef int    GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;
    typedef Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps LocalMatOps;
#undef MUELU_AMESOSSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors / destructors
    //@{

    /*! @brief Constructor

        Creates a MueLu interface to the direct solvers in the Amesos package.  The options are those specified in
        the Amesos user's manual.

        @param type solver type
        @param list options for the particular solver type

        Here is how to select the more commonly used direct solvers:

        - KLU (serial sparse direct solver)
            - <tt>type</tt> = <tt>Amesos-KLU</tt>
            - parameter list options
                - none required

        - SuperLU (serial sparse super-nodal direct solver)
            - <tt>type</tt> = <tt>Amesos-SuperLU</tt>
            - parameter list options
                - none required

        If you are using type=="", then either SuperLU or KLU are used by default.

        See also Amesos_Klu and Amesos_Superlu.

    */

    AmesosSmoother(std::string const & type = "", Teuchos::ParameterList const & paramList = Teuchos::ParameterList());

    //! Destructor
    virtual ~AmesosSmoother();

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    //! @name Setup and Apply methods.
    //@{

    /*! @brief Set up the direct solver.
      This creates the underlying Amesos solver object according to the parameter list options passed into the
      AmesosSmoother constructor.  This includes doing a numeric factorization of the matrix.
    */
    void Setup(Level &currentLevel);

    /*! @brief Apply the direct solver.

        Solves the linear system <tt>AX=B</tt> using the constructed solver.

        @param X initial guess
        @param B right-hand side
        @param InitialGuessIsZero This option has no effect with this smoother
    */
    void Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero = false) const;

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

    // Important note:
    // linearProblem_ must be destroyed before A_, because destructor of linearProblem_ is using A_.
    // In C++, destructor of member objects are called in the reverse order they appear within the declaration for the class.
    // ==18029== Invalid read of size 8
    // ==18029==    at 0xC0780A: Epetra_LinearProblem::GetMatrix() const (Epetra_LinearProblem.h:173)
    // ==18029==    by 0xC5EC27: Amesos_Superlu::PrintTiming() const (Amesos_Superlu.cpp:664)
    // ==18029==    by 0xC628C6: Amesos_Superlu::~Amesos_Superlu() (Amesos_Superlu.cpp:108)

    //! amesos-specific key phrase that denote smoother type
    std::string type_;

    //! parameter list that is used by Amesos internally
    Teuchos::ParameterList paramList_;

    //! Matrix. Not used directly, but held inside of linearProblem_. So we have to keep an RCP pointer to it!
    RCP<Matrix> A_;

    //! Problem that Amesos uses internally.
    RCP<Epetra_LinearProblem> linearProblem_;

    //! pointer to Amesos solver object
    RCP<Amesos_BaseSolver> prec_;

  }; // class AmesosSmoother

  //! Non-member templated function GetAmesosSmoother() returns a new AmesosSmoother object when <Scalar, LocalOrdinal, GlobalOrdinal> == <double, int, int>. Otherwise, an exception is thrown.
  //! This function simplifies the usage of AmesosSmoother objects inside of templates as templates do not have to be specialized for <double, int, int> (see DirectSolver for an example).
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > GetAmesosSmoother(std::string const & type = "", Teuchos::ParameterList const & paramList = Teuchos::ParameterList(), RCP<FactoryBase> AFact = Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "AmesosSmoother cannot be used with Scalar != double, LocalOrdinal != int, GlobalOrdinal != int");
    return Teuchos::null;
  }
  //
  template <>
  inline RCP<MueLu::SmootherPrototype<double, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps> > GetAmesosSmoother<double, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps>(std::string const & type, Teuchos::ParameterList const & paramList, RCP<FactoryBase> AFact) {
    RCP<AmesosSmoother> smoo = rcp( new AmesosSmoother(type, paramList) );
    smoo->SetFactory("A", AFact);
    return smoo;
  }

} // namespace MueLu

#define MUELU_AMESOSSMOOTHER_SHORT

#endif // HAVE_MUELU_AMESOS
#endif // MUELU_AMESOSSMOOTHER_HPP
