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
#ifndef MUELU_MATLABSMOOTHER_DECL_HPP
#define MUELU_MATLABSMOOTHER_DECL_HPP

#include <Teuchos_ParameterList.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_MATLAB)
#include <Tpetra_CrsMatrix.hpp>
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_MatlabSmoother_fwd.hpp"
#include "MueLu_MatlabUtils_decl.hpp"

namespace MueLu {

  /*!
    @class MatlabSmoother
    @ingroup MueMexClasses 
    @brief Class that encapsulates Matlab smoothers.

    //   This class creates an Matlab preconditioner factory. The factory creates a smoother based on the
    //   type and ParameterList passed into the constructor. See the constructor for more information.
    */

  template <class Scalar = SmootherPrototype<>::scalar_type,
            class LocalOrdinal = typename SmootherPrototype<Scalar>::local_ordinal_type,
            class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
            class Node = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class MatlabSmoother : public SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>
  {
#undef MUELU_MATLABSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors / destructors
    //@{
    //TODO: update doc for Matlab.
    /*! @brief Constructor

      ADD DOCUMENTATION HERE

    */

#ifndef _MSC_VER
    // Avoid error C3772: invalid friend template declaration
    template<class Scalar2, class LocalOrdinal2, class GlobalOrdinal2, class Node2>
    friend class MatlabSmoother;
#endif

    MatlabSmoother(const Teuchos::ParameterList& paramList = Teuchos::ParameterList());

    //! Destructor
    virtual ~MatlabSmoother() { }

    //@}

    void SetParameterList(const Teuchos::ParameterList& paramList);

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    //! @name Computational methods.
    //@{

    /*! @brief Set up the smoother.

    This creates the underlying Matlab smoother object, copies any parameter list options
    supplied to the constructor to the Matlab object, and computes the preconditioner.

    TODO The eigenvalue estimate should come from A_, not the Matlab parameter list.
    */
    void Setup(Level &currentLevel);

    /*! @brief Apply the preconditioner.

    Solves the linear system <tt>AX=B</tt> using the constructed smoother.

    @param X initial guess
    @param B right-hand side
    @param InitialGuessIsZero (optional) If false, some work can be avoided. Whether this actually saves any work depends on the underlying Matlab implementation.
    */
    void Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero = false) const;

    //@}

    //! @name Utilities
    //@{

    RCP<SmootherPrototype> Copy() const;

    //@}

    //! Clone the smoother to a different node type
    template<typename Node2>
    RCP<MueLu::MatlabSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node2> >
    clone(const RCP<Node2>& node2, const Teuchos::RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node2> >& A_newnode) const;

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

    size_t getNodeSmootherComplexity() const {return Teuchos::OrdinalTraits<size_t>::invalid();}


    //@}

  private:

    //! List of arguments to the MATLAB setup function besides "A", in order
    mutable std::string needsSetup_;

    //! Amount of solve data (besides A, LHS & RHS)
    size_t solveDataSize_;

    //! List of data generated by setup which will be sent to solve after "A", "LHS" and "RHS"
    std::vector<Teuchos::RCP<MuemexArg> > solveData_;
    
    //! Matlab setup function
    std::string setupFunction_;

    //! Matlab solve function
    std::string solveFunction_;

    //! Matrix, (maybe) used in apply 
    mutable RCP<Matrix> A_;
    
  }; // class MatlabSmoother

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  template<typename Node2>
  RCP<MueLu::MatlabSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node2> >
  MatlabSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node>::clone(const RCP<Node2>& node2, const RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node2> >& A_newnode) const {
    const ParameterList& paramList = this->GetParameterList();

    RCP<MatlabSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node2> > cloneSmoother =
      rcp(new MatlabSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node2>(paramList));   

    cloneSmoother->needsSetup_    = needsSetup_;
    cloneSmoother->setupFunction_ = setupFunction_;
    cloneSmoother->solveFunction_  = solveFunction_;
    cloneSmoother->A_             = A_;    
    
    for(size_t i=0; i< solveData_.size(); i++)
      cloneSmoother->solveData_->push_back(solveData_[i]);
    cloneSmoother->SetParameterList(paramList);
    cloneSmoother->IsSetup(this->IsSetup());
    return cloneSmoother;
  }


} // namespace MueLu

#define MUELU_MATLABSMOOTHER_SHORT
#endif // HAVE_MUELU_TPETRA && HAVE_MUELU_MATLAB
#endif // MUELU_MATLABSMOOTHER_DECL_HPP
