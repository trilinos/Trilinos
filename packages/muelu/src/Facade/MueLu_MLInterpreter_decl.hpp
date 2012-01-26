/*
 * MueLu_Interpreter_decl.hpp
 *
 *  Created on: Dec 7, 2011
 *      Author: wiesner
 */

#ifndef MUELU_MLINTERPRETER_DECL_HPP
#define MUELU_MLINTERPRETER_DECL_HPP

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_MLInterpreter_fwd.hpp"

#include "MueLu_Hierarchy_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"

namespace MueLu {

  /*!
    @class MLInterpreter class.
    @brief Class that accepts ML-style parameters and builds a MueLu preconditioner.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class MLInterpreter : public BaseClass {
#undef MUELU_MLINTERPRETER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    MLInterpreter();

    //! Destructor.
    virtual ~MLInterpreter();

    //@}

    //! @name Build methods.
    //@{

    //! Build a MueLu::Hierarchy from an ML parameter list
    static RCP<Hierarchy> Setup(const Teuchos::ParameterList & params, const RCP<Operator> & A, const RCP<MultiVector> nsp = Teuchos::null); // TODO: should it be renamed Build() ?

    //@}

    //! @name Helper functions translating parameter list to factories
    //@{

    //! Read coarse solver options and build the corresponding smoother factory
    static RCP<SmootherFactory> GetCoarsestSolverFactory(const Teuchos::ParameterList & params);

    //! Read smoother options and build the corresponding smoother factory
    static RCP<SmootherFactory> GetSmootherFactory(const Teuchos::ParameterList & params, int level);

    //@}

    //! Build an example of valid ML parameter list
    static void FillMLParameterList(Teuchos::ParameterList & params);

  private:


  }; // class MLInterpreter

} // namespace MueLu

#define MUELU_MLINTERPRETER_SHORT
#endif /* MUELU_MLINTERPRETER_DECL_HPP */
