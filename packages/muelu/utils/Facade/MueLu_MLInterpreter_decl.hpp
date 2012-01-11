/*
 * MueLu_Interpreter_decl.hpp
 *
 *  Created on: Dec 7, 2011
 *      Author: wiesner
 */

#ifndef MUELU_MLINTERPRETER_DECL_HPP_
#define MUELU_MLINTERPRETER_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_MLInterpreter_fwd.hpp"

#include "MueLu_Hierarchy_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"
#include <Xpetra_Operator.hpp>
#include <Teuchos_ParameterList.hpp>

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

    static RCP<Hierarchy> Setup(const Teuchos::ParameterList & params, const RCP<Operator> & A, const RCP<MultiVector> nsp = Teuchos::null);

    static void FillMLParameterList(Teuchos::ParameterList & params);

    static RCP<SmootherFactory> GetCoarsestSolverFactory(const Teuchos::ParameterList & params);

    static RCP<SmootherFactory> GetSmootherFactory(const Teuchos::ParameterList & params, int level);


  private:


  }; // class MLInterpreter

} // namespace MueLu

#define MUELU_MLINTERPRETER_SHORT
#endif /* MUELU_MLINTERPRETER_DECL_HPP_ */
