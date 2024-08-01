// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SG_PRECONDITIONER_HPP
#define SG_PRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace Kokkos {
namespace Example {

  /*!
   * \brief An abstract class to represent a generic stochastic Galerkin
   * preconditioner
   */
  template<class S, class LO, class GO, class N>
  class SGPreconditioner {
  public:

    //! Constructor
    SGPreconditioner() {}

    //! Destructor
    virtual ~SGPreconditioner() {}

    //! Setup preconditioner
    virtual
    Teuchos::RCP<Tpetra::Operator<S,LO,GO,N> >
    setupPreconditioner(
      const Teuchos::RCP<Tpetra::CrsMatrix<S,LO,GO,N> >& A,
      const Teuchos::RCP<Teuchos::ParameterList>& precParams,
      const Teuchos::RCP<Tpetra::MultiVector<double,LO,GO,N> >& coords) = 0;

  private:

    // Prohibit copying
    SGPreconditioner(const SGPreconditioner&);

    // Prohibit Assignment
    SGPreconditioner& operator=(const SGPreconditioner& b);
  }; // class SGPreconditioner

}
}
#endif // SG_PRECONDITIONER_HPP
