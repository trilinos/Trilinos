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
#ifndef MUELU_EPETRAOPERATOR_HPP
#define MUELU_EPETRAOPERATOR_HPP

// Turns a MueLu::Hierarchy into a Epetra_Operator. 
// It allows to use MueLu as a preconditionner for AztecOO (for instance).

// TODO: Is MueLu::EpetraOperator a good name for this adapter? (
//       For the MueLu/Belos adapter, the name of the class is Belos::MueLuOp)

#include <Epetra_Operator.h>
#include "MueLu_Hierarchy.hpp"
//TODO: Kokkos headers

namespace MueLu {

  class EpetraOperator
    : public Epetra_Operator
  {
    
    typedef Kokkos::DefaultNode::DefaultNodeType Node;
    typedef Kokkos::DefaultKernels<double,int,Node>::SparseOps LocalMatOps;

    typedef Xpetra::Matrix<double, int, int, Node, LocalMatOps> Matrix;
    typedef MueLu::Utils<double, int, int, Node, LocalMatOps>     Utils;

  public:
    
    //! @name Constructor/Destructor
    //@{ 
    
    //! Constructor
    EpetraOperator(const RCP<MueLu::Hierarchy<double, int, int, Node, LocalMatOps> > & H) : Hierarchy_(H) { }
    
    //! Destructor.
    virtual ~EpetraOperator() { }

    //@}

    int SetUseTranspose(bool UseTransposeBool);
  
    //! @name Mathematical functions
    //@{ 

    //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*! 
      \param In
      X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
      \param Out
      Y -A Epetra_MultiVector of dimension NumVectors containing result.

      \return Integer error code, set to 0 if successful.
    */
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    /*! 
      \param In
      X - A Epetra_MultiVector of dimension NumVectors to solve for.
      \param Out
      Y -A Epetra_MultiVector of dimension NumVectors containing result.

      \return Integer error code, set to 0 if successful.

      \warning In order to work with AztecOO, any implementation of this method must 
      support the case where X and Y are the same object.
    */
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Returns the infinity norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

       \warning This method must not be called unless HasNormInf() returns true.
    */ 
    double NormInf() const;
    //@}
  
    //! @name Attribute access functions
    //@{ 

    //! Returns a character string describing the operator
    const char * Label() const;

    //! Returns the current UseTranspose setting.
    bool UseTranspose() const;

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    bool HasNormInf() const;

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    const Epetra_Comm & Comm() const;

    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map & OperatorDomainMap() const;

    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map & OperatorRangeMap() const;

    //! Direct access to the underlying MueLu::Hierarchy.
    RCP<MueLu::Hierarchy<double, int, int, Node, LocalMatOps> > GetHierarchy() const;    

    //@}

  private:
    
    RCP<MueLu::Hierarchy<double, int, int, Node, LocalMatOps> > Hierarchy_;

  };

} // namespace

#endif // MUELU_EPETRAOPERATOR_HPP
