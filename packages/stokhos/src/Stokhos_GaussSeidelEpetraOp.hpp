// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_GAUSS_SEIDEL_EPETRA_OP_HPP
#define STOKHOS_GAUSS_SEIDEL_EPETRA_OP_HPP

#include "Teuchos_RCP.hpp"

#include "NOX_Epetra_LinearSystem.H"    // base class
#include "Stokhos.hpp"

#include "Epetra_Operator.h"
#include "Epetra_Map.h"

namespace Stokhos {
    
  /*! 
   * \brief An Epetra operator representing applying the mean in a block
   * stochastic Galerkin expansion.
   */
  class GaussSeidelEpetraOp : public Epetra_Operator {
      
  public:

    //! Constructor 
  GaussSeidelEpetraOp(
  const Teuchos::RCP<const Epetra_Map>& base_map_,
  const Teuchos::RCP<const Epetra_Map>& sg_map_,
  unsigned int num_blocks_,
  Teuchos::ParameterList& linearSolverParams,
  const Teuchos::RCP<NOX::Epetra::LinearSystem>& detsolve_,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_,
  const Teuchos::RCP<Epetra_Operator>& J);
    
    //! Destructor
    virtual ~GaussSeidelEpetraOp();

    void setOperatorAndConstructPreconditioner(
      const Teuchos::RCP<Epetra_Operator>& J, const Epetra_Vector& x);

    //! Set to true if the transpose of the operator is requested
    virtual int SetUseTranspose(bool UseTranspose);
    
    /*! 
     * \brief Returns the result of a Epetra_Operator applied to a 
     * Epetra_MultiVector Input in Result as described above.
     */
    virtual int Apply(const Epetra_MultiVector& Input, 
                      Epetra_MultiVector& Result) const;

    /*! 
     * \brief Returns the result of the inverse of the operator applied to a 
     * Epetra_MultiVector Input in Result as described above.
     */
    virtual int ApplyInverse(const Epetra_MultiVector& X, 
                             Epetra_MultiVector& Y) const;
    
    //! Returns an approximate infinity norm of the operator matrix.
    virtual double NormInf() const;
    
    //! Returns a character string describing the operator
    virtual const char* Label () const;
  
    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const;
    
    /*! 
     * \brief Returns true if the \e this object can provide an 
     * approximate Inf-norm, false otherwise.
     */
    virtual bool HasNormInf() const;

    /*! 
     * \brief Returns a reference to the Epetra_Comm communicator 
     * associated with this operator.
     */
    virtual const Epetra_Comm & Comm() const;

    /*!
     * \brief Returns the Epetra_Map object associated with the 
     * domain of this matrix operator.
     */
    virtual const Epetra_Map& OperatorDomainMap () const;

    /*! 
     * \brief Returns the Epetra_Map object associated with the 
     * range of this matrix operator.
     */
    virtual const Epetra_Map& OperatorRangeMap () const;

  private:
    
    //! Private to prohibit copying
    GaussSeidelEpetraOp(const GaussSeidelEpetraOp&);
    
    //! Private to prohibit copying
    GaussSeidelEpetraOp& operator=(const GaussSeidelEpetraOp&);
    
  protected:
    
    //! Label for operator
    std::string label;
    
    //! Stores base map
    Teuchos::RCP<const Epetra_Map> base_map;

    //! Stores SG map
    Teuchos::RCP<const Epetra_Map> sg_map;

    //! Flag indicating whether transpose was selected
    bool useTranspose;

    //! Number of blocks
    unsigned int num_blocks;
  
     //! Pointer to deterministic solver
    Teuchos::RCP<NOX::Epetra::LinearSystem> detsolve;

    //! Pointer to Cijk
    Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk;

    //! Pointer to the Jacobian operator.
    mutable Teuchos::RCP<Epetra_Operator> jacPtr;

    //! Pointer to the Stokhos matrixfree epetra operator.
    mutable Teuchos::RCP<Stokhos::MatrixFreeEpetraOp> stokhos_op;

    //! Pointer to the PCE expansion of Jacobian.
    mutable Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_Operator> > sg_J_poly;

    //! Pointer to the deterministic vector 
    Teuchos::RCP<Epetra_Vector> detvec;

   //! Parameter list
   Teuchos::ParameterList& params;

  }; // class GaussSeidelEpetraOp
  
} // namespace Stokhos

#endif // STOKHOS_GAUSS_SEIDEL_EPETRA_OP_HPP
