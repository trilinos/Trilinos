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

#ifndef STOKHOS_APPROX_SCHUR_COMPLEMENT_PRECONDITIONER_HPP
#define STOKHOS_APPROX_SCHUR_COMPLEMENT_PRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"

#include "Stokhos_SGPreconditioner.hpp"
#include "EpetraExt_MultiComm.h"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_EpetraSparse3Tensor.hpp"
#include "Epetra_Map.h"
#include "Stokhos_AbstractPreconditionerFactory.hpp"
#include "Teuchos_ParameterList.hpp"
#include "EpetraExt_BlockMultiVector.h"

namespace Stokhos {
    
  /*! 
   * \brief A stochastic preconditioner based on applying one iteration of
   * approximate Gauss-Seidel.
   */
  class ApproxSchurComplementPreconditioner : public Stokhos::SGPreconditioner {
      
  public:

    //! Constructor 
    ApproxSchurComplementPreconditioner(
      const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm,
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis,
      const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk,
      const Teuchos::RCP<const Epetra_Map>& base_map,
      const Teuchos::RCP<const Epetra_Map>& sg_map,
      const Teuchos::RCP<Stokhos::AbstractPreconditionerFactory>& prec_factory,
      const Teuchos::RCP<Teuchos::ParameterList>& params);
    
    //! Destructor
    virtual ~ApproxSchurComplementPreconditioner();

    /** \name Stokhos::SGPreconditioner methods */
    //@{

    //! Setup preconditioner
    virtual void 
    setupPreconditioner(const Teuchos::RCP<Stokhos::SGOperator>& sg_op, 
			const Epetra_Vector& x);

    //@}

    /** \name Epetra_Operator methods */
    //@{

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

    //@}

  protected:

    void compute_index(int l, int& n_begin, int& n_end) const;

    void multiply_block(int row_begin, int row_end, int col_begin, int col_end,
			double alpha,
			const EpetraExt::BlockMultiVector& Input, 
			EpetraExt::BlockMultiVector& Result) const;

    void divide_diagonal_block(int row_begin, int row_end,
			       const EpetraExt::BlockMultiVector& Input, 
			       EpetraExt::BlockMultiVector& Result) const;

  private:
    
    //! Private to prohibit copying
    ApproxSchurComplementPreconditioner(const ApproxSchurComplementPreconditioner&);
    
    //! Private to prohibit copying
    ApproxSchurComplementPreconditioner& operator=(const ApproxSchurComplementPreconditioner&);
    
  protected:
    
    //! Label for operator
    std::string label;

    //! Stores SG parallel communicator
    Teuchos::RCP<const EpetraExt::MultiComm> sg_comm;

    //! Stochastic Galerking basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > sg_basis;

    //! Stores Epetra Cijk tensor
    Teuchos::RCP<const Stokhos::EpetraSparse3Tensor> epetraCijk;
    
    //! Stores base map
    Teuchos::RCP<const Epetra_Map> base_map;

    //! Stores SG map
    Teuchos::RCP<const Epetra_Map> sg_map;

    //! Stores factory for building mean preconditioner
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> prec_factory;

    //! Stores mean preconditioner
    Teuchos::RCP<Epetra_Operator> mean_prec;

    //! Flag indicating whether transpose was selected
    bool useTranspose;

    //! Pointer to the SG operator.
    Teuchos::RCP<Stokhos::SGOperator> sg_op;

    //! Pointer to the PCE expansion of Jacobian.
    Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly > sg_poly;
    
    //! Short-hand for Cijk
    typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;

    //! Pointer to triple product
    Teuchos::RCP<const Cijk_type > Cijk;

    //! Flag indicating whether operator be scaled with <\psi_i^2>
    bool scale_op;

    //! Use symmetric Gauss-Seidel
    bool symmetric;

    //! Limit Gauss-Seidel loop to linear terms
    bool only_use_linear;

    //! Starting k iterator
    Cijk_type::k_iterator k_begin;

    //! Ending k iterator
    Cijk_type::k_iterator k_end;

    //! Maximum number of matvecs in Apply
    int max_num_mat_vec;

    //! Temporary vector for storing matrix-vector products
    mutable Teuchos::RCP<Epetra_MultiVector> tmp;

    //! Temporary vector for storing rhs in Gauss-Seidel loop
    mutable Teuchos::RCP<EpetraExt::BlockMultiVector> rhs_block;

  }; // class ApproxSchurComplementPreconditioner
  
} // namespace Stokhos

#endif // STOKHOS_APPROX_SCHUR_COMPLEMENT_PRECONDITIONER_HPP
