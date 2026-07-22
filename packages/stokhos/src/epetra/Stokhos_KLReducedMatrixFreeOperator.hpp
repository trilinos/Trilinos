// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_KL_REDUCED_MATRIX_FREE_OPERATOR_HPP
#define STOKHOS_KL_REDUCED_MATRIX_FREE_OPERATOR_HPP

#include "Stokhos_SGOperator.hpp"
#include "EpetraExt_MultiComm.h"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_EpetraSparse3Tensor.hpp"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "Epetra_MultiVector.h"
#include "Stokhos_VectorOrthogPoly.hpp"
#include "Stokhos_VectorOrthogPolyTraitsEpetra.hpp"
#include "Stokhos_EpetraVectorOrthogPoly.hpp"
#include "Stokhos_MatrixFreeOperator.hpp"

namespace Stokhos {
    
  /*! 
   * \brief An Epetra operator representing the block stochastic Galerkin
   * operator.
   */
  class KLReducedMatrixFreeOperator : public Stokhos::SGOperator {
      
  public:

    //! Constructor 
    KLReducedMatrixFreeOperator(
      const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm,
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis,
      const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk,
      const Teuchos::RCP<const Epetra_Map>& domain_base_map,
      const Teuchos::RCP<const Epetra_Map>& range_base_map,
      const Teuchos::RCP<const Epetra_Map>& domain_sg_map,
      const Teuchos::RCP<const Epetra_Map>& range_sg_map,
      const Teuchos::RCP<Teuchos::ParameterList>& params);
    
    //! Destructor
    virtual ~KLReducedMatrixFreeOperator();

    /** \name Stokhos::SGOperator methods */
    //@{

    //! Setup operator
    virtual void setupOperator(
      const Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly >& poly);

    //! Get SG polynomial
    virtual Teuchos::RCP< Stokhos::EpetraOperatorOrthogPoly > 
    getSGPolynomial();

    //! Get SG polynomial
    virtual Teuchos::RCP<const Stokhos::EpetraOperatorOrthogPoly > 
    getSGPolynomial() const;

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

    //! Setup KL blocks
    void setup();

  private:
    
    //! Private to prohibit copying
    KLReducedMatrixFreeOperator(const KLReducedMatrixFreeOperator&);
    
    //! Private to prohibit copying
    KLReducedMatrixFreeOperator& operator=(const KLReducedMatrixFreeOperator&);
    
  protected:
    
    //! Label for operator
    std::string label;

    //! Stores SG parallel communicator
    Teuchos::RCP<const EpetraExt::MultiComm> sg_comm;

    //! Stochastic Galerking basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > sg_basis;

    //! Stores Epetra Cijk tensor
    Teuchos::RCP<const Stokhos::EpetraSparse3Tensor> epetraCijk;
    
    //! Stores domain base map
    Teuchos::RCP<const Epetra_Map> domain_base_map;

    //! Stores range base map
    Teuchos::RCP<const Epetra_Map> range_base_map;

    //! Stores domain SG map
    Teuchos::RCP<const Epetra_Map> domain_sg_map;

    //! Stores range SG map
    Teuchos::RCP<const Epetra_Map> range_sg_map;

    //! Short-hand for Cijk
    typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;

    //! Stores triple product tensor
    Teuchos::RCP<const Cijk_type> Cijk;

    //! Stores operators
    Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly > block_ops;

    //! Algorithmic parameters
    Teuchos::RCP<Teuchos::ParameterList> params;

    //! Flag indicating whether transpose was selected
    bool useTranspose;

    //! Number of terms in expansion
    int expansion_size;

    //! Number of blocks
    int num_blocks;

    //! Number of KL terms
    int num_KL;

    //! Number of computed KL terms
    int num_KL_computed;

    //! Mean block
    Teuchos::RCP<Epetra_CrsMatrix> mean;

    //! Block map for vectorized-matrices
    Teuchos::RCP<Epetra_Map> block_vec_map;

    //! Polynomial sorting vectorized matrix coefficients
    Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly> block_vec_poly;

    //! Dot products of KL eigenvectors and Jacobian blocks
    Teuchos::Array< Teuchos::Array<double> > dot_products;

    //! Sparse KL coefficients
    Teuchos::RCP< Stokhos::Sparse3Tensor<int,double> > sparse_kl_coeffs;

    //! KL blocks
    Teuchos::Array< Teuchos::RCP<Epetra_CrsMatrix> > kl_blocks;

    //! KL blocks as operators
    Teuchos::RCP< Stokhos::EpetraOperatorOrthogPoly > kl_ops;

    //! Matrix-Free operator using KL operators
    Teuchos::RCP< Stokhos::MatrixFreeOperator > kl_mat_free_op;

    //! Tolerance for dropping entries in sparse 3 tensor
    double drop_tolerance;

    //! Whether to do KL error tests (can be expensive)
    bool do_error_tests;

  }; // class KLReducedMatrixFreeOperator
  
} // namespace Stokhos

#endif // STOKHOS_KL_REDUCED_MATRIX_FREE_OPERATOR_HPP
