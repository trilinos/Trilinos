// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_KLMATRIX_FREE_OPERATOR_HPP
#define STOKHOS_KLMATRIX_FREE_OPERATOR_HPP

#include "Stokhos_SGOperator.hpp"
#include "EpetraExt_MultiComm.h"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_EpetraSparse3Tensor.hpp"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

namespace Stokhos {
    
  /*! 
   * \brief An Epetra operator representing the block stochastic Galerkin
   * operator.
   */
  class KLMatrixFreeOperator : public Stokhos::SGOperator {
      
  public:

    //! Constructor 
    KLMatrixFreeOperator(
      const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm,
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis,
      const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk,
      const Teuchos::RCP<const Epetra_Map>& domain_base_map,
      const Teuchos::RCP<const Epetra_Map>& range_base_map,
      const Teuchos::RCP<const Epetra_Map>& domain_sg_map,
      const Teuchos::RCP<const Epetra_Map>& range_sg_map,
      const Teuchos::RCP<Teuchos::ParameterList>& params);
    
    //! Destructor
    virtual ~KLMatrixFreeOperator();

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
    
    //! Returns a character std::string describing the operator
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

  private:
    
    //! Private to prohibit copying
    KLMatrixFreeOperator(const KLMatrixFreeOperator&);
    
    //! Private to prohibit copying
    KLMatrixFreeOperator& operator=(const KLMatrixFreeOperator&);
    
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

    //! Whether we have parallelism over stochastic blocks
    bool is_stoch_parallel;

    //! Stores operator column SG map
    Teuchos::RCP<Epetra_Map> global_col_map;

    //! Stores operator column SG map for transpose
    Teuchos::RCP<Epetra_Map> global_col_map_trans;

    //! Stores stochastic part of column map
    Teuchos::RCP<const Epetra_BlockMap> stoch_col_map;

    //! Importer from domain map to column map
    Teuchos::RCP<Epetra_Import> col_importer;

    //! Importer from range map to column map
    Teuchos::RCP<Epetra_Import> col_importer_trans;

    //! Short-hand for Cijk
    typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;

    //! Stores triple product tensor
    Teuchos::RCP<const Cijk_type> Cijk;

    //! Stores operators
    Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly > block_ops;

    //! Flag indicating whether operator be scaled with <\psi_i^2>
    bool scale_op;

    //! Flag indicating whether to include mean term
    bool include_mean;

    //! Flag indicating whether transpose was selected
    bool useTranspose;

     //! Number of terms in expansion
    int expansion_size;

    //! Number of Jacobian blocks (not necessarily equal to expansion_size)
    int num_blocks;

    //! Maximum number of matvecs in Apply
    int max_num_mat_vec;

    //! Temporary to store result of importing input into column map
    mutable Teuchos::RCP<Epetra_MultiVector> input_col;

    //! Temporary to store result of importing input into column map (transpose)
    mutable Teuchos::RCP<Epetra_MultiVector> input_col_trans;

    //! MultiVectors for each block for Apply() input
    mutable Teuchos::Array< Teuchos::RCP<const Epetra_MultiVector> > input_block;

    //! MultiVectors for each block for Apply() result
    mutable Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> > result_block;   

    //! Temporary multivector used in Apply()
    mutable Teuchos::RCP<Epetra_MultiVector> tmp;

    //! Temporary multivector used in Apply() for transpose
    mutable Teuchos::RCP<Epetra_MultiVector> tmp_trans;

    //! Starting k iterator
    Cijk_type::k_iterator k_begin;

    //! Ending k iterator
    Cijk_type::k_iterator k_end;

  }; // class KLMatrixFreeOperator
  
} // namespace Stokhos

#endif // STOKHOS_KLMATRIX_FREE_OPERATOR_HPP
