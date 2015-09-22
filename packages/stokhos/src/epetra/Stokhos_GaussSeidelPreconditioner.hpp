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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_GAUSS_SEIDEL_PRECONDITIONER_HPP
#define STOKHOS_GAUSS_SEIDEL_PRECONDITIONER_HPP

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_NOX

#include "Teuchos_RCP.hpp"

#include "Stokhos_SGPreconditioner.hpp"
#include "EpetraExt_MultiComm.h"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_EpetraSparse3Tensor.hpp"
#include "Epetra_Map.h"
#include "Epetra_Export.h"
#include "NOX_Epetra_LinearSystem.H"
#include "Teuchos_ParameterList.hpp"
#include "EpetraExt_BlockMultiVector.h"

namespace Stokhos {
    
  /*! 
   * \brief An Epetra operator representing applying the mean in a block
   * stochastic Galerkin expansion.
   */
  class GaussSeidelPreconditioner : public Stokhos::SGPreconditioner {
      
  public:

    //! Constructor 
    GaussSeidelPreconditioner(
      const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm,
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis,
      const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk,
      const Teuchos::RCP<const Epetra_Map>& base_map,
      const Teuchos::RCP<const Epetra_Map>& sg_map,
      const Teuchos::RCP<NOX::Epetra::LinearSystem>& det_solver,
      const Teuchos::RCP<Teuchos::ParameterList>& params);
    
    //! Destructor
    virtual ~GaussSeidelPreconditioner();

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

  private:
    
    //! Private to prohibit copying
    GaussSeidelPreconditioner(const GaussSeidelPreconditioner&);
    
    //! Private to prohibit copying
    GaussSeidelPreconditioner& operator=(const GaussSeidelPreconditioner&);
    
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

    //! Whether we have parallelism over stochastic blocks
    bool is_stoch_parallel;

    //! Stores stochastic part of row map
    Teuchos::RCP<const Epetra_BlockMap> stoch_row_map;

    //! Deterministic solver
    Teuchos::RCP<NOX::Epetra::LinearSystem> det_solver;

    //! Preconditioner and solver parameters
    Teuchos::RCP<Teuchos::ParameterList> params;

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

    //! Stores block vector of right-hand-sides
      mutable Teuchos::RCP<EpetraExt::BlockMultiVector> sg_df_block;

      //! Stores block residual vector to compute residual norm
      mutable Teuchos::RCP<EpetraExt::BlockMultiVector> sg_y_block;
      
      //! Stores K_0*x for the most recently computed x
      mutable Teuchos::RCP<Epetra_MultiVector> kx;

      //! Flag indicating whether stochastic blocks are distributed
      bool is_parallel;

      //! Stores global column map
      Teuchos::RCP<const Epetra_BlockMap> sg_col_map;

      //! Stores exporter from column map to row map
      Teuchos::RCP<Epetra_Export> col_exporter;

      //! Stores off-processor contributions to right-hand-sides
      mutable Teuchos::RCP<EpetraExt::BlockMultiVector> sg_df_col;

      //! Stores summed off-processor contributions
      mutable Teuchos::RCP<EpetraExt::BlockMultiVector> sg_df_tmp;


  }; // class GaussSeidelPreconditioner
  
} // namespace Stokhos

#endif

#endif // STOKHOS_GAUSS_SEIDEL_PRECONDITIONER_HPP
