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

#ifndef STOKHOS_DIAG_EPETRA_OP_HPP
#define STOKHOS_DIAG_EPETRA_OP_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Time.hpp"

#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Stokhos_VectorOrthogPoly.hpp"
#include "Stokhos_VectorOrthogPolyTraitsEpetra.hpp"

namespace Stokhos {
    
  /*! 
   * \brief An Epetra operator representing the block stochastic Galerkin
   * operator.
   */
  class DiagEpetraOp : public Epetra_Operator {
      
  public:

    //! Constructor 
    DiagEpetraOp(
      const Teuchos::RCP<const Epetra_Map>& domain_base_map_,
      const Teuchos::RCP<const Epetra_Map>& range_base_map_,
      const Teuchos::RCP<const Epetra_Map>& domain_sg_map_,
      const Teuchos::RCP<const Epetra_Map>& range_sg_map_,
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk,
      const Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly >& ops);
    
    //! Destructor
    virtual ~DiagEpetraOp();

    //! Reset operator blocks
    virtual void 
    reset(const Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly >& ops);

    //! Get operator blocks
    virtual Teuchos::RCP<const Stokhos::EpetraOperatorOrthogPoly >
    getOperatorBlocks() const;

    //! Get operator blocks
    virtual Teuchos::RCP< Stokhos::EpetraOperatorOrthogPoly >
    getOperatorBlocks();
    
    //! Set to true if the transpose of the operator is requested
    virtual int SetUseTranspose(bool UseTranspose);
    
    /*! 
     * \brief Returns Diagonal blocks of SG matrix when
     * PC coefficients of the SG matrix are given
     */
    virtual int Apply(std::vector< Teuchos::RCP< const Epetra_CrsMatrix> >& sg_J_all, std::vector< Teuchos::RCP< Epetra_CrsMatrix> >& sg_Kkk_all) const;

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

    /*!
     * \brief Returns the time spent applying this operator
     */
    virtual const double ApplyTime() const{
      return this->ApplyTimer->totalElapsedTime(false);};

  private:
    
    //! Private to prohibit copying
    DiagEpetraOp(const DiagEpetraOp&);
    
    //! Private to prohibit copying
    DiagEpetraOp& operator=(const DiagEpetraOp&);
    
  protected:
    
    //! Label for operator
    std::string label;
    
    //! Stores domain base map
    Teuchos::RCP<const Epetra_Map> domain_base_map;

    //! Stores range base map
    Teuchos::RCP<const Epetra_Map> range_base_map;

    //! Stores domain SG map
    Teuchos::RCP<const Epetra_Map> domain_sg_map;

    //! Stores range SG map
    Teuchos::RCP<const Epetra_Map> range_sg_map;

    //! Stochastic Galerking basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > sg_basis;

    //! Stores triple product tensor
    Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk;

    //! Stores operators
    Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly > block_ops;

    //! Flag indicating whether transpose was selected
    bool useTranspose;

    //! Number of terms in expansion
    int expansion_size;

    //! Number of Jacobian blocks (not necessarily equal to expansion_size)
    int num_blocks;

    //! MultiVectors for each block for Apply() input
    mutable Teuchos::Array< Teuchos::RCP<const Epetra_MultiVector> > input_block;

    //! MultiVectors for each block for Apply() result
    mutable Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> > result_block;

    //! Temporary multivector used in Apply()
    mutable Teuchos::RCP<Epetra_MultiVector> tmp;

    //! Temporary multivector used in Apply() for transpose
    mutable Teuchos::RCP<Epetra_MultiVector> tmp_trans;

    //! Operation Timer
    Teuchos::RCP<Teuchos::Time> ApplyTimer;

  }; // class DiagEpetraOp
  
} // namespace Stokhos

#endif // STOKHOS_DIAG_EPETRA_OP_HPP
