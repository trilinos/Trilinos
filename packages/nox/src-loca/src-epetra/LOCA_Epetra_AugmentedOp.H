//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
                                                                                
#ifndef LOCA_EPETRA_AUGMENTEDOP_H
#define LOCA_EPETRA_AUGMENTEDOP_H

#include "Teuchos_RCP.hpp"

#include "Epetra_Operator.h"	       // base class
#include "Epetra_LocalMap.h"           // class data element
#include "Epetra_MultiVector.h"        // class data element

#include "NOX_Common.H"                // for std::string data member
#include "NOX_Abstract_MultiVector.H"  // for DenseMatrix

#include "Teuchos_LAPACK.hpp" // class data element

// Forward declarations
namespace LOCA {
  class GlobalData;
}

namespace LOCA {

  namespace Epetra {

    /*!
     * \brief %Epetra operator representing a \f$n+m\f$ bordered matrix.
     */
    /*!
     * The %LOCA::Epetra::AugmentedOp is an Epetra_Operator representing the
     * \f$n+m\f$ bordered matrix
     * \f[
     *     \begin{bmatrix}
     *        J & A \\
     *        B^T & C
     *     \end{bmatrix}
     * \f]
     * where \f$J\f$ is an Epetra_Operator representing an \f$n\times n\f$
     * matrix, and \f$A\f$ and \f$B\f$ are length \f$n\f$ Epetra_MultiVector's
     * with \f$m\f$ columns, and \f$C\f$ is an \f$m\times m\f$ dense matrix.
     * It is assumed the Epetra_Map's for \f$A\f$, \f$B\f$, and
     * \f$J\f$ are the same and the corresponding map for the bordered
     * matrix is constructed from this map by storing the additional components
     * on processor 0.  The buildEpetraAugmentedMultiVec() method can be used 
     * to construct an Epetra_MultiVector using this map, a supplied 
     * length \f$n\f$ Epetra_MultiVector and an \f$m\times m\f$ matrix, while 
     * setEpetraAugmentedMultiVec() splits an extended multivector into its 
     * length \f$n\f$ and \f$m\f$ components.  The 
     * Apply() method performs the \f$n+m\times n+m\f$ matrix multiplication
     * while ApplyInverse() uses a block-elimination algorithm to compute
     * the inverse using the ApplyInverse() method of the underlying operator
     * \f$J\f$.  In this way, linear systems of the form
     * \f[
     *     \begin{bmatrix}
     *        J & A \\
     *        B^T & C
     *     \end{bmatrix}
     *     \begin{bmatrix}
     *       X \\
     *       Y
     *     \end{bmatrix} =
     *     \begin{bmatrix}
     *       F \\
     *       G
     *     \end{bmatrix}
     * \f]
     * can be solved in a matrix-free mode using the Apply() method.  This 
     * operator can also represent a preconditioner of the form
     * \f[
     *     \begin{bmatrix}
     *        M & A \\
     *        B^T & C
     *     \end{bmatrix}
     * \f]
     * using the ApplyInvese() method, where \f$M\f$ is a preconditioner 
     * for \f$J\f$.  Note that if \f$J\f$ is nearly singular, the 
     * preconditioner should not be too good because otherwise the
     * preconditining operation represented by ApplyInverse() becomes unstable.
     */
    class AugmentedOp : public Epetra_Operator {
      
    public:

      //! Constructor 
      /*!
       * Builds the bordered operator using the supplied operator \em jac 
       * and Epetra_Vector's \em a and \em b.  It is assumed \em a, \em b,
       * and \em jac all have the same map.
       */
      AugmentedOp(const Teuchos::RCP<LOCA::GlobalData>& global_data,
		  const Teuchos::RCP<Epetra_Operator>& jac, 
		  const Teuchos::RCP<const Epetra_MultiVector>& a, 
		  const Teuchos::RCP<const Epetra_MultiVector>& b,
                  const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> c);

      //! Destructor
      virtual ~AugmentedOp();

      //! If set true, transpose of this operator will be applied.
      /*! 
       * Note that is only valid if the underlying operator \f$J\f$
       * supports a transpose.
       */
      virtual int SetUseTranspose(bool UseTranspose);

      /*! 
       * \brief Returns the result of a Epetra_Operator applied to a 
       * Epetra_MultiVector Input in Result.
       */
      /*!
       * Computes the extended matrix-vector product
       * \f[
       *     \begin{bmatrix}
       *        J & A \\
       *        B^T & C
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *       X \\
       *       Y
       *     \end{bmatrix} =
       *     \begin{bmatrix}
       *       JX + AY \\
       *       B^T X + CY
       *     \end{bmatrix}
       * \f]
       * or its transpose if UseTranpose() is \em true.
       */
      virtual int Apply(const Epetra_MultiVector& Input, 
			Epetra_MultiVector& Result) const;

      /*!
       * \brief Returns the result of a Epetra_Operator inverse applied to 
       * an Epetra_MultiVector Input in Result.
       */
      /*!
       * Solves the extended system
       * \f[
       *     \begin{bmatrix}
       *        J & A \\
       *        B^T & C
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *       X \\
       *       Y
       *     \end{bmatrix} =
       *     \begin{bmatrix}
       *       F \\
       *       G
       *     \end{bmatrix}
       * \f]
       * using the following block-elimination algorithm:  
       * \f[
       *  \begin{split}
       *    X_1 &= J^{-1} F, \\
       *    X_2 &= J^{-1} A, \\
       *    Y &= (C-B^T X_2)^{-1}(G-B^T X_1), \\ 
       *    X &= X_1 - X_2 Y
       *  \end{split}
       * \f]
       * If UseTranpose() is \em true, the tranpose of the system is solved.
       */
      virtual int ApplyInverse(const Epetra_MultiVector& X, 
			       Epetra_MultiVector& Y) const;
  
      //! Returns the infinity norm of the bordered matrix.
      /*!
       * This is defined only if %NormInf() of the underlying operator \f$J\f$
       * is defined and is given by 
       * \f$\|J\|_\infty+\|A\|_\infty+\|B\|_\infty\f$.
       */
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

      //! Initialiazes operator for a solve
      virtual void init(const Epetra_MultiVector& x);

      /*!
       * \brief Builds an extended  vector from components
       */
      /*!
       * Builds an extended vector using the map representing the bordered
       * matrix.  If \em doCopy is \em true, the contents of \em x are
       * copied into the extended vector, otherwise only space for the 
       * extended vector is created.
       */
      virtual Teuchos::RCP<Epetra_MultiVector> 
      buildEpetraAugmentedMultiVec(
			     const Epetra_MultiVector& x,
			     const NOX::Abstract::MultiVector::DenseMatrix *y, 
			     bool doCopy) const;

      /*!
       * \brief Sets components from extended vector
       */
      /*!
       * Splits the extended vector \em augMultiVec into components \em x and
       * \em y by copying values out of \em extVec.
       */
      virtual void 
      setEpetraAugmentedMultiVec(Epetra_MultiVector& x, 
				 NOX::Abstract::MultiVector::DenseMatrix& y, 
				 const Epetra_MultiVector& augMultiVec) const;

    protected:

      //! Builds extended domain, range maps
      void buildExtendedMap(const Epetra_BlockMap& map,
			    Epetra_Map*& extMapPtr, bool buildImporter, 
			    bool haveParam);
	
      //! Converts a block map to an equivalent point map
      int blockMap2PointMap(const Epetra_BlockMap& BlockMap,
			    Epetra_Map*& PointMap) const;

    private:

      //! Private to prohibit copying
      AugmentedOp(const AugmentedOp&);

      //! Private to prohibit copying
      AugmentedOp& operator=(const AugmentedOp&);

    protected:

      //! LOCA global data object
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Label for operator
      std::string label;
    
      //! Stores operator representing \f$J\f$
      Teuchos::RCP<Epetra_Operator> jacOperator;

       //! Stores underlying domain map
      const Epetra_BlockMap& underlyingMap;

      //! Stores comm
      const Epetra_Comm& underlyingComm;

      //! Local map for generating Epetra matrices
      Epetra_LocalMap localMap;
  
      //! Stores pointer to a multivector
      Teuchos::RCP<const Epetra_MultiVector> a;

      //! Stores pointer to b multivector
      Teuchos::RCP<const Epetra_MultiVector> b;

      //! Stores c matrix
      Epetra_MultiVector c;

      //! Stores underlying vector local length
      int underlyingLength;

      //! Number of constraints
      int numConstraints;

      //! Stores extended domain map
      Epetra_Map* extendedMapPtr;

      //! Stores extended turning point map for importing param component
      Epetra_Map* extendedImportMapPtr;

      //! Stores importer object for importing param component
      Epetra_Import* extendedImporter;

      //! Stores imported input multivector
      Epetra_MultiVector* importedInput;

      //! Stores parameter component of result multivector
      Epetra_MultiVector* result_y;

      //! Stores temporary multivector used in ApplyInverse()
      Epetra_MultiVector* tmp;

      //! Flag indicating whether we have the parameter component
      bool haveParamComponent;

      //! Flag indicating whether to use transpose of operator
      bool useTranspose;

      //! LAPACK Wrappers
      Teuchos::LAPACK<int,double> dlapack;

    };
  } // namespace Epetra
} // namespace LOCA

#endif 
