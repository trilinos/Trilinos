//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER
                                                                                
#ifndef LOCA_EPETRA_BORDEREDOP_H
#define LOCA_EPETRA_BORDEREDOP_H

#include "Epetra_Operator.h"	// base class
#include "NOX_Common.H"         // for string data member

// Forward declarations
class Epetra_Vector;
class Epetra_Map;

namespace LOCA {

  namespace Epetra {

    /*!
     * \brief %Epetra operator representing a \f$n+1\f$ bordered matrix.
     */
    /*!
     * The %LOCA::Epetra::BorderedOp is an Epetra_Operator representing the
     * \f$n+1\f$ bordered matrix
     * \f[
     *     \begin{bmatrix}
     *        J & a \\
     *        b^T & 0
     *     \end{bmatrix}
     * \f]
     * where \f$J\f$ is an Epetra_Operator representing an \f$n\times n\f$
     * matrix, and \f$a\f$ and \f$b\f$ are length \f$n\f$ Epetra_Vector's
     * (The case where \f$a\f$ and/or \f$b\f$ are multi-vectors is not yet
     * supported).  It is assumed the Epetra_Map's for \f$a\f$, \f$b\f$, and
     * \f$J\f$ are the same and the corresponding map for the bordered
     * matrix is constructed from this map by storing the additional component
     * on processor 0.  The buildEpetraExtendedVec() method can be used to
     * construct an Epetra_Vector using this map and a supplied length \f$n\f$
     * Epetra_Vector and a scalar, while setEpetraExtendedVec() splits an
     * extended vector into its length \f$n\f$ and scalar components.  The 
     * Apply() method performs the \f$n+1\times n+1\f$ matrix multiplication
     * while ApplyInverse() uses a block-elimination algorithm to compute
     * the inverse using the ApplyInverse() method of the underlying operator
     * \f$J\f$.  In this way, linear systems of the form
     * \f[
     *     \begin{bmatrix}
     *        J & a \\
     *        b^T & 0
     *     \end{bmatrix}
     *     \begin{bmatrix}
     *       x \\
     *       y
     *     \end{bmatrix} =
     *     \begin{bmatrix}
     *       c \\
     *       d
     *     \end{bmatrix}
     * \f]
     * can be solved in a matrix-free mode using the Apply() method.  This 
     * operator can also represent a preconditioner of the form
     * \f[
     *     \begin{bmatrix}
     *        M & a \\
     *        b^T & 0
     *     \end{bmatrix}
     * \f]
     * using the ApplyInvese() method, where \f$M\f$ is a preconditioner 
     * for \f$J\f$.  Note that if \f$J\f$ is nearly singular, the 
     * preconditioner should not be too good because otherwise the
     * preconditining operation represented by ApplyInverse() becomes unstable.
     */
    class BorderedOp : public Epetra_Operator {
      
    public:

      //! Constructor 
      /*!
       * Builds the bordered operator using the supplied operator \em jac 
       * and Epetra_Vector's \em a and \em b.  It is assumed \em a, \em b,
       * and \em jac all have the same map.
       */
      BorderedOp(Epetra_Operator& jac, 
		 const Epetra_Vector& a, 
		 const Epetra_Vector& b);

      //! Destructor
      virtual ~BorderedOp();

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
       *        J & a \\
       *        b^T & 0
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *       x \\
       *       y
       *     \end{bmatrix} =
       *     \begin{bmatrix}
       *       Jx + ya \\
       *       b^T x
       *     \end{bmatrix}
       * \f]
       * or its transpose if UseTranpose() is \em true.  Currently, only the
       * case where \em Input and \em Result are Epetra_Vector's is supported.
       */
      virtual int Apply(const Epetra_MultiVector& Input, 
			Epetra_MultiVector& Result) const;

      /*!
       * \brief Returns the result of a Epetra_Operator inverse applied to 
       * an Epetra_MultiVector X in Y.
       */
      /*!
       * Solves the extended system
       * \f[
       *     \begin{bmatrix}
       *        J & a \\
       *        b^T & 0
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *       x \\
       *       y
       *     \end{bmatrix} =
       *     \begin{bmatrix}
       *       c \\
       *       d
       *     \end{bmatrix}
       * \f]
       * using the following block-elimination algorithm:  
       * \f[
       *    \tilde{c} = J^{-1} c, \quad \tilde{a} = J^{-1} a, \quad y = \frac{b^T\tilde{c} - d}{b^T\tilde{a}}, \quad x = \tilde{c} - y\tilde{a}
       * \f]
       * Currently, only the case where \em Input and \em Result are 
       * Epetra_Vector's is supported.  If UseTranpose() is \em true, the
       * tranpose of the system is solved.
       */
      virtual int ApplyInverse(const Epetra_MultiVector& X, 
			       Epetra_MultiVector& Y) const;
  
      //! Returns the infinity norm of the bordered matrix.
      /*!
       * This is defined only if %NormInf() of the underlying operator \f$J\f$
       * is defined and is given by 
       * \f$\|J\|_\infty+\|a\|_\infty+\|b\|_\infty\f$.
       */
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
       * \brief Builds an extended  vector from components
       */
      /*!
       * Builds an extended vector using the map representing the bordered
       * matrix.  If \em doCopy is \em true, the contents of \em x are
       * copied into the extended vector, otherwise only space for the 
       * extended vector is created.
       */
      virtual Epetra_Vector* buildEpetraExtendedVec(Epetra_Vector& x,
						    double p, bool doCopy) const;

      /*!
       * \brief Sets components from extended vector
       */
      /*!
       * Splits the extended vector \em extVec into components \em x and
       * \em p by copying values out of \em extVec.
       */
      virtual void setEpetraExtendedVec(Epetra_Vector& x, double& p, 
					Epetra_Vector& extVec) const;

    protected:

      //! Builds extended domain, range maps
      void buildExtendedMap(const Epetra_BlockMap& map,
			    Epetra_Map*& extMapPtr, bool buildImporter, 
			    bool haveParam);
	
      //! Converts a block map to an equivalent point map
      int blockMap2PointMap(const Epetra_BlockMap& BlockMap,
			    Epetra_Map*& PointMap) const;

    protected:

      //! Label for operator
      string label;
    
      //! Stores operator representing \f$J\f$
      Epetra_Operator& jacOperator;
  
      //! Stores pointer to a vector
      const Epetra_Vector* aVecPtr;

      //! Stores pointer to b vector
      const Epetra_Vector* bVecPtr;

      //! Stores underlying domain map
      const Epetra_BlockMap* underlyingMapPtr;

      //! Stores comm
      const Epetra_Comm* underlyingCommPtr;

      //! Stores underlying vector local length
      int underlyingLength;

      //! Stores extended domain map
      Epetra_Map* extendedMapPtr;

      //! Stores extended turning point map for importing param component
      Epetra_Map* extendedImportMapPtr;

      //! Stores importer object for importing param component
      Epetra_Import* extendedImporter;

      //! Flag indicating whether we have the parameter component
      bool haveParamComponent;

      //! Flag indicating whether to use transpose of operator
      bool useTranspose;

    };
  } // namespace Epetra
} // namespace LOCA

#endif 
