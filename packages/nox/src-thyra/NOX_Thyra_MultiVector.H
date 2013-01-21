// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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

#ifndef NOX_THYRA_MULTIVECTOR_H
#define NOX_THYRA_MULTIVECTOR_H

#include "NOX_Abstract_MultiVector.H" // base class
#include "Teuchos_RCP.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include <vector>

namespace NOX {
  namespace Thyra {
    class Vector;
  }
}

// The NOX::Thyra::MultiVector is not currently used because 
// Thyra::MultiVectorBase views are not compatible with the way LOCA
// uses views

namespace NOX {

  namespace Thyra {

    //! Implementation of NOX::Abstract::MultiVector for %Thyra multi-vectors.
    class MultiVector : public virtual NOX::Abstract::MultiVector {

    public:			
      
      //! \brief Constructor that creates a VIEW of the Thyra multivector.
      MultiVector(
       const Teuchos::RCP< ::Thyra::MultiVectorBase<double> >& source);

      //! Construct from a given Thyra multivector
      MultiVector(const ::Thyra::MultiVectorBase<double>& source);

      //! Copy constructor.
      MultiVector(const NOX::Thyra::MultiVector& source, 
		  NOX::CopyType type = NOX::DeepCopy);

      //! Destruct MultiVector.
      ~MultiVector();

      //@{ \name Access to underlying Thyra multivector.

      //! Get RCP to underlying Thyra vector.
      virtual Teuchos::RCP< ::Thyra::MultiVectorBase<double> >
      getThyraMultiVector();

      //! Get const RCP to underlying Thyra vector.
      virtual Teuchos::RCP<const ::Thyra::MultiVectorBase<double> >
      getThyraMultiVector() const;
 
      //@}

      //@{ \name Initialization methods.

      //! Initialize every element of this multi-vector with \c gamma.
      virtual NOX::Abstract::MultiVector& 
      init(double value);

      //! Initialize each element of this multi-vector with a random value
      virtual NOX::Abstract::MultiVector& 
      random(bool useSeed = false, int seed = 1);

      //! Copy source multi-vector \c source into this multi-vector.
      virtual NOX::Abstract::MultiVector& 
      operator=(const NOX::Abstract::MultiVector& source);

      /*! 
       * \brief Copy the vectors in \c source to a set of vectors in \c *this.
       * The \c index.size() vectors in \c source are copied to a subset of 
       * vectors in \c *this indicated by the indices given in \c index.
       */
      virtual NOX::Abstract::MultiVector&
      setBlock(const NOX::Abstract::MultiVector& source, 
	       const std::vector<int>& index);

      //! Append the vectors in \c source to \c *this
      virtual NOX::Abstract::MultiVector&
      augment(const NOX::Abstract::MultiVector& source);

      //@}

      //@{ \name Accessor methods

      /*!
       * \brief Return a reference to the i-th column of the multivector as
       * an abstract vector.  
       */
      virtual NOX::Abstract::Vector& operator [] (int i);
    
      /*!
       * \brief Return a const reference to the i-th column of the 
       * multivector as an abstract vector.  
       */
      virtual const NOX::Abstract::Vector& operator [] (int i) const;

      //@}

      //@{ \name Update methods.

      //! Scale each element of this multivector by \c gamma.
      virtual NOX::Abstract::MultiVector& 
      scale(double gamma);

      /*! 
       * \brief Compute x = (alpha * a) + (gamma * x) where \c a is a 
       * multi-vector and \c x = \c *this.
       */
      virtual NOX::Abstract::MultiVector& 
      update(double alpha, const NOX::Abstract::MultiVector& a, 
	     double gamma = 0.0);

      /*! 
       * \brief Compute x = (alpha * a) + (beta * b) + (gamma * x) where
       * \c a and \c b are multi-vectors and \c x = \c *this.
       */
      virtual NOX::Abstract::MultiVector& 
      update(double alpha, const NOX::Abstract::MultiVector& a, 
	     double beta, const NOX::Abstract::MultiVector& b,
	     double gamma = 0.0);

      /*! 
       * \brief Compute x = (alpha * a * b) + (gamma * x) where 
       * \c a is a multivector, \c b is a dense matrix, \c x = \c *this,
       * and op(b) = b if transb = Teuchos::NO_TRANS and op(b) is b transpose
       * if transb = Teuchos::TRANS.
       */
      virtual NOX::Abstract::MultiVector& 
      update(Teuchos::ETransp transb, double alpha, 
	     const NOX::Abstract::MultiVector& a, 
	     const NOX::Abstract::MultiVector::DenseMatrix& b, 
	     double gamma = 0.0);

      //@}

      //@{ \name Creating new MultiVectors. 

      /*! 
       * \brief Create a new %Vector of the same underlying type by
       * cloning "this", and return a pointer to the new vector.  
       *
       * If type is NOX::DeepCopy, then we need to create an exact replica
       * of "this". Otherwise, if type is NOX::ShapeCopy, we need only
       * replicate the shape of "this". Note that there is 
       * <em>no assumption</em> that a vector created by ShapeCopy is 
       * initialized to zeros.
       * 
       * \return Pointer to newly created vector or NULL if clone is not 
       * supported. 
       */
      virtual Teuchos::RCP<NOX::Abstract::MultiVector> 
      clone(CopyType type = DeepCopy) const;

      //! Creates a new multi-vector with \c numvecs columns
      virtual Teuchos::RCP<NOX::Abstract::MultiVector>
      clone(int numvecs) const;

      /*! 
       * \brief Creates a new multi-vector with \c index.size() columns whose 
       * columns are copies of the columns of \c *this given by \c index.
       */
      virtual Teuchos::RCP<NOX::Abstract::MultiVector> 
      subCopy(const std::vector<int>& index) const;

      /*! 
       * \brief Creates a new multi-vector with \c ndex.size() columns that 
       * shares the columns of \c *this given by \c index.
       */
      virtual Teuchos::RCP<NOX::Abstract::MultiVector> 
      subView(const std::vector<int>& index) const;

      //@}

      //@{ \name Norms.

      //! Norm.
      virtual void
      norm(std::vector<double>& result, NOX::Abstract::Vector::NormType type = 
	   NOX::Abstract::Vector::TwoNorm) const;

      //@}

      //@{ \name Matrix products.

      /*! 
       * \brief Computes the matrix-matrix product \f$\alpha * y^T * (*this)\f$
       */
      virtual void multiply(double alpha, const NOX::Abstract::MultiVector& y,
			    NOX::Abstract::MultiVector::DenseMatrix& b) const;
      
      //@}

      //! Return the length of multi-vector.
      virtual int length() const;

      //! Return the number of vectors in the multi-vector.
      virtual int numVectors() const;

      //! Print the vector.  This is meant for debugging purposes only 
      virtual void print(std::ostream& stream) const;

    protected:

      //! Check whether an index array is contiguous
      bool isContiguous(const std::vector<int>& index) const;

    protected:

       //! Pointer to petra vector owned by this object
      Teuchos::RCP< ::Thyra::MultiVectorBase<double> > thyraMultiVec;

      //! NOX::Thyra::Vector's for each column of the multivector
      /*!
       * Each Thyra_Vector in the NOX::Thyra::Vector has a view into a 
       * column of the multivector and get filled in as needed by operator[].
       */
      mutable std::vector< Teuchos::RCP<NOX::Thyra::Vector> > noxThyraVectors;

    };
  } // namespace Thyra
} // namespace NOX

#endif // NOX_THYRA_MULTIVECTOR_H
