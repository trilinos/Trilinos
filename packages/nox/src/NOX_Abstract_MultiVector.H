// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_ABSTRACT_MULTIVECTOR_H
#define NOX_ABSTRACT_MULTIVECTOR_H

#include "NOX_Common.H"  // for NOX_Config.h

#include "Teuchos_SerialDenseMatrix.hpp"

#include "NOX_Abstract_Vector.H"

namespace NOX {

  namespace Abstract {

    /*!
     * \brief Abstract interface for multi-vectors used by %NOX.
     */
    class MultiVector {

    public:

      //! Typename of dense matrices
      typedef Teuchos::SerialDenseMatrix<int,double> DenseMatrix;

      //! Default constructor.  Does nothing.
      MultiVector() {};

      //! Destructor.  Does nothing.
      virtual ~MultiVector() {};

      //@{ \name Initialization methods.

      //! Initialize every element of this multi-vector with \c gamma.
      virtual NOX::Abstract::MultiVector&
      init(double gamma) = 0;

      //! Initialize each element of this multi-vector with a random value
      virtual NOX::Abstract::MultiVector&
      random(bool useSeed = false, int seed = 1) = 0;

      //! Copy source multi-vector \c source into this multi-vector.
      virtual NOX::Abstract::MultiVector&
      operator=(const NOX::Abstract::MultiVector& source) = 0;

      /*!
       * \brief Copy the vectors in \c source to a set of vectors in \c *this.
       * The \c index.size() vectors in \c source are copied to a subset of
       * vectors in \c *this indicated by the indices given in \c index.
       */
      virtual NOX::Abstract::MultiVector&
      setBlock(const NOX::Abstract::MultiVector& source,
           const std::vector<int>& index) = 0;

      //! Append the vectors in \c source to \c *this
      virtual NOX::Abstract::MultiVector&
      augment(const NOX::Abstract::MultiVector& source) = 0;

      //@}

      //@{ \name Accessor methods

      /*!
       * \brief Return a reference to the i-th column of the multivector as
       * an abstract vector.
       */
      /*!
       * It is assumed that any implementation of this method returns a
       * vector that has a view of the underlying data in the multivector
       */
      virtual NOX::Abstract::Vector& operator [] (int i) = 0;

      /*!
       * \brief Return a const reference to the i-th column of the
       * multivector as an abstract vector.
       */
      /*!
       * It is assumed that any implementation of this method returns a
       * vector that has a view of the underlying data in the multivector.
       */
      virtual const NOX::Abstract::Vector& operator [] (int i) const = 0;

      //@}

      //@{ \name Update methods.

      //! Scale each element of this multivector by \c gamma.
      virtual NOX::Abstract::MultiVector&
      scale(double gamma) = 0;

      /*!
       * \brief Compute x = (alpha * a) + (gamma * x) where \c a is a
       * multi-vector and \c x = \c *this.
       */
      virtual NOX::Abstract::MultiVector&
      update(double alpha, const NOX::Abstract::MultiVector& a,
         double gamma = 0.0) = 0;

      /*!
       * \brief Compute x = (alpha * a) + (beta * b) + (gamma * x) where
       * \c a and \c b are multi-vectors and \c x = \c *this.
       */
      virtual NOX::Abstract::MultiVector&
      update(double alpha, const NOX::Abstract::MultiVector& a,
         double beta, const NOX::Abstract::MultiVector& b,
         double gamma = 0.0) = 0;

      /*!
       * \brief Compute x = (alpha * a * op(b)) + (gamma * x) where
       * \c a is a multivector, \c b is a dense matrix, \c x = \c *this,
       * and op(b) = b if transb = Teuchos::NO_TRANS and op(b) is b transpose
       * if transb = Teuchos::TRANS.
       */
      virtual NOX::Abstract::MultiVector&
      update(Teuchos::ETransp transb, double alpha,
         const NOX::Abstract::MultiVector& a,
         const DenseMatrix& b, double gamma = 0.0) = 0;

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
      clone(NOX::CopyType type = NOX::DeepCopy) const = 0;

      //! Creates a new multi-vector with \c numvecs columns
      virtual Teuchos::RCP<NOX::Abstract::MultiVector>
      clone(int numvecs) const = 0;

      /*!
       * \brief Creates a new multi-vector with \c index.size() columns whose
       * columns are copies of the columns of \c *this given by \c index.
       */
      virtual Teuchos::RCP<NOX::Abstract::MultiVector>
      subCopy(const std::vector<int>& index) const = 0;

      /*!
       * \brief Creates a new multi-vector with \c ndex.size() columns that
       * shares the columns of \c *this given by \c index.
       */
      virtual Teuchos::RCP<NOX::Abstract::MultiVector>
      subView(const std::vector<int>& index) const = 0;

      //@}

      //@{ \name Norms.

      //! Norm.
      virtual void
      norm(std::vector<double>& result,
       NOX::Abstract::Vector::NormType type =
       NOX::Abstract::Vector::TwoNorm) const = 0;

      //@}

      //@{ \name Matrix products.

      /*!
       * \brief Computes the matrix-matrix product \f$\alpha * y^T * (*this)\f$
       */
      virtual void multiply(double alpha, const NOX::Abstract::MultiVector& y,
                DenseMatrix& b) const = 0;

      //@}

      //! Return the length of multi-vector.
      virtual NOX::size_type length() const = 0;

      //! Return the number of vectors in the multi-vector.
      virtual int numVectors() const = 0;

      //! Print the vector.  This is meant for debugging purposes only
      virtual void print(std::ostream& stream) const = 0;

    }; // class MultiVector
  } // namespace Abstract
} // namespace NOX

#endif // NOX_ABSTRACT_MULTIVECTOR_H
