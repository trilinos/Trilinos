// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_MULTIVECTOR_H
#define NOX_MULTIVECTOR_H

#include "NOX_Abstract_MultiVector.H" // base class
#include "Teuchos_RCP.hpp"
#include <vector>

namespace NOX {

  /*!
   * \brief Default implementation for NOX::Abstract::MultiVector using
   * an array of NOX::Abstract::MultiVector's.
   */
  class MultiVector : public virtual NOX::Abstract::MultiVector {

  public:

    /*!
     * \brief Create %MultiVector with \c numVecs columns out of a single
     * NOX::Abstract::Vector.
     */
    MultiVector(const NOX::Abstract::Vector& v, int numVecs = 1,
        NOX::CopyType type = NOX::DeepCopy);

    //! Create %MultiVector out of array of NOX::Abstract::Vector's
    MultiVector(const NOX::Abstract::Vector* const* vs, int numVecs,
        NOX::CopyType type = NOX::DeepCopy);

    //! Copy constructor
    MultiVector(const MultiVector& source,
        NOX::CopyType type = NOX::DeepCopy);

    //! Destructor.
    virtual ~MultiVector();

    //@{ \name Initialization methods.

    //! Initialize every element of this multi-vector with \c gamma.
    virtual NOX::Abstract::MultiVector&
    init(double gamma);

    //! Initialize each element of this multi-vector with a random value
    virtual NOX::Abstract::MultiVector&
    random(bool useSeed = false, int seed = 1);

    //! Copy source multi-vector \c source into this multi-vector.
    virtual NOX::Abstract::MultiVector&
    operator=(const NOX::Abstract::MultiVector& source);

    //! Copy source multi-vector \c source into this multi-vector.
    virtual NOX::Abstract::MultiVector&
    operator=(const NOX::MultiVector& source);

    /*!
     * \brief Copy the vectors in \c source to a set of vectors in \c *this.
     * The \c index.size() vectors in \c source are copied to a subset of
     * vectors in \c *this indicated by the indices given in \c index.
     */
    virtual NOX::Abstract::MultiVector&
    setBlock(const NOX::Abstract::MultiVector& source,
         const std::vector<int>& index);
    virtual NOX::Abstract::MultiVector&
    setBlock(const NOX::MultiVector& source, const std::vector<int>& index);

    //! Append the vectors in \c source to \c *this
    virtual NOX::Abstract::MultiVector&
    augment(const NOX::Abstract::MultiVector& source);
    virtual NOX::Abstract::MultiVector&
    augment(const NOX::MultiVector& source);

    //@}

    //@{ \name Accessor methods

    /*!
     * \brief Return a reference to the i-th column of the multivector as
     * an abstract vector.
     */
    virtual NOX::Abstract::Vector& operator [] (int i);

    /*!
     * \brief Return a const reference to the i-th column of the multivector as
     * an abstract vector.
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
    virtual NOX::Abstract::MultiVector&
    update(double alpha, const NOX::MultiVector& a,
       double gamma = 0.0);

    /*!
     * \brief Compute x = (alpha * a) + (beta * b) + (gamma * x) where
     * \c a and \c b are multi-vectors and \c x = \c *this.
     */
    virtual NOX::Abstract::MultiVector&
    update(double alpha, const NOX::Abstract::MultiVector& a,
       double beta, const NOX::Abstract::MultiVector& b,
       double gamma = 0.0);
    virtual NOX::Abstract::MultiVector&
    update(double alpha, const NOX::MultiVector& a,
       double beta, const NOX::MultiVector& b,
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
    virtual NOX::Abstract::MultiVector&
    update(Teuchos::ETransp transb, double alpha,
       const NOX::MultiVector& a,
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
    clone(NOX::CopyType type = NOX::DeepCopy) const;

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
     * \brief Creates a new multi-vector with \c index.size() columns that
     * shares the columns of \c *this given by \c index.
     */
    virtual Teuchos::RCP<NOX::Abstract::MultiVector>
    subView(const std::vector<int>& index) const;

    //@}

    //@{ \name Norms.

    //! Norm.
    virtual void
    norm(std::vector<double>& result,
     NOX::Abstract::Vector::NormType type =
     NOX::Abstract::Vector::TwoNorm) const;

    //@}

    //@{ \name Matrix products.

    /*!
     * \brief Computes the matrix-matrix product \f$\alpha * y^T * (*this)\f$
     */
    virtual void multiply(double alpha, const NOX::Abstract::MultiVector& y,
              NOX::Abstract::MultiVector::DenseMatrix& b) const;
    virtual void multiply(double alpha, const NOX::MultiVector& y,
              NOX::Abstract::MultiVector::DenseMatrix& b) const;

    //@}

    //! Return the length of multi-vector.
    virtual NOX::size_type length() const;

    //! Return the number of vectors in the multi-vector.
    virtual int numVectors() const;

    //! Print the vector.  This is meant for debugging purposes only
    virtual void print(std::ostream& stream) const;

  protected:

    //! Stores pointers to NOX::Abstract::Vector's
    std::vector< Teuchos::RCP<NOX::Abstract::Vector> > vecs;

    //! Create %MultiVector holding \c numvecs vectors
    MultiVector(int numvecs);

    //! Checks whether an index is valid.  Throws an error if invalid
    void checkIndex(int idx) const;

    /*!
     * \brief Checks whether the given multivector size is compatible with this
     * multivector.  Throws an error if invalid.
     */
    void checkSize(int sz) const;

  }; // class MultiVector
} // namespace NOX

#endif // NOX_MULTIVECTOR_H
