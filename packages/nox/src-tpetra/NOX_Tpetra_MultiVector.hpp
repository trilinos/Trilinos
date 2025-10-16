// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_TPETRA_MULTIVECTOR_HPP
#define NOX_TPETRA_MULTIVECTOR_HPP

#include "Tpetra_Map_fwd.hpp"
#include "Tpetra_MultiVector_fwd.hpp"
#include "Tpetra_Vector_fwd.hpp"

#include "NOX_Abstract_ImplicitWeighting.H"
#include "NOX_Abstract_MultiVector.H"

#include "NOX_Tpetra_Vector.hpp"

namespace NOX::Tpetra
{

/**
* @brief Implementation of @ref NOX::Abstract::MultiVector for @ref ::Tpetra multivectors.
*
* This class is templated on the same quadruple of template parameters as the
* wrapped @ref ::Tpetra::MultiVector.
*
* In particular, this class encapsulates a shared pointer holding a @ref ::Tpetra::MultiVector.
* In @ref NOX, a @c MultiVector object also has shared pointers to objects for using weights
* in norm and inner product computations. We follow the design of @ref NOX::Thyra::MultiVector
* using @ref NOX::Abstract::ImplicitWeighting. As in @ref NOX::Thyra::MultiVector, this class
* also holds a container of @ref NOX::Tpetra::Vector objects that encapsulate shared pointers
* with views of the individual columns of the wrapped @ref ::Tpetra::MultiVector.
*/
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
class MultiVector : public NOX::Abstract::MultiVector,
                    public NOX::Abstract::ImplicitWeighting
{
public:
    //! Type of the wrapped @ref ::Tpetra::MultiVector.
    using multivector_type = ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    //! Type of a @ref NOX::Tpetra::Vector instantiated with the same quadruple of template parameters.
    using nox_tpetra_vector_type = Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    /**
     * @brief Type of the @ref ::Tpetra::Vector wrapped in a @ref NOX::Tpetra::Vector
     *        instantiated with the same quadruple of template parameters.
     */
    using vector_type = typename nox_tpetra_vector_type::vector_type;

public:
    //! Constructor. Analogous to @ref Vector::Vector(const Teuchos::RCP<vector_type>&).
    MultiVector(const Teuchos::RCP<multivector_type>& source);

    //! Constructor. Analogous to @ref Vector::Vector(const vector_type&, CopyType).
    MultiVector(const multivector_type& source, CopyType type = DeepCopy);

    //! Copy constructor. Analogous to @ref Vector::Vector(const Vector&, CopyType).
    MultiVector(const MultiVector& source, CopyType type = DeepCopy);

    //! Returns the encapsulated shared pointer holding the @ref ::Tpetra::MultiVector.
    Teuchos::RCP<multivector_type> getTpetraMultiVector() const;

    //! Copy assignment operator. Analogous to @ref Vector::operator=(const Teuchos::RCP<const multivector_type>&).
    Abstract::MultiVector& operator=(const Teuchos::RCP<const multivector_type>& source);

    //! Copy assignment operator. Analogous to @ref Vector::operator=(const MultiVector&).
    Abstract::MultiVector& operator=(const MultiVector& source);

    //! Copy assignment operator. Analogous to @ref Vector::operator=(Abstract::MultiVector&).
    Abstract::MultiVector& operator=(const Abstract::MultiVector& source) override;

    //! Analogous to @ref Vector::init().
    Abstract::MultiVector& init(double gamma) override;

    //! Analogous to @ref Vector::random(bool, int).
    Abstract::MultiVector& random(bool useSeed = false, int seed = 1) override;

    //! Analogous to @ref Vector::scale(double).
    Abstract::MultiVector& scale(double gamma) override;

    //! Analogous to @ref Vector::update(double, const MultiVector&, double).
    Abstract::MultiVector& update(double alpha, const MultiVector& a, double gamma = 0.);

    //! Analogous to @ref Vector::update(double, const Abstract::MultiVector&, double).
    Abstract::MultiVector& update(double alpha, const Abstract::MultiVector& a, double gamma = 0.) override;

    //! Analogous to @ref Vector::update(double, const MultiVector&, double, const MultiVector&, double).
    Abstract::MultiVector& update(
        double alpha,
        const MultiVector& a,
        double beta,
        const MultiVector& b,
        double gamma = 0.
    );

    //! Analogous to @ref Vector::update(double, const Abstract::MultiVector&, double, const Abstract::MultiVector&, double).
    Abstract::MultiVector& update(
        double alpha,
        const Abstract::MultiVector& a,
        double beta,
        const Abstract::MultiVector& b,
        double gamma = 0.
    ) override;

    /**
     * @brief Updates the elements of the wrapped @ref ::Tpetra::MultiVector with the elements of the
     *        @ref ::Tpetra::Vector objects wrapped in the referenced @ref NOX::Tpetra::MultiVector @p a
     *        and the elements of the referenced @ref Teuchos::SerialDenseMatrix @p b.
     *
     * This operation is such that @c x becomes <tt> alpha * a * op(b) + gamma * x </tt>, with @c x this vector
     * and @c op(b) equal to itself if <tt> transb = Teuchos::NO_TRANS </tt> and equal to the transpose
     * of @p b if <tt> transb = Teuchos::TRANS </tt>.
     */
    Abstract::MultiVector& update(
        Teuchos::ETransp transb,
        double alpha,
        const MultiVector& a,
        const Abstract::MultiVector::DenseMatrix& b,
        double gamma = 0.
    );

    //! See above.
    Abstract::MultiVector& update(
        Teuchos::ETransp transb,
        double alpha,
        const Abstract::MultiVector& a,
        const Abstract::MultiVector::DenseMatrix& b,
        double gamma = 0.
    ) override;

    //! Analogous to @ref Vector::clone(CopyType).
    Teuchos::RCP<Abstract::MultiVector> clone(CopyType type = DeepCopy) const override;

    /**
     * @brief Creates a new @ref NOX::Tpetra::MultiVector that wraps a @ref ::Tpetra::MultiVector
     *        constructed in newly allocated memory with @p numVecs columns and the map of the
     *        encapsulated @ref ::Tpetra::MultiVector.
     */
    Teuchos::RCP<Abstract::MultiVector> clone(int numVecs) const override;

    //! Compute the norms of the individual columns of the wrapped @ref ::Tpetra::MultiVector.
    void norm(
        std::vector<double>& result,
        Abstract::Vector::NormType type = Abstract::Vector::TwoNorm
    ) const override;

    //! Returns the number of vectors held in the wrapped @ref ::Tpetra::MultiVector.
    int numVectors() const override;
 
    /**
     *  @brief Copies the elements of the @c indices.size() columns of the @ref ::Tpetra::MultiVector
     *         wrapped in @p source to the elements of the columns of the wrapped @ref ::Tpetra::MultiVector
     *         indexed by the indices in @p indices.
     */
    Abstract::MultiVector& setBlock(
        const MultiVector& source,
        const std::vector<int>& indices
    );

    //! See above.
    Abstract::MultiVector& setBlock(
        const Abstract::MultiVector& source,
        const std::vector<int>& indices
    ) override;

    /**
     * @brief Appends the columns of the @ref ::Tpetra::MultiVector wrapped in @p source to the columns
     *        of the wrapped @ref ::Tpetra::MultiVector.
     *
     * This function constructs a @ref ::Tpetra::MultiVector in newly allocated memory to hold
     * the augmented set of columns.
     */
    Abstract::MultiVector& augment(const MultiVector& source);

    //! See above.
    Abstract::MultiVector& augment(const Abstract::MultiVector& source) override;

    /**
     * @brief Returns a reference to a @ref NOX::Tpetra::Vector that encapsulates a shared pointer
     *        holding a @ref ::Tpetra::Vector object with a view of the @p i th column of
     *        the wrapped @ref ::Tpetra::MultiVector.
     *
     * As in @ref NOX::Thyra::MultiVector, the weight vector, if any, is not passed to the referenced
     * @ref NOX::Tpetra::Vector.
     */
    NOX::Abstract::Vector& operator[](int i) override;

    //! See above.
    const NOX::Abstract::Vector& operator[](int i) const override;

    /**
     * @brief Creates a new @ref NOX::Tpetra::MultiVector object with @c indices.size() columns whose
     *        elements are copies of the elements of the columns of the wrapped @ref ::Tpetra::MultiVector
     *        indexed by @p indices.
     *
     * As in @ref NOX::Thyra::MultiVector, the weight vector, if any, is not passed to the new
     * @ref NOX::Tpetra::MultiVector.
     */
    Teuchos::RCP<Abstract::MultiVector> subCopy(const std::vector<int>& indices) const override;

    /**
     * @brief Creates a new @ref NOX::Tpetra::MultiVector object with @c indices.size() columns that
     *        are views of the columns of the wrapped @ref ::Tpetra::MultiVector indexed by @p indices.
     *
     * As in @ref NOX::Thyra::MultiVector, the weight vector, if any, is not passed to the new
     * @ref NOX::Tpetra::MultiVector.
     */
    Teuchos::RCP<Abstract::MultiVector> subView(const std::vector<int>& indices) const override;

    /**
     * @brief Matrix-matrix product.
     *
     * This operation returns <tt> b = alpha * transp(y) * x </tt>, where @c x is this multivector.
     */
    void multiply(
        double alpha,
        const MultiVector& y,
        Abstract::MultiVector::DenseMatrix& b
    ) const;

    //! See above.
    void multiply(
        double alpha,
        const Abstract::MultiVector& y,
        Abstract::MultiVector::DenseMatrix& b
    ) const override;

    //! Prints the wrapped @ref ::Tpetra::MultiVector.
    void print(std::ostream& stream) const override;

    //! Analogous to @ref Vector::setWeightVector().
    void setWeightVector(const Teuchos::RCP<const vector_type>& weightVec_);

    //! Analogous to @ref Vector::hasWeightVector().
    bool hasWeightVector() const;

    //! Analogous to @ref Vector::getWeightVector().
    Teuchos::RCP<const vector_type> getWeightVector() const;

    //! Returns the number of rows in the wrapped @ref ::Tpetra::MultiVector.
    NOX::size_type length() const override;

    //! Analogous to @ref Vector::getImplicitWeighting().
    bool getImplicitWeighting() const override;

    //! Analogous to @ref Vector::setImplicitWeighting().
    void setImplicitWeighting(bool do_implicit_weighting_) override;

private:
    //! Shared pointer holding the @ref ::Tpetra::MultiVector wrapped in this object.
    Teuchos::RCP<multivector_type> tpetraMultiVec;

    /**
     * @brief Container of @ref NOX::Tpetra::Vector objects that encapsulate shared pointers
     *        holding @ref ::Tpetra::Vector objects with views of the individual columns of
     *        the wrapped @ref ::Tpetra::MultiVector.
     *
     * @note The @ref NOX::Tpetra::Vector objects are constructed lazily, i.e., only when the
     *       corresponding column is accessed for the first time.
     */
    mutable std::vector<std::optional<nox_tpetra_vector_type>> noxTpetraVecs;

    //! Shared pointer holding the @ref ::Tpetra::Vector used for weighting inner products and norms.
    Teuchos::RCP<const vector_type> weightVec;

    /**
     * @brief Flag that controls whether the @ref ::Tpetra::Vector for weighting
     *        inner products and norms (if not null) is used.
     */
    bool do_implicit_weighting = false;
};

} // namespace NOX::Tpetra

#endif // NOX_TPETRA_MULTIVECTOR_HPP
