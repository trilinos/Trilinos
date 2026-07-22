// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_TPETRA_VECTOR_HPP
#define NOX_TPETRA_VECTOR_HPP

#include "Tpetra_Map_fwd.hpp"
#include "Tpetra_MultiVector_fwd.hpp"
#include "Tpetra_Vector_fwd.hpp"

#include "NOX_Abstract_ImplicitWeighting.H"
#include "NOX_Abstract_Vector.H"

namespace NOX::Tpetra
{

 /**
  * @brief Implementation of @ref NOX::Abstract::Vector for @ref ::Tpetra vectors.
  *
  * This class is templated on the same quadruple of template parameters as the
  * wrapped @ref ::Tpetra::Vector.
  *
  * In particular, this class encapsulates a shared pointer holding a @ref ::Tpetra::Vector.
  * In @ref NOX, a @c Vector object also has shared pointers to objects for using weights
  * in norm and inner product computations. We follow the design of @ref NOX::Thyra::Vector
  * using @ref NOX::Abstract::ImplicitWeighting.
  */
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
class Vector : public Abstract::Vector,
               public NOX::Abstract::ImplicitWeighting
{
public:
    //! Type of the wrapped @ref ::Tpetra::Vector.
    using vector_type = ::Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    //! Type of a @ref ::Tpetra::MultiVector instantiated with the same quadruple of template parameters.
    using multivector_type = ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

public:
    //! Default constructor.
    Vector() = default;

    /**
     * @brief Constructor that constructs a @ref NOX::Tpetra::Vector object by encapsulating
     *        the passed shared pointer holding the @ref ::Tpetra::Vector.
     *
     * This constructor is similar to the first constructor in @ref NOX::Thyra::Vector in that
     * no new memory is allocated because the passed shared pointer is encapsulated. As in
     * @ref NOX::Thyra::Vector, this constructor sets @ref do_implicit_weighting to @c false.
     */
    Vector(const Teuchos::RCP<vector_type>& source);

    /**
     * @brief Constructor that constructs a @ref NOX::Tpetra::Vector that wraps a @ref ::Tpetra::Vector
     *        constructed in newly allocated memory by using the map and, if flagged, copying the elements
     *        of the referenced @ref ::Tpetra::Vector.
     *
     * The parameter @p type is a copy-control flag. If <tt> type = NOX::ShapeCopy </tt>, only the
     * map is used. If <tt> type = NOX::DeepCopy </tt>, the elements are copied.
     *
     * This constructor is similar to the second constructor in @ref NOX::Thyra::Vector in that new
     * memory is allocated. As in @ref NOX::Thyra::Vector, it sets @ref do_implicit_weighting to @c false.
     */
    Vector(const vector_type& source, CopyType type = DeepCopy);

    /**
     * @brief Copy constructor. It constructs a @ref NOX::Tpetra::Vector that wraps a @ref ::Tpetra::Vector
     *        constructed in newly allocated memory by using the map and, if flagged, copying the elements
     *        of the @ref ::Tpetra::Vector wrapped by the referenced @ref NOX::Tpetra::Vector.
     *
     * The parameter @p type is a copy-control flag. If <tt> type = NOX::ShapeCopy </tt>, only the
     * map is used. If <tt> type = NOX::DeepCopy </tt>, the elements are copied.
     *
     * This copy constructor is similar to the copy constructor in @ref NOX::Thyra::Vector. As in
     * @ref NOX::Thyra::Vector, the pointer to the @ref ::Tpetra::Vector used for weighting inner products
     * and norms and the flag @ref do_implicit_weighting are copied if this pointer is non-null.
     */
    Vector(const Vector& source, CopyType type = DeepCopy);

    //! Returns the encapsulated shared pointer holding the @ref ::Tpetra::Vector.
    Teuchos::RCP<vector_type> getTpetraVector() const;

    /**
     * @brief Copies the elements of the @ref ::Tpetra::Vector held by the passed shared pointer into the
     *        elements of the wrapped @ref ::Tpetra::Vector. The maps of the two objects must be compatible.
     *
     * As in @ref NOX::Thyra::Vector, the weighting objects are not copied. They are copied only by the
     * copy constructor and the @ref clone function.
     */
    Abstract::Vector& operator=(const Teuchos::RCP<const vector_type>& source);

    /**
     * @brief Copies the elements of the @ref ::Tpetra::Vector wrapped in the referenced
     *        @ref NOX::Tpetra::Vector into the elements of the wrapped @ref ::Tpetra::Vector.
     *        The maps of the two objects must be compatible.
     *
     * See above for the behavior with regards to the weighting objects.
     */
    Abstract::Vector& operator=(const Vector& source);

    //! See above.
    Abstract::Vector& operator=(const Abstract::Vector& source) override;

    //! Initializes every element with @p gamma.
    Abstract::Vector& init(double gamma) override;

    /**
     * @brief Puts the absolute values of the elements of the @ref ::Tpetra::Vector wrapped in the referenced
     *        @ref NOX::Tpetra::Vector into the the elements of the wrapped @ref ::Tpetra::Vector.
     */
    Abstract::Vector& abs(const Vector& source);

    //! See above.
    Abstract::Vector& abs(const Abstract::Vector& source) override;

    /**
     * @brief Initializes each element of the wrapped @ref ::Tpetra::Vector with a random value.
     *
     * As in @ref NOX::Thyra::Vector, the randomly generated values are between @c -1 and @c 1. The
     * parameter @p seed is unused because the underlying @ref ::Tpetra function generates its own seed.
     */
    Abstract::Vector& random(bool useSeed = false, int seed = 1) override;

    /**
     * @brief Puts the inverses of the elements of the @ref ::Tpetra::Vector wrapped in the referenced
     *        @ref NOX::Tpetra::Vector into the the elements of the wrapped @ref ::Tpetra::Vector.
     */
    Abstract::Vector& reciprocal(const Vector& source);

    //! See above.
    Abstract::Vector& reciprocal(const Abstract::Vector& source) override;

    //! Multiplies the elements of the wrapped @ref ::Tpetra::Vector with @p gamma.
    Abstract::Vector& scale(double gamma) override;

    /**
     * @brief Multiplies the elements of the wrapped @ref ::Tpetra::Vector elementwise with the
     *        elements of the @ref ::Tpetra::Vector wrapped in the referenced @ref NOX::Tpetra::Vector.
     */
    Abstract::Vector& scale(const Vector& a);

    //! See above.
    Abstract::Vector& scale(const Abstract::Vector& a) override;

    /**
     * @brief Updates the elements of the wrapped @ref ::Tpetra::Vector with the elements of the
     *        @ref ::Tpetra::Vector wrapped in the referenced @ref NOX::Tpetra::Vector.
     *
     * This operation is such that @c x becomes <tt> alpha * a + gamma * x </tt>, with @c x this vector.
     */
    Abstract::Vector& update(double alpha, const Vector& a, double gamma = 0.);

    //! See above.
    Abstract::Vector& update(double alpha, const Abstract::Vector& a, double gamma = 0.) override;

    /**
     * @brief Updates the elements of the wrapped @ref ::Tpetra::Vector with the elements of the
     *        @ref ::Tpetra::Vector objects wrapped in the referenced @ref NOX::Tpetra::Vector objects.
     *
     * This operation is such that @c x becomes <tt> alpha * a + beta * b + gamma * x </tt>, with @c x this vector.
     */
    Abstract::Vector& update(
        double alpha,
        const Vector& a,
        double beta,
        const Vector& b,
        double gamma = 0.
    );

    //! See above.
    Abstract::Vector& update(
        double alpha,
        const Abstract::Vector& a,
        double beta,
        const Abstract::Vector&b,
        double gamma = 0.
    ) override;

    /**
     * @brief Clone function that constructs a @ref NOX::Tpetra::Vector that wraps a @ref ::Tpetra::Vector
     *        constructed in newly allocated memory by using the map and, if flagged, copying the elements
     *        of the encapsulated @ref ::Tpetra::Vector.
     *
     * The parameter @p type is a copy-control flag. If <tt> type = NOX::ShapeCopy </tt>, only the map is used.
     */
    Teuchos::RCP<Abstract::Vector> clone(CopyType type = DeepCopy) const override;

    /**
     * @brief Create function that constructs a @c NOX::Tpetra::MultiVector that wraps a @ref ::Tpetra::MultiVector
     *        with <tt> numVecs + 1 </tt> vectors constructed in newly allocated memory from the wrapped
     *        @ref ::Tpetra::Vector and the @ref ::Tpetra::Vector objects wrapped in the passed array
     *        of @ref NOX::Tpetra::Vector objects.
     *
     * The parameter @p type is a copy-control flag. If <tt> type = NOX::ShapeCopy </tt>, only the map is used.
     * If <tt> type = NOX::DeepCopy </tt>, the elements are copied: the elements of the wrapped
     * @ref ::Tpetra::Vector are copied to the first column, with the elements given by @p vecs copied
     * to the remaining columns.
     *
     * As in @ref NOX::Thyra::Vector, weighting objects are *not* copied into the multivector.
     */
    Teuchos::RCP<Abstract::MultiVector> createMultiVector(
        const Abstract::Vector* const* vecs,
        int numVecs,
        CopyType type = DeepCopy
    ) const override;

    /**
     * @brief Create function that constructs a @ref NOX::Tpetra::MultiVector that wraps a @ref ::Tpetra::MultiVector
     *        with <tt> numVecs </tt> vectors constructed in newly allocated memory from the wrapped @ref ::Tpetra::Vector.
     *
     * The parameter @p type is a copy-control flag. If <tt> type = NOX::ShapeCopy </tt>, only the map is used.
     * If <tt> type = NOX::DeepCopy </tt>, the elements are copied: the elements of the wrapped
     * @ref ::Tpetra::Vector are copied to each column of the created multivector.
     *
     * As in @ref NOX::Thyra::Vector, weighting objects (if not null) *are* copied into the multivector.
     */
    Teuchos::RCP<Abstract::MultiVector> createMultiVector(int numVecs, CopyType type = DeepCopy) const override;

    /**
     * @brief Returns the norm of the wrapped @ref ::Tpetra::Vector.
     *
     * The parameter @p type determines the type of norm. The one-norm is returned for <tt> type = NOX::OneNorm </tt>,
     * the two-norm for <tt> type = NOX::TwoNorm </tt>, and the infinity-norm for <tt> type = NOX::MaxNorm </tt>.
     *
     * As in @ref NOX::Thyra::Vector, the flag @ref do_implicit_weighting controls whether the wrapped @ref ::Tpetra::Vector
     * for weighting inner products and norms (if not null) is used.
     */
    double norm(Abstract::Vector::NormType type = TwoNorm) const override;

    /**
     * @brief Returns the norm of the wrapped @ref ::Tpetra::Vector.
     *
     * As in @ref NOX::Thyra::Vector, the flag @ref do_implicit_weighting controls whether, in addition to the weights
     * in @p weights, the wrapped @ref ::Tpetra::Vector for weighting inner products and norms (if not null) is used.
     */
    double norm(const Vector& weights) const;

    //! See above.
    double norm(const Abstract::Vector& weights) const override;

    /**
     * @brief Returns the inner product of the wrapped @ref ::Tpetra::Vector vector with
     *        the @ref ::Tpetra::Vector wrapped in the referenced @ref NOX:Tpetra::Vector.
     *
     * As in @ref NOX::Thyra::Vector, the flag @ref do_implicit_weighting controls whether the @ref ::Tpetra::Vector
     * for weighting inner products and norms (if not null) is used.
     */
    double innerProduct(const Vector& y) const;

    //! See above.
    double innerProduct(const Abstract::Vector& y) const override;

    /**
     * @brief Sets the @ref ::Tpetra::Vector used for weighting inner products and norms.
     *
     * Setting a vector of weights results in norms being computed as
     *     - \f$\|x\|_{1} = \sum |w_{i} x_{i})|\f$,
     *     - \f$\|x\|_{2} = \sqrt{\sum (w_{i} x_{i})^{2}}\f$,
     *     - \f$\|x\|_{\infty} = \max |w_{i} x_{i})|\f$,
     * and inner products as
     *     - \f$(x, y) = \sum w_{i} x_{i}) w_{i} y(i)\f$.
     */
    void setWeightVector(const Teuchos::RCP<const vector_type>& weightVec_);

    /**
     * @brief Returns flag to indicate whether @ref setWeightVector has been called
     *        to set a @ref ::Tpetra::Vector used for weighting inner products and norms.
     */
    bool hasWeightVector() const;

    //! Returns the shared pointer holding the @ref ::Tpetra::Vector used for weighting inner products and norms.
    Teuchos::RCP<const vector_type> getWeightVector() const;

    //! Returns the number of elements in the wrapped @ref ::Tpetra::Vector.
    NOX::size_type length() const override;

    /**
     * @brief Returns the flag that controls whether the wrapped @ref ::Tpetra::Vector for weighting
     *        inner products and norms (if not null) is used.
     */
    bool getImplicitWeighting() const override;

    /**
     * @brief Sets the flag that controls whether the wrapped @ref ::Tpetra::Vector for weighting
     *        inner products and norms (if not null) is used.
     */
    void setImplicitWeighting(bool do_implicit_weighting_) override;

private:
    //! Shared pointer holding the @ref ::Tpetra::Vector wrapped in this object.
    Teuchos::RCP<vector_type> tpetraVec;

    //! Shared pointer holding the @ref ::Tpetra::Vector used for weighting inner products and norms.
    Teuchos::RCP<const vector_type> weightVec;

    /**
     * @brief Flag that controls whether the @ref ::Tpetra::Vector for weighting
     *        inner products and norms (if not null) is used.
     */
    bool do_implicit_weighting = false;
};

} // namespace NOX::Tpetra

#endif // NOX_TPETRA_VECTOR_HPP
