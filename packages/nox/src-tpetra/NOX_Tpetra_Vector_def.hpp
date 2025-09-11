// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_TPETRA_VECTOR_DEF_HPP
#define NOX_TPETRA_VECTOR_DEF_HPP

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"

namespace NOX::Tpetra
{

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
Vector(const Teuchos::RCP<vector_type>& source)
  : tpetraVec(source), do_implicit_weighting(false) {}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
Vector(const vector_type& source, CopyType type)
  : do_implicit_weighting(false)
{
    switch (type)
    {
        case DeepCopy:
            tpetraVec = Teuchos::make_rcp<vector_type>(source, Teuchos::Copy);
            break;
        case ShapeCopy:
            tpetraVec = Teuchos::make_rcp<vector_type>(source.getMap());
            break;
    }
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
Vector(const Vector& source, CopyType type)
  : do_implicit_weighting(false)
{
    switch (type)
    {
        case DeepCopy:
            tpetraVec = Teuchos::make_rcp<vector_type>(*source.getTpetraVector(), Teuchos::Copy);
            break;
        case ShapeCopy:
            tpetraVec = Teuchos::make_rcp<vector_type>(source.getTpetraVector()->getMap());
            break;
    }
    if (!source.getWeightVector().is_null())
    {
        weightVec = source.getWeightVector();
        do_implicit_weighting = source.do_implicit_weighting;
    }
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<typename Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::vector_type>
Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getTpetraVector() const {
    return tpetraVec;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
init(double gamma) {
    tpetraVec->putScalar(gamma);
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator=(const Teuchos::RCP<const vector_type>& source)
{
    tpetraVec->assign(*source);
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator=(const Vector& source)
{
    tpetraVec->assign(*source.getTpetraVector());
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator=(const Abstract::Vector& source) {
    return operator=(dynamic_cast<const Vector&>(source));
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
abs(const Vector& source)
{
    tpetraVec->abs(*source.getTpetraVector());
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
abs(const Abstract::Vector& source) {
    return abs(dynamic_cast<const Vector&>(source));
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
random(bool useSeed, int seed)
{
    if (!tpetraVec->isDistributed())
    {
        if (tpetraVec->getMap()->getComm()->getRank() == 0) {
            tpetraVec->randomize(-Teuchos::ScalarTraits<Scalar>::one(), Teuchos::ScalarTraits<Scalar>::one());
        }
        else {
            tpetraVec->putScalar(Teuchos::ScalarTraits<Scalar>::zero());
        }
        tpetraVec->reduce();
    }
    else {
        tpetraVec->randomize(-Teuchos::ScalarTraits<Scalar>::one(), Teuchos::ScalarTraits<Scalar>::one());
    }
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
reciprocal(const Vector& source)
{
    tpetraVec->reciprocal(*source.getTpetraVector());
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
reciprocal(const Abstract::Vector& source) {
    return reciprocal(dynamic_cast<const Vector&>(source));
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
scale(double gamma)
{
    tpetraVec->scale(gamma);
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
scale(const Vector& a)
{
    tpetraVec->elementWiseMultiply(1., *a.getTpetraVector(), *tpetraVec, 0.);
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
scale(const Abstract::Vector& a) {
    return scale(dynamic_cast<const Vector&>(a));
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
update(double alpha, const Vector& a, double gamma)
{
    tpetraVec->update(alpha, *a.getTpetraVector(), gamma);
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
update(double alpha, const Abstract::Vector& a, double gamma) {
        return update(alpha, dynamic_cast<const Vector&>(a), gamma);
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
update(double alpha, const Vector& a, double beta, const Vector& b, double gamma)
{
    tpetraVec->update(alpha, *a.getTpetraVector(), beta, *b.getTpetraVector(), gamma);
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
update(double alpha, const Abstract::Vector& a, double beta, const Abstract::Vector& b, double gamma) {
    return update(alpha, dynamic_cast<const Vector&>(a), beta, dynamic_cast<const Vector&>(b), gamma);
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<Abstract::Vector> Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
clone(CopyType type) const {
    return Teuchos::make_rcp<Vector>(*this, type);
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<Abstract::MultiVector> Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
createMultiVector(const Abstract::Vector *const *vecs, int numVecs, CopyType type) const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented yet.");
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<Abstract::MultiVector> Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
createMultiVector(int numVecs, CopyType type) const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented yet.");
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
double Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
norm(Abstract::Vector::NormType type) const
{
    double norm;

    if (weightVec.is_null() || !do_implicit_weighting)
    {
        switch (type)
        {
            case OneNorm:
                norm = tpetraVec->norm1();
                break;
            case MaxNorm:
                norm = tpetraVec->normInf();
                break;
            case TwoNorm:
            default:
                norm = tpetraVec->norm2();
                break;
        }
    }
    else
    {
        const auto tmpVec = Teuchos::make_rcp<vector_type>(tpetraVec->getMap());

        tmpVec->elementWiseMultiply(1., *weightVec, *tpetraVec, 0.);

        switch (type)
        {
            case OneNorm:
                norm = tmpVec->norm1();
                break;
            case MaxNorm:
                norm = tmpVec->normInf();
                break;
            case TwoNorm:
            default:
                norm = tmpVec->norm2();
                break;
        }
    }

    return norm;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
double Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
norm(const Vector& weights) const
{
    const auto tmpVec = Teuchos::make_rcp<vector_type>(*tpetraVec, Teuchos::Copy);

    if (!weightVec.is_null() && do_implicit_weighting) {
        tmpVec->elementWiseMultiply(1., *weightVec, *tpetraVec, 0.);
    }

    tmpVec->elementWiseMultiply(1., *weights.getTpetraVector(), *tpetraVec, 0.);

    return std::sqrt(tmpVec->dot(*tpetraVec));
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
double Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
norm(const Abstract::Vector& weights) const {
    return norm(dynamic_cast<const Vector&>(weights));
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
double Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
innerProduct(const Vector& y) const
{
    double dot;

    if (weightVec.is_null() || !do_implicit_weighting) {
        dot = tpetraVec->dot(*y.getTpetraVector());
    }
    else
    {
        const auto tmpVec = Teuchos::make_rcp<vector_type>(tpetraVec->getMap());

        tmpVec->elementWiseMultiply(1., *weightVec, *tpetraVec, 0.);
        tmpVec->elementWiseMultiply(1., *weightVec, *tmpVec, 0.);

        dot = tmpVec->dot(*y.getTpetraVector());
    }

    return dot;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
double Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
innerProduct(const Abstract::Vector& y) const {
    return innerProduct(dynamic_cast<const Vector&>(y));
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
setWeightVector(const Teuchos::RCP<const vector_type>& weightVec_)
{
    weightVec = weightVec_;
    do_implicit_weighting = true;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
bool Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
hasWeightVector() const {
    return !weightVec.is_null();
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const typename Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::vector_type>
Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getWeightVector() const {
    return weightVec;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
NOX::size_type Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
length() const {
    return tpetraVec->getGlobalLength();
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
bool Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getImplicitWeighting() const {
    return do_implicit_weighting;

}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
setImplicitWeighting(bool do_implicit_weighting_) {
    do_implicit_weighting = do_implicit_weighting_;
}

} // namespace NOX::Tpetra

//! Explicit template instantiation macro for @ref NOX::Tpetra::Vector.
#define NOX_TPETRA_VECTOR_INSTANT(S, L, G, N) template class NOX::Tpetra::Vector<S, L, G, N>;

#endif // NOX_TPETRA_VECTOR_DEF_HPP
