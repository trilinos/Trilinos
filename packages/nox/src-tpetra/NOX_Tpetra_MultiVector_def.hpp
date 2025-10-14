// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_TPETRA_MULTIVECTOR_DEF_HPP
#define NOX_TPETRA_MULTIVECTOR_DEF_HPP

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"

namespace NOX::Tpetra
{

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
MultiVector(const Teuchos::RCP<multivector_type>& source)
  : tpetraMultiVec(source),
    noxTpetraVecs(source->getNumVectors()) {}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
MultiVector(const multivector_type& source, CopyType type)
{
    switch (type)
    {
        case DeepCopy:
            tpetraMultiVec = Teuchos::make_rcp<multivector_type>(source, Teuchos::Copy);
            break;
        case ShapeCopy:
            tpetraMultiVec = Teuchos::make_rcp<multivector_type>(source.getMap(), source.getNumVectors());
            break;
    }

    noxTpetraVecs = std::vector<std::optional<nox_tpetra_vector_type>>(source.getNumVectors());
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
MultiVector(const MultiVector& source, CopyType type)
{
    switch (type)
    {
        case DeepCopy:
            tpetraMultiVec = Teuchos::make_rcp<multivector_type>(*source.getTpetraMultiVector(), Teuchos::Copy);
            break;
        case ShapeCopy:
            tpetraMultiVec = Teuchos::make_rcp<multivector_type>(
                source.getTpetraMultiVector()->getMap(),
                source.getTpetraMultiVector()->getNumVectors()
            );
            break;
    }
    if (!source.getWeightVector().is_null())
    {
        weightVec = source.getWeightVector();
        do_implicit_weighting = source.do_implicit_weighting;
    }

    noxTpetraVecs = std::vector<std::optional<nox_tpetra_vector_type>>(source.numVectors());
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::multivector_type>
MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getTpetraMultiVector() const {
    return tpetraMultiVec;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::MultiVector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator=(const Teuchos::RCP<const multivector_type>& source)
{
    tpetraMultiVec->assign(*source);
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::MultiVector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator=(const MultiVector& source)
{
    tpetraMultiVec->assign(*source.getTpetraMultiVector());
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::MultiVector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator=(const Abstract::MultiVector& source)
{
    return operator=(dynamic_cast<const MultiVector&>(source));
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::MultiVector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
init(double gamma)
{
    tpetraMultiVec->putScalar(gamma);
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::MultiVector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
random(bool useSeed, int seed)
{
    if (!tpetraMultiVec->isDistributed())
#if KOKKOS_VERSION >= 40799
    {
        if (tpetraMultiVec->getMap()->getComm()->getRank() == 0) {
            tpetraMultiVec->randomize(-KokkosKernels::ArithTraits<Scalar>::one(), KokkosKernels::ArithTraits<Scalar>::one());
        }
        else {
            tpetraMultiVec->putScalar(KokkosKernels::ArithTraits<Scalar>::zero());
        }
        tpetraMultiVec->reduce();
    }
    else {
        tpetraMultiVec->randomize(-KokkosKernels::ArithTraits<Scalar>::one(), KokkosKernels::ArithTraits<Scalar>::one());
    }
#else
    {
        if (tpetraMultiVec->getMap()->getComm()->getRank() == 0) {
            tpetraMultiVec->randomize(-Kokkos::ArithTraits<Scalar>::one(), Kokkos::ArithTraits<Scalar>::one());
        }
        else {
            tpetraMultiVec->putScalar(Kokkos::ArithTraits<Scalar>::zero());
        }
        tpetraMultiVec->reduce();
    }
    else {
        tpetraMultiVec->randomize(-Kokkos::ArithTraits<Scalar>::one(), Kokkos::ArithTraits<Scalar>::one());
    }
#endif
    return *this;
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::MultiVector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
scale(double gamma)
{
    tpetraMultiVec->scale(gamma);
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::MultiVector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
update(double alpha, const MultiVector& a, double gamma)
{
    tpetraMultiVec->update(alpha, *a.getTpetraMultiVector(), gamma);
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::MultiVector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
update(double alpha, const Abstract::MultiVector& a, double gamma)
{
    return update(alpha, dynamic_cast<const MultiVector&>(a), gamma);
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::MultiVector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
update(
    double alpha,
    const MultiVector& a,
    double beta,
    const MultiVector& b,
    double gamma
)
{
    tpetraMultiVec->update(alpha, *a.getTpetraMultiVector(), beta, *b.getTpetraMultiVector(), gamma);
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::MultiVector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
update(
    double alpha,
    const Abstract::MultiVector& a,
    double beta,
    const Abstract::MultiVector& b,
    double gamma
)
{
    return update(
        alpha,
        dynamic_cast<const MultiVector&>(a),
        beta,
        dynamic_cast<const MultiVector&>(b),
        gamma
    );
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::MultiVector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
update(
    Teuchos::ETransp transb,
    double alpha,
    const MultiVector& a,
    const Abstract::MultiVector::DenseMatrix& b,
    double gamma
)
{
    using map_type = ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;

    const auto localMap = Teuchos::make_rcp<map_type>(b.numRows(), 0, tpetraMultiVec->getMap()->getComm(), ::Tpetra::LocallyReplicated);

    const auto B = Teuchos::make_rcp<multivector_type>(localMap, b.numCols());

    {
        const auto local_B = B->getLocalViewHost(::Tpetra::Access::OverwriteAll);
        // The dense matrix is assumed to be small, so that copying the values in serial loops is acceptable.
        for (int irow = 0; irow < b.numRows(); ++irow)
        {
            for (int icol = 0; icol < b.numCols(); ++icol) {
                local_B(irow, icol) = static_cast<typename multivector_type::impl_scalar_type>(b(irow, icol));
            }
        }
    }

    tpetraMultiVec->multiply(Teuchos::NO_TRANS, transb, alpha, *a.getTpetraMultiVector(), *B, gamma);

    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::MultiVector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
update(
    Teuchos::ETransp transb,
    double alpha,
    const Abstract::MultiVector& a,
    const Abstract::MultiVector::DenseMatrix& b,
    double gamma
)
{
    return update(
        transb,
        alpha,
        dynamic_cast<const MultiVector&>(a),
        b,
        gamma
    );
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<Abstract::MultiVector> MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
clone(CopyType type) const
{
    return Teuchos::make_rcp<MultiVector>(*this, type);
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<Abstract::MultiVector> MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
clone(int numVecs) const
{
    auto newTpetraMultiVec = Teuchos::make_rcp<multivector_type>(tpetraMultiVec->getMap(), numVecs);
    return Teuchos::make_rcp<MultiVector>(std::move(newTpetraMultiVec));
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
norm(std::vector<double>& result, Abstract::Vector::NormType type) const
{
    using double_view_um_t = Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    if (weightVec.is_null() || !do_implicit_weighting)
    {
        switch (type)
        {
            case Abstract::Vector::OneNorm:
                tpetraMultiVec->norm1(double_view_um_t(result.data(), result.size()));
                break;
            case Abstract::Vector::MaxNorm:
                tpetraMultiVec->normInf(double_view_um_t(result.data(), result.size()));
                break;
            case Abstract::Vector::TwoNorm:
            default:
                tpetraMultiVec->norm2(double_view_um_t(result.data(), result.size()));
                break;
        }
    }
    else
    {
        const auto tmpMultiVec = Teuchos::make_rcp<multivector_type>(
            tpetraMultiVec->getMap(), tpetraMultiVec->getNumVectors()
        );

        tmpMultiVec->elementWiseMultiply(1., *weightVec, *tpetraMultiVec, 0.);

        switch (type)
        {
            case Abstract::Vector::OneNorm:
                tmpMultiVec->norm1(double_view_um_t(result.data(), result.size()));
                break;
            case Abstract::Vector::MaxNorm:
                tmpMultiVec->normInf(double_view_um_t(result.data(), result.size()));
                break;
            case Abstract::Vector::TwoNorm:
            default:
                tmpMultiVec->norm2(double_view_um_t(result.data(), result.size()));
                break;
        }
    }
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::MultiVector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
setBlock(const MultiVector& source, const std::vector<int>& indices)
{
    for (size_t idx = 0; idx < indices.size(); ++idx) {
        tpetraMultiVec->getVectorNonConst(indices[idx])->assign(*(source.getTpetraMultiVector()->getVector(idx)));
    }
    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::MultiVector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
setBlock(const Abstract::MultiVector& source, const std::vector<int>& indices)
{
    return setBlock(dynamic_cast<const MultiVector&>(source), indices);
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::MultiVector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
augment(const MultiVector& source)
{
    auto newTpetraMultiVec = Teuchos::make_rcp<multivector_type>(
        tpetraMultiVec->getMap(), tpetraMultiVec->getNumVectors() + source.numVectors()
    );

    for (size_t idx = 0; idx < tpetraMultiVec->getNumVectors(); ++idx) {
        newTpetraMultiVec->getVectorNonConst(idx)->assign(*(tpetraMultiVec->getVector(idx)));
    }

    for (size_t idx = 0; idx < static_cast<size_t>(source.numVectors()); ++idx) {
        newTpetraMultiVec->getVectorNonConst(tpetraMultiVec->getNumVectors() + idx)->assign(*(source.getTpetraMultiVector()->getVector(idx)));
    }

    tpetraMultiVec = std::move(newTpetraMultiVec);

    noxTpetraVecs = std::vector<std::optional<nox_tpetra_vector_type>>(tpetraMultiVec->getNumVectors() + source.numVectors());

    return *this;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::MultiVector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
augment(const Abstract::MultiVector& source)
{
    return augment(dynamic_cast<const MultiVector&>(source));
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Abstract::Vector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator[](int i)
{
    if (!noxTpetraVecs[i].has_value()) {
        noxTpetraVecs[i].emplace(tpetraMultiVec->getVectorNonConst(i));
    }
    return noxTpetraVecs[i].value();
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
const Abstract::Vector& MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator[](int i) const
{
    if (!noxTpetraVecs[i].has_value()) {
        noxTpetraVecs[i].emplace(tpetraMultiVec->getVectorNonConst(i));
    }
    return noxTpetraVecs[i].value();
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<Abstract::MultiVector> MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
subCopy(const std::vector<int>& indices) const
{
    std::vector<size_t> indices_transformed(indices.size());
    std::transform(
        indices.cbegin(), indices.cend(),
        indices_transformed.begin(),
        [](int idx) { return static_cast<size_t>(idx); }
    );

    auto newTpetraMultiVec = tpetraMultiVec->subCopy(indices_transformed);

    return Teuchos::make_rcp<MultiVector>(std::move(newTpetraMultiVec));
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<Abstract::MultiVector> MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
subView(const std::vector<int>& indices) const
{
    std::vector<size_t> indices_transformed(indices.size());
    std::transform(
        indices.cbegin(), indices.cend(),
        indices_transformed.begin(),
        [](int idx) { return static_cast<size_t>(idx); }
    );

    auto newTpetraMultiVec = tpetraMultiVec->subViewNonConst(indices_transformed);

    return Teuchos::make_rcp<MultiVector>(std::move(newTpetraMultiVec));
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
multiply(double alpha, const MultiVector& y, Abstract::MultiVector::DenseMatrix& b) const
{
    for (int irow = 0; irow < y.numVectors(); ++irow)
    {
        for (int icol = 0; icol < static_cast<int>(tpetraMultiVec->getNumVectors()); ++icol) {
            b(irow, icol) = alpha * y.getTpetraMultiVector()->getVector(irow)->dot(*(tpetraMultiVec->getVector(icol)));
        }
    }
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
multiply(double alpha, const Abstract::MultiVector& y, Abstract::MultiVector::DenseMatrix& b) const
{
    multiply(alpha, dynamic_cast<const MultiVector&>(y), b);
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
int MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
numVectors() const
{
    return tpetraMultiVec->getNumVectors();
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
print(std::ostream& stream) const
{
    tpetraMultiVec->describe(
        *Teuchos::getFancyOStream(Teuchos::rcpFromRef(stream)),
        Teuchos::VERB_EXTREME
    );
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
setWeightVector(const Teuchos::RCP<const vector_type>& weightVec_)
{
    weightVec = weightVec_;
    do_implicit_weighting = true;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
bool MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
hasWeightVector() const {
    return !weightVec.is_null();
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::vector_type>
MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getWeightVector() const {
    return weightVec;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
NOX::size_type MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
length() const {
    return tpetraMultiVec->getGlobalLength();
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
bool MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getImplicitWeighting() const {
    return do_implicit_weighting;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
setImplicitWeighting(bool do_implicit_weighting_)
{
    do_implicit_weighting = do_implicit_weighting_;
}

} // namespace NOX::Tpetra

//! Explicit template instantiation macro for @ref NOX::Tpetra::MultiVector.
#define NOX_TPETRA_MULTIVECTOR_INSTANT(S, L, G, N) template class NOX::Tpetra::MultiVector<S, L, G, N>;

#endif // NOX_TPETRA_MULTIVECTOR_DEF_HPP
