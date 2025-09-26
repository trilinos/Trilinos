// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"

#include "NOX.H"
#include "NOX_Tpetra_MultiVector.hpp"
#include "NOX_Tpetra_Vector.hpp"

#include "Teuchos_UnitTestHarness.hpp"

template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
using tpetra_map_t = ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
using tpetra_vector_t = ::Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
using tpetra_multivector_t = ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
using nox_tpetra_vector_t = NOX::Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
using nox_tpetra_multivector_t = NOX::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

template <typename Scalar>
using magnitude_t = typename Kokkos::ArithTraits<Scalar>::magnitudeType;

//! Function to check solution.
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
bool checkMultiVectors(
    const Teuchos::RCP<tpetra_multivector_t<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& a,
    const Teuchos::RCP<tpetra_multivector_t<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& b,
    const Kokkos::View<magnitude_t<Scalar>*, Kokkos::HostSpace>& expt,
    Teuchos::FancyOStream& out,
    bool success
)
{
    constexpr auto tol = 10 * Kokkos::Experimental::epsilon<magnitude_t<Scalar>>::value;

    TEUCHOS_TEST_EQUALITY(a->getNumVectors(), b->getNumVectors(), out, success);

    b->update(-1 * Kokkos::ArithTraits<Scalar>::one(), *a, Kokkos::ArithTraits<Scalar>::one());
    const Kokkos::View<magnitude_t<Scalar>*, Kokkos::HostSpace> norms_b("norms b after update", b->getNumVectors());
    b->norm2(norms_b);
    for (size_t idx = 0; idx < norms_b.extent(0); ++idx) {
        TEUCHOS_TEST_COMPARE(norms_b(idx), <, tol, out, success);
    }

    const Kokkos::View<magnitude_t<Scalar>*, Kokkos::HostSpace> norms_a("norms a", a->getNumVectors());
    a->norm2(norms_a);
    for (size_t idx = 0; idx < norms_a.extent(0); ++idx) {
        TEUCHOS_TEST_FLOATING_EQUALITY(norms_a(idx), expt(idx), tol, out, success);
    }

    return success;
}

constexpr size_t numLocalElements = 1000;
constexpr size_t numCols = 10;

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, Traits, S, LO, GO, N)
{
    static_assert( ! std::is_default_constructible_v<nox_tpetra_multivector_t<S, LO, GO, N>>);
    static_assert(   std::is_destructible_v<         nox_tpetra_multivector_t<S, LO, GO, N>>);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, ConstructorFromRCP, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);

    // Construct NOX multivector that wraps the Tpetra multivector and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox(x);

    TEST_EQUALITY(x_nox.getTpetraMultiVector(), x);
    TEST_EQUALITY(x.strong_count(), 2);

    TEST_ASSERT( ! x_nox.getImplicitWeighting());

    TEST_EQUALITY(x_nox.length(), static_cast<NOX::size_type>(numGlobalElements));
    TEST_EQUALITY(x_nox.numVectors(), static_cast<int>(numCols));
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, ConstructorFromRef, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);

    x->putScalar(Kokkos::ArithTraits<S>::one());

    // Construct NOX multivector with copy of Tpetra multivector using the flag NOX::DeepCopy and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox_deep_copy(*x, NOX::DeepCopy);

    TEUCHOS_TEST_INEQUALITY(x_nox_deep_copy.getTpetraMultiVector(), x, out, success);
    TEUCHOS_TEST_EQUALITY(x.strong_count(), 1, out, success);

    TEUCHOS_TEST_ASSERT( ! x_nox_deep_copy.getImplicitWeighting(), out, success);

    const Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    Kokkos::deep_copy(expt, Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements));
    success = checkMultiVectors(x_nox_deep_copy.getTpetraMultiVector(), x, expt, out, success);

    // Construct NOX vector with copy of Tpetra vector using the flag NOX::ShapeCopy and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox_shape_copy(*x, NOX::ShapeCopy);

    TEUCHOS_TEST_EQUALITY(x_nox_shape_copy.getTpetraMultiVector()->getMap(), map, out, success);

    TEUCHOS_TEST_ASSERT( ! x_nox_shape_copy.getImplicitWeighting(), out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, CopyConstructor, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);

    x->putScalar(Kokkos::ArithTraits<S>::one());

    // Construct NOX vector that wraps the Tpetra multivector.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox(x);

    TEUCHOS_TEST_EQUALITY(x.strong_count(), 2, out, success);

    // Copy construct NOX multivector using the flag NOX::DeepCopy and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox_deep_copy(x_nox, NOX::DeepCopy);

    TEUCHOS_TEST_INEQUALITY(x_nox_deep_copy.getTpetraMultiVector(), x, out, success);
    TEUCHOS_TEST_EQUALITY(x.strong_count(), 2, out, success);

    TEUCHOS_TEST_ASSERT( ! x_nox_deep_copy.getImplicitWeighting(), out, success);

    const Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    Kokkos::deep_copy(expt, Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements));
    success = checkMultiVectors(x_nox_deep_copy.getTpetraMultiVector(), x, expt, out, success);

    // Copy construct NOX multivector using the flag NOX::ShapeCopy and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox_shape_copy(x_nox, NOX::ShapeCopy);

    TEUCHOS_TEST_ASSERT( ! x_nox_shape_copy.getImplicitWeighting(), out, success);

    TEUCHOS_TEST_EQUALITY(x_nox_shape_copy.getTpetraMultiVector()->getMap(), map, out, success);

    // After setting a weight vector, copy construct NOX vector using the flag NOX::DeepCopy and check.
    const auto w = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    nox_tpetra_multivector_t<S, LO, GO, N> x_nox_weighted(x_nox, NOX::DeepCopy);
    x_nox_weighted.setWeightVector(w);
    TEUCHOS_TEST_ASSERT(x_nox_weighted.hasWeightVector(), out, success);
    TEUCHOS_TEST_EQUALITY(x_nox_weighted.getWeightVector(), w, out, success);
    TEUCHOS_TEST_ASSERT(x_nox_weighted.getImplicitWeighting(), out, success);

    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox_weighted_deep_copy(x_nox_weighted, NOX::DeepCopy);
    TEUCHOS_TEST_ASSERT(x_nox_weighted_deep_copy.hasWeightVector(), out, success);
    TEUCHOS_TEST_EQUALITY(x_nox_weighted_deep_copy.getWeightVector(), w, out, success);
    TEUCHOS_TEST_ASSERT(x_nox_weighted_deep_copy.getImplicitWeighting(), out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, CopyAssignment, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);
    const auto y = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);

    x->putScalar(Kokkos::ArithTraits<S>::one());
    y->putScalar(Kokkos::ArithTraits<S>::zero());

    // Copy x into y through NOX interface.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox(x);
          nox_tpetra_multivector_t<S, LO, GO, N> y_nox(y);

    y_nox = x_nox;

    // Check for correct answer.
    const Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    Kokkos::deep_copy(expt, Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements));
    success = checkMultiVectors(y_nox.getTpetraMultiVector(), x_nox.getTpetraMultiVector(), expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, Init, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);
    const auto y = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);

    // Perform operation with Tpetra directly.
    x->putScalar(2 * Kokkos::ArithTraits<S>::one());

    // Perform this operation through NOX interface.
    nox_tpetra_multivector_t<S, LO, GO, N> y_nox(y);
    y_nox.init(2 * Kokkos::ArithTraits<S>::one());

    // Check for correct answer.
    const Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    Kokkos::deep_copy(expt, Kokkos::ArithTraits<magnitude_t<S>>::squareroot(4 * numGlobalElements));
    success = checkMultiVectors(x, y, expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, Random, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);

    // Randomize through NOX interface.
    nox_tpetra_multivector_t<S, LO, GO, N> x_nox(x);
    x_nox.random();

    // Check that the randomly generated entries are within the expected range.
    std::vector<double> norms(numCols);
    x_nox.norm(norms, NOX::Abstract::Vector::NormType::MaxNorm);
    for (size_t icol = 0; icol < numCols; ++icol) {
        TEST_ASSERT(norms[icol] <= 1.0);
    }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, Scale, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);
    const auto y = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);

    x->putScalar(Kokkos::ArithTraits<S>::one());
    y->putScalar(Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    x->scale(2 * Kokkos::ArithTraits<S>::one());

    // Perform this operation through NOX interface and check.
    nox_tpetra_multivector_t<S, LO, GO, N> y_nox(y);
    y_nox.scale(2 * Kokkos::ArithTraits<S>::one());

    const Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    Kokkos::deep_copy(expt, Kokkos::ArithTraits<magnitude_t<S>>::squareroot(4 * numGlobalElements));
    success = checkMultiVectors(x, y, expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, Update_1, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);
    const auto y = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);
    const auto z = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);

    x->putScalar(Kokkos::ArithTraits<S>::one());
    y->putScalar(Kokkos::ArithTraits<S>::one());
    z->putScalar(Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    x->update(2 * Kokkos::ArithTraits<S>::one(), *z, Kokkos::ArithTraits<S>::one());

    // Perform this operation through NOX interface and check.
          nox_tpetra_multivector_t<S, LO, GO, N> y_nox(y);
    const nox_tpetra_multivector_t<S, LO, GO, N> z_nox(z);

    y_nox.update(2 * Kokkos::ArithTraits<S>::one(), z_nox, Kokkos::ArithTraits<S>::one());

    const Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    Kokkos::deep_copy(expt, Kokkos::ArithTraits<magnitude_t<S>>::squareroot(9 * numGlobalElements));
    success = checkMultiVectors(x, y, expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, Update_2, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto w = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);
    const auto y = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);
    const auto z = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);

    w->putScalar(Kokkos::ArithTraits<S>::one());
    x->putScalar(Kokkos::ArithTraits<S>::one());
    y->putScalar(Kokkos::ArithTraits<S>::one());
    z->putScalar(Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    x->update(2 * Kokkos::ArithTraits<S>::one(), *w, 2 * Kokkos::ArithTraits<S>::one(), *z, Kokkos::ArithTraits<S>::one());

    // Perform this operation through NOX interface and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> w_nox(w);
          nox_tpetra_multivector_t<S, LO, GO, N> y_nox(y);
    const nox_tpetra_multivector_t<S, LO, GO, N> z_nox(z);

    y_nox.update(2 * Kokkos::ArithTraits<S>::one(), w_nox, 2 * Kokkos::ArithTraits<S>::one(), z_nox, Kokkos::ArithTraits<S>::one());

    const Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    Kokkos::deep_copy(expt, Kokkos::ArithTraits<magnitude_t<S>>::squareroot(25 * numGlobalElements));
    success = checkMultiVectors(x, y, expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, Update_3, S, LO, GO, N)
{
    constexpr size_t otherNumCols = 5;

    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map,      numCols);
    const auto y = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map,      numCols);
    const auto z = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, otherNumCols);

    x->putScalar(Kokkos::ArithTraits<S>::one());
    y->putScalar(Kokkos::ArithTraits<S>::one());
    z->putScalar(Kokkos::ArithTraits<S>::one());

    // Create matrix for the update.
    NOX::Abstract::MultiVector::DenseMatrix mat(otherNumCols, numCols);
    for (int icol = 0; icol < static_cast<int>(numCols); ++icol)
    {
        for (int irow = 0; irow < static_cast<int>(otherNumCols); ++irow) {
            mat(irow, icol) = icol;
        }
    }

    // Perform operation with Tpetra directly.
    const auto localMap = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(otherNumCols, 0, comm, ::Tpetra::LocallyReplicated);
    const Teuchos::ArrayView<double> mat_as_array_view(mat.values(), otherNumCols * numCols);
    const auto mat_as_mv = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(localMap, mat_as_array_view, otherNumCols, numCols);

    x->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, Kokkos::ArithTraits<S>::one(), *z, *mat_as_mv, Kokkos::ArithTraits<S>::one());

    // Perform this operation through NOX interface and check.
          nox_tpetra_multivector_t<S, LO, GO, N> y_nox(y);
    const nox_tpetra_multivector_t<S, LO, GO, N> z_nox(z);

    y_nox.update(Teuchos::NO_TRANS, Kokkos::ArithTraits<S>::one(), z_nox, mat, Kokkos::ArithTraits<S>::one());

    const Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    for (size_t idx = 0; idx < numCols; ++idx) {
        expt(idx) = (1 + 5 * idx) * Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements);
    }
    success = checkMultiVectors(x, y, expt, out, success);

    // Perform this operation over again, this time using the transpose, and check.
    NOX::Abstract::MultiVector::DenseMatrix mat_transpose(mat, Teuchos::TRANS);

    y->putScalar(Kokkos::ArithTraits<S>::one());
    y_nox.update(Teuchos::TRANS, Kokkos::ArithTraits<S>::one(), z_nox, mat_transpose, Kokkos::ArithTraits<S>::one());
    success = checkMultiVectors(x, y, expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, OneNorm, S, LO, GO, N)
{
    using manitude_view_um_t = Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    constexpr auto tol = 10 * Kokkos::Experimental::epsilon<magnitude_t<S>>::value;

    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);

    x->putScalar(Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    std::vector<magnitude_t<S>> expt(numCols);
    x->norm1(manitude_view_um_t(expt.data(), numCols));

    // Perform this operation through NOX interface and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox(x);
    std::vector<double> norms(numCols);
    x_nox.norm(norms, NOX::Abstract::Vector::OneNorm);

    TEST_ABSOLUTE_COMPARE_FLOATING_ARRAYS(norms, expt, tol);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, TwoNorm, S, LO, GO, N)
{
    using manitude_view_um_t = Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    constexpr auto tol = 10 * Kokkos::Experimental::epsilon<magnitude_t<S>>::value;

    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);

    x->putScalar(Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    std::vector<magnitude_t<S>> expt(numCols);
    x->norm2(manitude_view_um_t(expt.data(), numCols));

    // Perform this operation through NOX interface and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox(x);
    std::vector<double> norms(numCols);
    x_nox.norm(norms, NOX::Abstract::Vector::TwoNorm);

    TEST_ABSOLUTE_COMPARE_FLOATING_ARRAYS(norms, expt, tol);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, MaxNorm, S, LO, GO, N)
{
    using manitude_view_um_t = Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    constexpr auto tol = 10 * Kokkos::Experimental::epsilon<magnitude_t<S>>::value;

    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);

    x->putScalar(Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    std::vector<magnitude_t<S>> expt(numCols);
    x->normInf(manitude_view_um_t(expt.data(), numCols));

    // Perform this operation through NOX interface and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox(x);
    std::vector<double> norms(numCols);
    x_nox.norm(norms, NOX::Abstract::Vector::MaxNorm);

    TEST_ABSOLUTE_COMPARE_FLOATING_ARRAYS(norms, expt, tol);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, SetBlock, S, LO, GO, N)
{
    constexpr size_t otherNumCols = 2;

    const std::vector<int>    idxs            {2, 4};
    const std::vector<size_t> idxs_transformed{2, 4};

    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map,      numCols);
    const auto y = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map,      numCols);
    const auto z = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, otherNumCols);

    x->putScalar(    Kokkos::ArithTraits<S>::one());
    y->putScalar(    Kokkos::ArithTraits<S>::one());
    z->putScalar(2 * Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    x->subViewNonConst(idxs_transformed)->assign(*z);

    // Perform this operation through NOX interface and check.
          nox_tpetra_multivector_t<S, LO, GO, N> y_nox(y);
    const nox_tpetra_multivector_t<S, LO, GO, N> z_nox(z);

    y_nox.setBlock(z_nox, idxs);

    const Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    for (size_t idx = 0; idx < numCols; ++idx)
    {
        if (idx == 2 || idx == 4) {
            expt(idx) = 2 * Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements);
        } else {
            expt(idx) = Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements);
        }
    }
    success = checkMultiVectors(x, y, expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, Augment, S, LO, GO, N)
{
    constexpr auto tol = 10 * Kokkos::Experimental::epsilon<magnitude_t<S>>::value;

    constexpr size_t otherNumCols = 2;

    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map,      numCols);
    const auto y = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map,      numCols);
    const auto z = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, otherNumCols);

    x->putScalar(    Kokkos::ArithTraits<S>::one());
    y->putScalar(    Kokkos::ArithTraits<S>::one());
    z->putScalar(2 * Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    const auto x_aug = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols + otherNumCols);
    x_aug->subViewNonConst(Teuchos::Range1D(0, numCols - 1))->assign(*x);
    x_aug->subViewNonConst(Teuchos::Range1D(numCols, numCols + otherNumCols - 1))->assign(*z);

    // Perform this operation through NOX interface and check.
          nox_tpetra_multivector_t<S, LO, GO, N> y_nox(y);
    const nox_tpetra_multivector_t<S, LO, GO, N> z_nox(z);
    y_nox.augment(z_nox);

    const Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace> expt("expected norms", numCols + otherNumCols);
    for (size_t idx = 0; idx < numCols + otherNumCols; ++idx)
    {
        if (idx >= numCols) {
            expt(idx) = 2 * Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements);
        } else {
            expt(idx) = Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements);
        }
    }

    for (size_t idx = 0; idx < numCols + otherNumCols; ++idx)
    {
        TEUCHOS_TEST_FLOATING_EQUALITY(
            (y_nox[idx].norm(NOX::Abstract::Vector::TwoNorm)),
            expt(idx),
            tol, out, success
        );
    }

    success = checkMultiVectors(x_aug, y_nox.getTpetraMultiVector(), expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, SubCopy, S, LO, GO, N)
{
    constexpr size_t otherNumCols = 2;

    const std::vector<int>    idxs            {2, 4};
    const std::vector<size_t> idxs_transformed{2, 4};

    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map,      numCols);
    const auto y = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, otherNumCols);

    x->putScalar(Kokkos::ArithTraits<S>::one());
    y->putScalar(Kokkos::ArithTraits<S>::one());

    // Perform this operation through NOX interface and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox(x);

    const auto x_sub_nox = x_nox.subCopy(idxs);
    const auto x_sub = Teuchos::rcp_dynamic_cast<nox_tpetra_multivector_t<S, LO, GO, N>>(x_sub_nox)->getTpetraMultiVector();

    TEUCHOS_TEST_EQUALITY(x    .strong_count(), 2, out, success);
    TEUCHOS_TEST_EQUALITY(x_sub.strong_count(), 2, out, success);

    const Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    Kokkos::deep_copy(expt, Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements));
    success = checkMultiVectors(x_sub, y, expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, SubView, S, LO, GO, N)
{
    constexpr size_t otherNumCols = 2;

    const std::vector<int>    idxs            {2, 4};
    const std::vector<size_t> idxs_transformed{2, 4};

    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map,      numCols);
    const auto y = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map,      numCols);
    const auto z = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, otherNumCols);

    x->putScalar(    Kokkos::ArithTraits<S>::one());
    y->putScalar(    Kokkos::ArithTraits<S>::one());
    z->putScalar(2 * Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    x->subViewNonConst(idxs_transformed)->assign(*z);

    // Perform this operation through NOX interface and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> y_nox(y);

    y_nox.subView(idxs)->scale(2 * Kokkos::ArithTraits<S>::one());

    const Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    for (size_t idx = 0; idx < numCols; ++idx)
    {
        if (idx == 2 || idx == 4) {
            expt(idx) = 2 * Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements);
        } else {
            expt(idx) = Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements);
        }
    }
    success = checkMultiVectors(x, y, expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, Multiply, S, LO, GO, N)
{
    constexpr auto tol = 10 * Kokkos::Experimental::epsilon<magnitude_t<S>>::value;

    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map,      numCols);
    const auto y = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map,      numCols);

    for (size_t idx = 0; idx <numCols; ++idx)
    {
        x->getVectorNonConst(idx)->putScalar(    idx * Kokkos::ArithTraits<S>::one());
        y->getVectorNonConst(idx)->putScalar(2 * idx * Kokkos::ArithTraits<S>::one());
    }

    // Perform this operation through NOX interface.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox(x);
    const nox_tpetra_multivector_t<S, LO, GO, N> y_nox(y);

    NOX::Abstract::MultiVector::DenseMatrix dots_nox(numCols, numCols);
    x_nox.multiply(Kokkos::ArithTraits<S>::one(), y_nox, dots_nox);

    // Check for correct answer.
    for (size_t idx = 0; idx < numCols; ++idx)
    {
        for (size_t jdx = 0; jdx < numCols; ++jdx)
        {
            TEUCHOS_TEST_FLOATING_EQUALITY(
                (dots_nox(idx, jdx)),
                (2 * jdx * idx * Kokkos::ArithTraits<S>::one() * numGlobalElements),
                tol, out, success
            );
        }
    }
}

#define UNIT_TEST_GROUP(S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, Traits,             S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, ConstructorFromRCP, S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, ConstructorFromRef, S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, CopyConstructor,    S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, CopyAssignment,     S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, Init,               S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, Random,             S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, Scale,              S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, Update_1,           S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, Update_2,           S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, Update_3,           S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, OneNorm,            S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, TwoNorm,            S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, MaxNorm,            S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, SetBlock,           S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, Augment,            S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, SubCopy,            S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, SubView,            S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, Multiply,           S, LO, GO, N)


#include "NOX_Tpetra_ETIHelperMacros.hpp"

TPETRA_ETI_MANGLING_TYPEDEFS()

NOX_TPETRA_INSTANTIATE_SLGN(UNIT_TEST_GROUP)
