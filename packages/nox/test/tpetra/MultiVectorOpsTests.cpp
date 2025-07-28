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
using magnitude_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;

//! Function to check solution.
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
bool checkMultiVectors(
    const Teuchos::RCP<tpetra_multivector_t<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& a,
    const Teuchos::RCP<tpetra_multivector_t<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& b,
    const Kokkos::View<magnitude_type<Scalar>*, Kokkos::HostSpace>& expt,
    Teuchos::FancyOStream& out,
    bool success
)
{
    constexpr magnitude_type<Scalar> tol = 1.0e-14;

    TEUCHOS_TEST_EQUALITY(a->getNumVectors(), b->getNumVectors(), out, success);

    b->update(-1 * Teuchos::ScalarTraits<Scalar>::one(), *a, Teuchos::ScalarTraits<Scalar>::one());
    const Kokkos::View<magnitude_type<Scalar>*, Kokkos::HostSpace> norms_b("norms b after update", b->getNumVectors());
    b->norm2(norms_b);
    for (size_t idx = 0; idx < norms_b.extent(0); ++idx) {
        TEUCHOS_TEST_COMPARE(norms_b(idx), <, tol, out, success);
    }

    Kokkos::View<magnitude_type<Scalar>*, Kokkos::HostSpace> norms_a("norms a", a->getNumVectors());
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
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, ConstructorFromRef, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);

    x->putScalar(Teuchos::ScalarTraits<S>::one());

    // Construct NOX multivector with copy of Tpetra multivector using the flag NOX::DeepCopy and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox_deep_copy(*x, NOX::DeepCopy);

    TEUCHOS_TEST_INEQUALITY(x_nox_deep_copy.getTpetraMultiVector(), x, out, success);
    TEUCHOS_TEST_EQUALITY(x.strong_count(), 1, out, success);

    TEUCHOS_TEST_ASSERT( ! x_nox_deep_copy.getImplicitWeighting(), out, success);

    const Kokkos::View<magnitude_type<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    Kokkos::deep_copy(expt, Teuchos::ScalarTraits<magnitude_type<S>>::squareroot(numGlobalElements));
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

    x->putScalar(Teuchos::ScalarTraits<S>::one());

    // Construct NOX vector that wraps the Tpetra multivector.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox(x);

    TEUCHOS_TEST_EQUALITY(x.strong_count(), 2, out, success);

    // Copy construct NOX multivector using the flag NOX::DeepCopy and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox_deep_copy(x_nox, NOX::DeepCopy);

    TEUCHOS_TEST_INEQUALITY(x_nox_deep_copy.getTpetraMultiVector(), x, out, success);
    TEUCHOS_TEST_EQUALITY(x.strong_count(), 2, out, success);

    TEUCHOS_TEST_ASSERT( ! x_nox_deep_copy.getImplicitWeighting(), out, success);

    const Kokkos::View<magnitude_type<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    Kokkos::deep_copy(expt, Teuchos::ScalarTraits<magnitude_type<S>>::squareroot(numGlobalElements));
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

    x->putScalar(Teuchos::ScalarTraits<S>::one());
    y->putScalar(Teuchos::ScalarTraits<S>::zero());

    // Copy x into y through NOX interface.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox(x);
          nox_tpetra_multivector_t<S, LO, GO, N> y_nox(y);

    y_nox = x_nox;

    // Check for correct answer.
    const Kokkos::View<magnitude_type<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    Kokkos::deep_copy(expt, Teuchos::ScalarTraits<magnitude_type<S>>::squareroot(numGlobalElements));
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
    x->putScalar(2 * Teuchos::ScalarTraits<S>::one());

    // Perform this operation through NOX interface.
    nox_tpetra_multivector_t<S, LO, GO, N> y_nox(y);
    y_nox.init(2 * Teuchos::ScalarTraits<S>::one());

    // Check for correct answer.
    const Kokkos::View<magnitude_type<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    Kokkos::deep_copy(expt, Teuchos::ScalarTraits<magnitude_type<S>>::squareroot(4 * numGlobalElements));
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

    x->putScalar(Teuchos::ScalarTraits<S>::one());
    y->putScalar(Teuchos::ScalarTraits<S>::one());

    // Perform operation with Tpetra directly.
    x->scale(2 * Teuchos::ScalarTraits<S>::one());

    // Perform this operation through NOX interface and check.
    nox_tpetra_multivector_t<S, LO, GO, N> y_nox(y);
    y_nox.scale(2 * Teuchos::ScalarTraits<S>::one());

    const Kokkos::View<magnitude_type<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    Kokkos::deep_copy(expt, Teuchos::ScalarTraits<magnitude_type<S>>::squareroot(4 * numGlobalElements));
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

    x->putScalar(Teuchos::ScalarTraits<S>::one());
    y->putScalar(Teuchos::ScalarTraits<S>::one());
    z->putScalar(Teuchos::ScalarTraits<S>::one());
    
    // Perform operation with Tpetra directly.
    x->update(2 * Teuchos::ScalarTraits<S>::one(), *z, Teuchos::ScalarTraits<S>::one());
    
    // Perform this operation through NOX interface and check.
          nox_tpetra_multivector_t<S, LO, GO, N> y_nox(y);
    const nox_tpetra_multivector_t<S, LO, GO, N> z_nox(z);

    y_nox.update(2 * Teuchos::ScalarTraits<S>::one(), z_nox, Teuchos::ScalarTraits<S>::one());
    
    const Kokkos::View<magnitude_type<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    Kokkos::deep_copy(expt, Teuchos::ScalarTraits<magnitude_type<S>>::squareroot(9 * numGlobalElements));
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

    w->putScalar(Teuchos::ScalarTraits<S>::one());
    x->putScalar(Teuchos::ScalarTraits<S>::one());
    y->putScalar(Teuchos::ScalarTraits<S>::one());
    z->putScalar(Teuchos::ScalarTraits<S>::one());
    
    // Perform operation with Tpetra directly.
    x->update(2 * Teuchos::ScalarTraits<S>::one(), *w, 2 * Teuchos::ScalarTraits<S>::one(), *z, Teuchos::ScalarTraits<S>::one());
    
    // Perform this operation through NOX interface and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> w_nox(w);
          nox_tpetra_multivector_t<S, LO, GO, N> y_nox(y);
    const nox_tpetra_multivector_t<S, LO, GO, N> z_nox(z);
    
    y_nox.update(2 * Teuchos::ScalarTraits<S>::one(), w_nox, 2 * Teuchos::ScalarTraits<S>::one(), z_nox, Teuchos::ScalarTraits<S>::one());
    
    const Kokkos::View<magnitude_type<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    Kokkos::deep_copy(expt, Teuchos::ScalarTraits<magnitude_type<S>>::squareroot(25 * numGlobalElements));
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

    x->putScalar(Teuchos::ScalarTraits<S>::one());
    y->putScalar(Teuchos::ScalarTraits<S>::one());
    z->putScalar(Teuchos::ScalarTraits<S>::one());
    
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

    x->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, Teuchos::ScalarTraits<S>::one(), *z, *mat_as_mv, Teuchos::ScalarTraits<S>::one());
    
    // Perform this operation through NOX interface and check.
          nox_tpetra_multivector_t<S, LO, GO, N> y_nox(y);
    const nox_tpetra_multivector_t<S, LO, GO, N> z_nox(z);
    
    y_nox.update(Teuchos::NO_TRANS, Teuchos::ScalarTraits<S>::one(), z_nox, mat, Teuchos::ScalarTraits<S>::one());

    const Kokkos::View<magnitude_type<S>*, Kokkos::HostSpace> expt("expected norms", numCols);
    for (size_t idx = 0; idx < numCols; ++idx) {
        expt(idx) = (1 + 5 * idx) * Teuchos::ScalarTraits<magnitude_type<S>>::squareroot(numGlobalElements);
    }
    success = checkMultiVectors(x, y, expt, out, success);

    // Perform this operation over again, this time using the transpose, and check.
    NOX::Abstract::MultiVector::DenseMatrix mat_transpose(mat, Teuchos::TRANS);

    y->putScalar(Teuchos::ScalarTraits<S>::one());
    y_nox.update(Teuchos::TRANS, Teuchos::ScalarTraits<S>::one(), z_nox, mat_transpose, Teuchos::ScalarTraits<S>::one());
    success = checkMultiVectors(x, y, expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, OneNorm, S, LO, GO, N)
{
    using manitude_view_um_t = Kokkos::View<magnitude_type<S>*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    
    constexpr magnitude_type<S> tol = 1.0e-14;
    
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);

    x->putScalar(Teuchos::ScalarTraits<S>::one());

    // Perform operation with Tpetra directly.
    std::vector<magnitude_type<S>> expt(numCols);
    x->norm1(manitude_view_um_t(expt.data(), numCols));

    // Perform this operation through NOX interface and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox(x);
    std::vector<double> norms(numCols);
    x_nox.norm(norms, NOX::Abstract::Vector::OneNorm);

    TEST_ABSOLUTE_COMPARE_FLOATING_ARRAYS(norms, expt, tol);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, TwoNorm, S, LO, GO, N)
{
    using manitude_view_um_t = Kokkos::View<magnitude_type<S>*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    
    constexpr magnitude_type<S> tol = 1.0e-14;
    
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);

    x->putScalar(Teuchos::ScalarTraits<S>::one());

    // Perform operation with Tpetra directly.
    std::vector<magnitude_type<S>> expt(numCols);
    x->norm2(manitude_view_um_t(expt.data(), numCols));

    // Perform this operation through NOX interface and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox(x);
    std::vector<double> norms(numCols);
    x_nox.norm(norms, NOX::Abstract::Vector::TwoNorm);

    TEST_ABSOLUTE_COMPARE_FLOATING_ARRAYS(norms, expt, tol);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MultiVectorOps, MaxNorm, S, LO, GO, N)
{
    using manitude_view_um_t = Kokkos::View<magnitude_type<S>*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    
    constexpr magnitude_type<S> tol = 1.0e-14;
    
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra multivector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_multivector_t<S, LO, GO, N>>(map, numCols);

    x->putScalar(Teuchos::ScalarTraits<S>::one());

    // Perform operation with Tpetra directly.
    std::vector<magnitude_type<S>> expt(numCols);
    x->normInf(manitude_view_um_t(expt.data(), numCols));

    // Perform this operation through NOX interface and check.
    const nox_tpetra_multivector_t<S, LO, GO, N> x_nox(x);
    std::vector<double> norms(numCols);
    x_nox.norm(norms, NOX::Abstract::Vector::MaxNorm);

    TEST_ABSOLUTE_COMPARE_FLOATING_ARRAYS(norms, expt, tol);
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
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MultiVectorOps, MaxNorm,            S, LO, GO, N)

#include "NOX_Tpetra_ETIHelperMacros.hpp"

TPETRA_ETI_MANGLING_TYPEDEFS()

NOX_TPETRA_INSTANTIATE_SLGN(UNIT_TEST_GROUP)
