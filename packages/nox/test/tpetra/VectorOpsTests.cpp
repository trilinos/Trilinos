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
using nox_tpetra_multivector_t = NOX::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
using nox_tpetra_vector_t = NOX::Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

template <typename Scalar>
using magnitude_t = typename Kokkos::ArithTraits<Scalar>::magnitudeType;

//! Function to check solution.
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
bool checkVectors(
    const Teuchos::RCP<tpetra_vector_t<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& a,
    const Teuchos::RCP<tpetra_vector_t<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& b,
    const magnitude_t<Scalar> expt,
    Teuchos::FancyOStream& out,
    bool success
)
{
    constexpr auto tol = 10 * Kokkos::Experimental::epsilon<magnitude_t<Scalar>>::value;
    b->update(-1 * Kokkos::ArithTraits<Scalar>::one(), *a, Kokkos::ArithTraits<Scalar>::one());
    TEUCHOS_TEST_COMPARE(b->norm2(), <, tol, out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY(a->norm2(), expt, tol, out, success);
    return success;
}

constexpr size_t numLocalElements = 1000;

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, Traits, S, LO, GO, N)
{
    static_assert(   std::is_default_constructible_v<nox_tpetra_vector_t<S, LO, GO, N>>);
    static_assert(   std::is_destructible_v<         nox_tpetra_vector_t<S, LO, GO, N>>);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, ConstructorFromRCP, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    // Construct NOX vector that wraps the Tpetra vector and check.
    const nox_tpetra_vector_t<S, LO, GO, N> x_nox(x);

    TEST_EQUALITY(x_nox.getTpetraVector(), x);
    TEST_EQUALITY(x.strong_count(), 2);

    TEST_ASSERT( ! x_nox.getImplicitWeighting());

    TEST_EQUALITY(x_nox.length(), static_cast<NOX::size_type>(numGlobalElements));
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, ConstructorFromRef, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    x->putScalar(Kokkos::ArithTraits<S>::one());

    // Construct NOX vector with copy of Tpetra vector using the flag NOX::DeepCopy and check.
    const nox_tpetra_vector_t<S, LO, GO, N> x_nox_deep_copy(*x, NOX::DeepCopy);

    TEUCHOS_TEST_INEQUALITY(x_nox_deep_copy.getTpetraVector(), x, out, success);
    TEUCHOS_TEST_EQUALITY(x.strong_count(), 1, out, success);

    TEUCHOS_TEST_ASSERT( ! x_nox_deep_copy.getImplicitWeighting(), out, success);

    const auto expt = Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements);
    success = checkVectors(x_nox_deep_copy.getTpetraVector(), x, expt, out, success);

    // Construct NOX vector with copy of Tpetra vector using the flag NOX::ShapeCopy and check.
    const nox_tpetra_vector_t<S, LO, GO, N> x_nox_shape_copy(*x, NOX::ShapeCopy);

    TEUCHOS_TEST_EQUALITY(x_nox_shape_copy.getTpetraVector()->getMap(), map, out, success);

    TEUCHOS_TEST_ASSERT( ! x_nox_shape_copy.getImplicitWeighting(), out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, CopyConstructor, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    x->putScalar(Kokkos::ArithTraits<S>::one());

    // Construct NOX vector that wraps the Tpetra vector.
    const nox_tpetra_vector_t<S, LO, GO, N> x_nox(x);

    TEUCHOS_TEST_EQUALITY(x.strong_count(), 2, out, success);

    // Copy construct NOX vector using the flag NOX::DeepCopy and check.
    const nox_tpetra_vector_t<S, LO, GO, N> x_nox_deep_copy(x_nox, NOX::DeepCopy);

    TEUCHOS_TEST_INEQUALITY(x_nox_deep_copy.getTpetraVector(), x, out, success);
    TEUCHOS_TEST_EQUALITY(x.strong_count(), 2, out, success);

    TEUCHOS_TEST_ASSERT( ! x_nox_deep_copy.getImplicitWeighting(), out, success);

    const auto expt = Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements);
    success = checkVectors(x_nox_deep_copy.getTpetraVector(), x, expt, out, success);

    // Copy construct NOX vector using the flag NOX::ShapeCopy and check.
    const nox_tpetra_vector_t<S, LO, GO, N> x_nox_shape_copy(x_nox, NOX::ShapeCopy);

    TEUCHOS_TEST_ASSERT( ! x_nox_shape_copy.getImplicitWeighting(), out, success);

    TEUCHOS_TEST_EQUALITY(x_nox_shape_copy.getTpetraVector()->getMap(), map, out, success);

    // After setting a weight vector, copy construct NOX vector using the flag NOX::DeepCopy and check.
    const auto w = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    nox_tpetra_vector_t<S, LO, GO, N> x_nox_weighted(x_nox, NOX::DeepCopy);
    x_nox_weighted.setWeightVector(w);
    TEUCHOS_TEST_ASSERT(x_nox_weighted.hasWeightVector(), out, success);
    TEUCHOS_TEST_EQUALITY(x_nox_weighted.getWeightVector(), w, out, success);
    TEUCHOS_TEST_ASSERT(x_nox_weighted.getImplicitWeighting(), out, success);

    const nox_tpetra_vector_t<S, LO, GO, N> x_nox_weighted_deep_copy(x_nox_weighted, NOX::DeepCopy);
    TEUCHOS_TEST_ASSERT(x_nox_weighted_deep_copy.hasWeightVector(), out, success);
    TEUCHOS_TEST_EQUALITY(x_nox_weighted_deep_copy.getWeightVector(), w, out, success);
    TEUCHOS_TEST_ASSERT(x_nox_weighted_deep_copy.getImplicitWeighting(), out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, CopyAssignment, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);
    const auto y = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    x->putScalar(Kokkos::ArithTraits<S>::one());
    y->putScalar(Kokkos::ArithTraits<S>::zero());

    // Copy x into y through NOX interface.
    const nox_tpetra_vector_t<S, LO, GO, N> x_nox(x);
          nox_tpetra_vector_t<S, LO, GO, N> y_nox(y);

    y_nox = x_nox;

    // Check for correct answer.
    const auto expt = Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements);
    success = checkVectors(y_nox.getTpetraVector(), x_nox.getTpetraVector(), expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, Init, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);
    const auto y = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    // Perform operation with Tpetra directly.
    x->putScalar(2 * Kokkos::ArithTraits<S>::one());

    // Construct NOX vector that wraps Tpetra vector and perform the init operation through NOX interface.
    nox_tpetra_vector_t<S, LO, GO, N> y_nox(y);
    y_nox.init(2 * Kokkos::ArithTraits<S>::one());

    // Check for correct answer.
    const auto expt = Kokkos::ArithTraits<magnitude_t<S>>::squareroot(4 * numGlobalElements);
    success = checkVectors(y_nox.getTpetraVector(), x, expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, Abs, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);
    const auto y = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    x->putScalar(-1 * Kokkos::ArithTraits<S>::one());
    y->putScalar(-1 * Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    x->abs(*x);

    // Perform the abs operation through NOX interface and check.
    nox_tpetra_vector_t<S, LO, GO, N> y_nox(y);
    y_nox.abs(y_nox);

    const auto expt = Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements);
    success = checkVectors(y_nox.getTpetraVector(), x, expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, Random, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    // Randomize through NOX interface.
    nox_tpetra_vector_t<S, LO, GO, N> x_nox(x);

    x_nox.random();

    // Check that the randomly generated entries are within the expected range.
    TEST_ASSERT(x_nox.norm(NOX::Abstract::Vector::NormType::MaxNorm) <= 1.0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, Reciprocal, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);
    const auto y = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    x->putScalar(2 * Kokkos::ArithTraits<S>::one());
    y->putScalar(2 * Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    x->reciprocal(*x);

    // Perform this operation through NOX interface and check.
    nox_tpetra_vector_t<S, LO, GO, N> y_nox(y);

    y_nox.reciprocal(y_nox);

    const auto expt = Kokkos::ArithTraits<magnitude_t<S>>::squareroot(static_cast<magnitude_t<S>>(0.25) * numGlobalElements);
    success = checkVectors(y_nox.getTpetraVector(), x, expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, Scale_Scalar, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);
    const auto y = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    x->putScalar(Kokkos::ArithTraits<S>::one());
    y->putScalar(Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    x->scale(2 * Kokkos::ArithTraits<S>::one());

    // Perform this operation through NOX interface and check.
    nox_tpetra_vector_t<S, LO, GO, N> y_nox(y);

    y_nox.scale(2 * Kokkos::ArithTraits<S>::one());

    const auto expt = Kokkos::ArithTraits<magnitude_t<S>>::squareroot(4 * numGlobalElements);
    success = checkVectors(y_nox.getTpetraVector(), x, expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, Scale_Vector, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);
    const auto y = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);
    const auto z = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    x->putScalar(    Kokkos::ArithTraits<S>::one());
    y->putScalar(    Kokkos::ArithTraits<S>::one());
    z->putScalar(2 * Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    x->elementWiseMultiply(Kokkos::ArithTraits<S>::one(), *z, *x, Kokkos::ArithTraits<S>::zero());

    // Perform this operation through NOX interface and check.
          nox_tpetra_vector_t<S, LO, GO, N> y_nox(y);
    const nox_tpetra_vector_t<S, LO, GO, N> z_nox(z);

    y_nox.scale(z_nox);

    const auto expt = Kokkos::ArithTraits<magnitude_t<S>>::squareroot(4 * numGlobalElements);
    success = checkVectors(y_nox.getTpetraVector(), x, expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, Update_1, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);
    const auto y = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);
    const auto z = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    x->putScalar(Kokkos::ArithTraits<S>::one());
    y->putScalar(Kokkos::ArithTraits<S>::one());
    z->putScalar(Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    x->update(2 * Kokkos::ArithTraits<S>::one(), *z, Kokkos::ArithTraits<S>::one());

    // Perform this operation through NOX interface and check.
          nox_tpetra_vector_t<S, LO, GO, N> y_nox(y);
    const nox_tpetra_vector_t<S, LO, GO, N> z_nox(z);

    y_nox.update(2 * Kokkos::ArithTraits<S>::one(), z_nox, Kokkos::ArithTraits<S>::one());

    const auto expt = Kokkos::ArithTraits<magnitude_t<S>>::squareroot(9 * numGlobalElements);
    success = checkVectors(y_nox.getTpetraVector(), x, expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, Update_2, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto w = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);
    const auto y = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);
    const auto z = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    w->putScalar(Kokkos::ArithTraits<S>::one());
    x->putScalar(Kokkos::ArithTraits<S>::one());
    y->putScalar(Kokkos::ArithTraits<S>::one());
    z->putScalar(Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    x->update(2 * Kokkos::ArithTraits<S>::one(), *w, 2 * Kokkos::ArithTraits<S>::one(), *z, Kokkos::ArithTraits<S>::one());

    // Perform this operation through NOX interface and check.
    const nox_tpetra_vector_t<S, LO, GO, N> w_nox(w);
          nox_tpetra_vector_t<S, LO, GO, N> y_nox(y);
    const nox_tpetra_vector_t<S, LO, GO, N> z_nox(z);

    y_nox.update(2 * Kokkos::ArithTraits<S>::one(), w_nox, 2 * Kokkos::ArithTraits<S>::one(), z_nox, Kokkos::ArithTraits<S>::one());

    const auto expt = Kokkos::ArithTraits<magnitude_t<S>>::squareroot(25 * numGlobalElements);
    success = checkVectors(y_nox.getTpetraVector(), x, expt, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, CreateMultiVector_1, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    constexpr size_t num_additional_vecs = 2;

    // Create Tpetra vectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);
    const auto y = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);
    const auto z = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    x->putScalar(    Kokkos::ArithTraits<S>::one());
    y->putScalar(2 * Kokkos::ArithTraits<S>::one());
    z->putScalar(4 * Kokkos::ArithTraits<S>::one());

    // Construct NOX vectors.
    const nox_tpetra_vector_t<S, LO, GO, N> x_nox(x);
    const nox_tpetra_vector_t<S, LO, GO, N> y_nox(y);
    const nox_tpetra_vector_t<S, LO, GO, N> z_nox(z);

    // Construct NOX multivector using the flag NOX::DeepCopy and check.
    const auto v_nox = Teuchos::rcp_dynamic_cast<nox_tpetra_multivector_t<S, LO, GO, N>>(
        x_nox.createMultiVector(
            std::vector<NOX::Abstract::Vector const*>{&y_nox, &z_nox}.data(),
            num_additional_vecs,
            NOX::DeepCopy
        )
    );

    TEST_EQUALITY(x.strong_count(), 2);
    TEST_EQUALITY(y.strong_count(), 2);
    TEST_EQUALITY(z.strong_count(), 2);

    TEST_EQUALITY(v_nox->getTpetraMultiVector()->getMap(), map);
    TEST_EQUALITY(v_nox->getTpetraMultiVector()->getNumVectors(), 1 + num_additional_vecs);

    TEST_ASSERT( ! v_nox->getImplicitWeighting());

    const Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace> norms("norms vectors", 1 + num_additional_vecs);
    v_nox->getTpetraMultiVector()->norm2(norms);
    constexpr auto tol = 10 * Kokkos::Experimental::epsilon<magnitude_t<S>>::value;
    TEST_FLOATING_EQUALITY(norms(0), Kokkos::ArithTraits<magnitude_t<S>>::squareroot(    numGlobalElements),  tol);
    TEST_FLOATING_EQUALITY(norms(1), Kokkos::ArithTraits<magnitude_t<S>>::squareroot(4 * numGlobalElements),  tol);
    TEST_FLOATING_EQUALITY(norms(2), Kokkos::ArithTraits<magnitude_t<S>>::squareroot(16 * numGlobalElements), tol);

    // Construct NOX multivector using the flag NOX::ShapeCopy and check.
    const auto w_nox = Teuchos::rcp_dynamic_cast<nox_tpetra_multivector_t<S, LO, GO, N>>(
        x_nox.createMultiVector(
            std::vector<NOX::Abstract::Vector const*>{&y_nox, &z_nox}.data(),
            num_additional_vecs,
            NOX::ShapeCopy
        )
    );

    TEST_EQUALITY(x.strong_count(), 2);
    TEST_EQUALITY(y.strong_count(), 2);
    TEST_EQUALITY(z.strong_count(), 2);

    TEST_EQUALITY(w_nox->getTpetraMultiVector()->getMap(), map);
    TEST_EQUALITY(w_nox->getTpetraMultiVector()->getNumVectors(), 1 + num_additional_vecs);

    TEST_ASSERT( ! w_nox->getImplicitWeighting());
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, CreateMultiVector_2, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    constexpr size_t num_vecs = 2;

    // Create Tpetra vector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    x->putScalar(Kokkos::ArithTraits<S>::one());

    // Construct NOX vector.
    const nox_tpetra_vector_t<S, LO, GO, N> x_nox(x);

    // Construct NOX multivector using the flag NOX::DeepCopy and check.
    const auto y_nox = Teuchos::rcp_dynamic_cast<nox_tpetra_multivector_t<S, LO, GO, N>>(
        x_nox.createMultiVector(num_vecs, NOX::DeepCopy)
    );

    TEST_EQUALITY(x.strong_count(), 2);

    TEST_EQUALITY(y_nox->getTpetraMultiVector()->getMap(), map);
    TEST_EQUALITY(y_nox->getTpetraMultiVector()->getNumVectors(), num_vecs);

    TEST_ASSERT( ! y_nox->getImplicitWeighting());

    const Kokkos::View<magnitude_t<S>*, Kokkos::HostSpace> norms("norms vectors", num_vecs);
    y_nox->getTpetraMultiVector()->norm2(norms);
    constexpr auto tol = 10 * Kokkos::Experimental::epsilon<magnitude_t<S>>::value;
    TEST_FLOATING_EQUALITY(norms(0), Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements), tol);
    TEST_FLOATING_EQUALITY(norms(1), Kokkos::ArithTraits<magnitude_t<S>>::squareroot(numGlobalElements), tol);

    // Construct NOX multivector using the flag NOX::ShapeCopy and check.
    const auto z_nox = Teuchos::rcp_dynamic_cast<nox_tpetra_multivector_t<S, LO, GO, N>>(
        x_nox.createMultiVector(num_vecs, NOX::ShapeCopy)
    );

    TEST_EQUALITY(x.strong_count(), 2);

    TEST_EQUALITY(z_nox->getTpetraMultiVector()->getMap(), map);
    TEST_EQUALITY(z_nox->getTpetraMultiVector()->getNumVectors(), num_vecs);

    TEST_ASSERT( ! z_nox->getImplicitWeighting());

    // After setting a weight vector, perform the operations again and check that the weighting object is copied.
    const auto w = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    nox_tpetra_vector_t<S, LO, GO, N> x_nox_weighted(x_nox, NOX::DeepCopy);
    x_nox_weighted.setWeightVector(w);

    const auto y_nox_weighted = Teuchos::rcp_dynamic_cast<nox_tpetra_multivector_t<S, LO, GO, N>>(
        x_nox_weighted.createMultiVector(num_vecs, NOX::DeepCopy)
    );

    TEST_ASSERT(y_nox_weighted->hasWeightVector());
    TEST_EQUALITY(y_nox_weighted->getWeightVector(), w);
    TEST_ASSERT(y_nox_weighted->getImplicitWeighting());

    const auto z_nox_weighted = Teuchos::rcp_dynamic_cast<nox_tpetra_multivector_t<S, LO, GO, N>>(
        x_nox_weighted.createMultiVector(num_vecs, NOX::ShapeCopy)
    );

    TEST_ASSERT(z_nox_weighted->hasWeightVector());
    TEST_EQUALITY(z_nox_weighted->getWeightVector(), w);
    TEST_ASSERT(z_nox_weighted->getImplicitWeighting());
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, OneNorm, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    x->putScalar(Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    const auto norm_tpetra = x->norm1();

    // Perform this operation through NOX interface and check.
    const nox_tpetra_vector_t<S, LO, GO, N> x_nox(x);

    const double norm_nox = x_nox.norm(NOX::Abstract::Vector::OneNorm);

    TEST_FLOATING_EQUALITY(norm_nox, norm_tpetra, 1.0e-14);

    // After setting a weight vector, perform the operation again through the NOX interface.
    // Note that because we use a weight vector that assigns equal weight to all entries,
    // the expected return value is the product of this weight and the unweighted norm.
    const auto w = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    w->putScalar(2 * Kokkos::ArithTraits<S>::one());

    nox_tpetra_vector_t<S, LO, GO, N> x_nox_weighted(x_nox, NOX::DeepCopy);
    x_nox_weighted.setWeightVector(w);

    const double norm_nox_weighted = x_nox_weighted.norm(NOX::Abstract::Vector::OneNorm);

    TEST_FLOATING_EQUALITY(norm_nox_weighted, 2 * norm_tpetra, 1.0e-14);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, TwoNorm, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    x->putScalar(Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    const auto norm_tpetra = x->norm2();

    // Perform this operation through NOX interface and check.
    const nox_tpetra_vector_t<S, LO, GO, N> x_nox(x);

    const double norm_nox = x_nox.norm(NOX::Abstract::Vector::TwoNorm);

    TEST_FLOATING_EQUALITY(norm_nox, norm_tpetra, 1.0e-14);

    // After setting a weight vector, perform the operation again through the NOX interface.
    // Note that because we use a weight vector that assigns equal weight to all entries,
    // the expected return value is the product of this weight and the unweighted norm.
    const auto w = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    w->putScalar(2 * Kokkos::ArithTraits<S>::one());

    nox_tpetra_vector_t<S, LO, GO, N> x_nox_weighted(x_nox, NOX::DeepCopy);
    x_nox_weighted.setWeightVector(w);

    const double norm_nox_weighted = x_nox_weighted.norm(NOX::Abstract::Vector::TwoNorm);

    TEST_FLOATING_EQUALITY(norm_nox_weighted, 2 * norm_tpetra, 1.0e-14);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, MaxNorm, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    x->putScalar(Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    const auto norm_tpetra = x->normInf();

    // Perform this operation through NOX interface and check.
    const nox_tpetra_vector_t<S, LO, GO, N> x_nox(x);

    const double norm_nox = x_nox.norm(NOX::Abstract::Vector::MaxNorm);

    TEST_FLOATING_EQUALITY(norm_nox, norm_tpetra, 1.0e-14);

    // After setting a weight vector, perform the operation again through the NOX interface.
    // Note that because we use a weight vector that assigns equal weight to all entries,
    // the expected return value is the product of this weight and the unweighted norm.
    const auto w = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    w->putScalar(2 * Kokkos::ArithTraits<S>::one());

    nox_tpetra_vector_t<S, LO, GO, N> x_nox_weighted(x_nox, NOX::DeepCopy);
    x_nox_weighted.setWeightVector(w);

    const double norm_nox_weighted = x_nox_weighted.norm(NOX::Abstract::Vector::MaxNorm);

    TEST_FLOATING_EQUALITY(norm_nox_weighted, 2 * norm_tpetra, 1.0e-14);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, NormWeighted, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vectors.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);
    const auto w = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    x->putScalar(    Kokkos::ArithTraits<S>::one());
    w->putScalar(2 * Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    const auto x_deep_copy = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(*x, Teuchos::Copy);
    x_deep_copy->elementWiseMultiply(Kokkos::ArithTraits<S>::one(), *w, *x_deep_copy, Kokkos::ArithTraits<S>::zero());
    const auto norm_tpetra = Kokkos::ArithTraits<S>::squareroot(x_deep_copy->dot(*x));

    // Perform this operation through NOX interface and check.
    const nox_tpetra_vector_t<S, LO, GO, N> x_nox(x);
    const nox_tpetra_vector_t<S, LO, GO, N> w_nox(w);

    const double norm_nox = x_nox.norm(w_nox);

    TEST_FLOATING_EQUALITY(norm_nox, norm_tpetra, 1.0e-14);

    // After setting a weight vector, perform the operation again through the NOX interface.
    // Note that because we use a weight vector that assigns equal weight to all entries,
    // the expected return value is the product of this weight and the unweighted norm.
    const auto v = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    v->putScalar(2 * Kokkos::ArithTraits<S>::one());

    nox_tpetra_vector_t<S, LO, GO, N> x_nox_weighted(x_nox, NOX::DeepCopy);
    x_nox_weighted.setWeightVector(v);

    const double norm_nox_weighted = x_nox_weighted.norm(w_nox);

    TEST_FLOATING_EQUALITY(norm_nox_weighted, 2 * norm_tpetra, 1.0e-14);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_VectorOps, InnerProduct, S, LO, GO, N)
{
    const auto comm = ::Tpetra::getDefaultComm();
    const auto numGlobalElements = static_cast<Tpetra::global_size_t>(comm->getSize() * numLocalElements);

    // Create Tpetra vector.
    const auto map = Teuchos::make_rcp<tpetra_map_t<LO, GO, N>>(numGlobalElements, numLocalElements, 0, comm);
    const auto x = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);
    const auto y = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    x->putScalar(    Kokkos::ArithTraits<S>::one());
    y->putScalar(2 * Kokkos::ArithTraits<S>::one());

    // Perform operation with Tpetra directly.
    const auto dot_tpetra = x->dot(*y);

    // Perform this operation through NOX interface and check.
    const nox_tpetra_vector_t<S, LO, GO, N> x_nox(x);
    const nox_tpetra_vector_t<S, LO, GO, N> y_nox(y);

    const double dot_nox = x_nox.innerProduct(y_nox);

    TEST_FLOATING_EQUALITY(dot_nox, dot_tpetra, 1.0e-14);

    // After setting a weight vector, perform the operation again through the NOX interface.
    // Note that because we use a weight vector that assigns equal weight to all entries,
    // the expected return value is the product of the square of this weight and
    // the unweighted inner product.
    const auto w = Teuchos::make_rcp<tpetra_vector_t<S, LO, GO, N>>(map);

    w->putScalar(2 * Kokkos::ArithTraits<S>::one());

    nox_tpetra_vector_t<S, LO, GO, N> x_nox_weighted(x_nox, NOX::DeepCopy);
    x_nox_weighted.setWeightVector(w);

    const double dot_nox_weighted = x_nox_weighted.innerProduct(y_nox);

    TEST_FLOATING_EQUALITY(dot_nox_weighted, 4 * dot_tpetra, 1.0e-14);
}

#define UNIT_TEST_GROUP(S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, Traits,              S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, ConstructorFromRCP,  S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, ConstructorFromRef,  S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, CopyConstructor,     S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, CopyAssignment,      S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, Init,                S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, Abs,                 S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, Random,              S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, Reciprocal,          S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, Scale_Scalar,        S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, Scale_Vector,        S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, Update_1,            S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, Update_2,            S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, OneNorm,             S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, CreateMultiVector_1, S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, CreateMultiVector_2, S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, TwoNorm,             S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, MaxNorm,             S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, NormWeighted,        S, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_VectorOps, InnerProduct,        S, LO, GO, N)

#include "NOX_Tpetra_ETIHelperMacros.hpp"

TPETRA_ETI_MANGLING_TYPEDEFS()

NOX_TPETRA_INSTANTIATE_SLGN(UNIT_TEST_GROUP)
