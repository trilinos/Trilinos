// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Thyra_VectorSpaceTester.hpp"
#include "Thyra_VectorStdOpsTester.hpp"

//#define THYRA_DEFAULT_SPMD_VECTOR_SPACE_UNIT_TESTS_DUMP

#ifdef THYRA_DEFAULT_SPMD_VECTOR_SPACE_UNIT_TESTS_DUMP
#include "RTOpPack_SPMD_apply_op_decl.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#endif  // THYRA_DEFAULT_SPMD_VECTOR_SPACE_UNIT_TESTS_DUMP

namespace {

//
// Helper code and declarations
//

using Teuchos::as;
using Teuchos::get_extra_data;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ScalarTraits;
using Thyra::ConstDetachedSpmdVectorView;
using Thyra::ConstDetachedVectorView;
using Thyra::createMember;
using Thyra::createMembers;
using Thyra::DefaultSpmdVectorSpace;
using Thyra::defaultSpmdVectorSpace;
using Thyra::DetachedSpmdVectorView;
using Thyra::DetachedVectorView;
using Thyra::locallyReplicatedDefaultSpmdVectorSpace;
using Thyra::MultiVectorBase;
using Thyra::VectorBase;
using Thyra::VectorSpaceBase;
typedef Thyra::Ordinal Ordinal;

const int g_localDim = 4;  // ToDo: Make variable!

template <class Scalar>
RCP<VectorSpaceBase<Scalar> >
createSpmdVectorSpace(const Teuchos_Ordinal localDim) {
  return defaultSpmdVectorSpace<Scalar>(
      Teuchos::DefaultComm<Teuchos_Ordinal>::getComm(),
      localDim, -1);
}

//
// Unit Tests
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DefaultSpmdVectorSpace, defaultConstruct,
                                  Scalar) {
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
           Thyra::defaultSpmdVectorSpace<Scalar>());
  TEST_EQUALITY(vs->getComm(), null);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(-1));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(-1));
  TEST_EQUALITY(vs->mapCode(), as<Ordinal>(-1));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(-1));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(DefaultSpmdVectorSpace,
                                                  defaultConstruct)

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DefaultSpmdVectorSpace, serialConstruct,
                                  Scalar) {
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
           defaultSpmdVectorSpace<Scalar>(g_localDim));
  TEST_EQUALITY(vs->getComm(), null);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(0));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->mapCode(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(g_localDim));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(DefaultSpmdVectorSpace,
                                                  serialConstruct)

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DefaultSpmdVectorSpace, serialConstructZeroSize,
                                  Scalar) {
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
           defaultSpmdVectorSpace<Scalar>(0));
  TEST_EQUALITY_CONST(vs->getComm(), null);
  TEST_EQUALITY_CONST(vs->localOffset(), as<Ordinal>(0));
  TEST_EQUALITY_CONST(vs->localSubDim(), as<Ordinal>(0));
  TEST_EQUALITY_CONST(vs->mapCode(), as<Ordinal>(0));
  TEST_EQUALITY_CONST(vs->dim(), as<Ordinal>(0));

  ECHO(const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs));
  out << "v = " << *v;

  TEST_ASSERT(vs->isCompatible(*v->space()));
  TEST_EQUALITY_CONST(v->space()->dim(), as<Ordinal>(0));

  ECHO(const RCP<MultiVectorBase<Scalar> > mv = createMembers<Scalar>(vs, 0));
  out << "mv = " << *mv;

  TEST_ASSERT(vs->isCompatible(*mv->range()));
  TEST_EQUALITY_CONST(mv->range()->dim(), as<Ordinal>(0));
  TEST_EQUALITY_CONST(mv->domain()->dim(), as<Ordinal>(0));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(DefaultSpmdVectorSpace,
                                                  serialConstructZeroSize)

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DefaultSpmdVectorSpace, parallelConstruct,
                                  Scalar) {
  ECHO(const RCP<const Teuchos::Comm<Ordinal> > comm =
           Teuchos::DefaultComm<Teuchos_Ordinal>::getComm());
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
           defaultSpmdVectorSpace<Scalar>(comm, g_localDim, -1));
  TEST_EQUALITY(vs->getComm(), comm);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(comm->getRank() * g_localDim));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(comm->getSize() * g_localDim));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(DefaultSpmdVectorSpace,
                                                  parallelConstruct)

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DefaultSpmdVectorSpace, parallelConstructGlobalDim,
                                  Scalar) {
  ECHO(const RCP<const Teuchos::Comm<Ordinal> > comm =
           Teuchos::DefaultComm<Teuchos_Ordinal>::getComm());
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
           defaultSpmdVectorSpace<Scalar>(comm, g_localDim, g_localDim * comm->getSize()));
  TEST_EQUALITY(vs->getComm(), comm);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(comm->getRank() * g_localDim));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  // TEST_EQUALITY_CONST(vs->isLocallyReplicated(), false);
  TEST_EQUALITY(vs->isLocallyReplicated(), (comm->getSize() == 1));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(comm->getSize() * g_localDim));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(DefaultSpmdVectorSpace,
                                                  parallelConstructGlobalDim)

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DefaultSpmdVectorSpace, locallyReplicatedParallelConstruct,
                                  Scalar) {
  ECHO(const RCP<const Teuchos::Comm<Ordinal> > comm =
           Teuchos::DefaultComm<Teuchos_Ordinal>::getComm());
  ECHO(const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
           locallyReplicatedDefaultSpmdVectorSpace<Scalar>(comm, g_localDim));
  TEST_EQUALITY(vs->getComm(), comm);
  TEST_EQUALITY_CONST(vs->localOffset(), as<Ordinal>(0));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  TEST_EQUALITY_CONST(vs->isLocallyReplicated(), true);
  TEST_EQUALITY(vs->dim(), as<Ordinal>(g_localDim));
  ECHO(const RCP<const VectorSpaceBase<Scalar> > vs_clone =
           vs->clone());
  TEST_EQUALITY(vs_clone->dim(), as<Ordinal>(g_localDim));
  TEST_ASSERT(vs->isCompatible(*vs_clone));
  TEST_ASSERT(vs_clone->isCompatible(*vs));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(DefaultSpmdVectorSpace,
                                                  locallyReplicatedParallelConstruct)

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DefaultSpmdVectorSpace,
                                  locallyReplicatedDefaultSpmdVectorSpace_explicit, Scalar) {
  ECHO(const RCP<const Teuchos::Comm<Ordinal> > comm =
           Teuchos::DefaultComm<Teuchos_Ordinal>::getComm());
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
           defaultSpmdVectorSpace<Scalar>(comm, g_localDim, g_localDim, true));
  TEST_EQUALITY(vs->getComm(), comm);
  TEST_EQUALITY_CONST(vs->localOffset(), as<Ordinal>(0));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  TEST_EQUALITY_CONST(vs->isLocallyReplicated(), true);
  TEST_EQUALITY(vs->dim(), as<Ordinal>(g_localDim));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(DefaultSpmdVectorSpace,
                                                  locallyReplicatedDefaultSpmdVectorSpace_explicit)

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DefaultSpmdVectorSpace,
                                  locallyReplicatedDefaultSpmdVectorSpace_implicit, Scalar) {
  ECHO(const RCP<const Teuchos::Comm<Ordinal> > comm =
           Teuchos::DefaultComm<Teuchos_Ordinal>::getComm());
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
           defaultSpmdVectorSpace<Scalar>(comm, g_localDim, g_localDim));
  TEST_EQUALITY(vs->getComm(), comm);
  TEST_EQUALITY_CONST(vs->localOffset(), as<Ordinal>(0));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  TEST_EQUALITY_CONST(vs->isLocallyReplicated(), true);
  TEST_EQUALITY(vs->dim(), as<Ordinal>(g_localDim));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(DefaultSpmdVectorSpace,
                                                  locallyReplicatedDefaultSpmdVectorSpace_implicit)

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DefaultSpmdVectorSpace,
                                  locallyReplicatedDefaultSpmdVectorSpace_explicit_bad, Scalar) {
  ECHO(const RCP<const Teuchos::Comm<Ordinal> > comm =
           Teuchos::DefaultComm<Teuchos_Ordinal>::getComm());
  ECHO(const Ordinal rank = comm->getRank());
  out << "rank = " << rank << "\n";
  ECHO(const Ordinal globalDim = g_localDim);
  ECHO(const Ordinal localDim = (rank == 0 ? globalDim + 1 : globalDim));
  out << "globalDim = " << globalDim << "\n";
  out << "localDim = " << localDim << "\n";
  if (localDim != globalDim) {
    // Throws when localDim != globalDim on a process!
    TEST_THROW(defaultSpmdVectorSpace<Scalar>(comm, localDim, globalDim, true),
               std::logic_error);
  } else {
    ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
             defaultSpmdVectorSpace<Scalar>(comm, localDim, g_localDim, true));
  }
  // NOTE: There should not be any hanging because no global reductions should
  // be done when you pass in isLocallyReplicated==true!
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(DefaultSpmdVectorSpace,
                                                  locallyReplicatedDefaultSpmdVectorSpace_explicit_bad)

//#ifndef THYRA_HIDE_DEPRECATED_CODE

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DefaultSpmdVectorSpace, deprecatedSerialConstruct,
                                  Scalar) {
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
           Thyra::defaultSpmdVectorSpace<Scalar>(g_localDim));
  TEST_EQUALITY(vs->getComm(), null);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(0));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->mapCode(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(g_localDim));
  ECHO(const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs));
  TEST_ASSERT(vs->isCompatible(*v->space()));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(DefaultSpmdVectorSpace,
                                                  deprecatedSerialConstruct)

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DefaultSpmdVectorSpace, deprecatedParallelConstruct,
                                  Scalar) {
  ECHO(const RCP<const Teuchos::Comm<Ordinal> > comm =
           Teuchos::DefaultComm<Teuchos_Ordinal>::getComm());
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
           Thyra::defaultSpmdVectorSpace<Scalar>(comm, g_localDim, -1));
  TEST_EQUALITY(vs->getComm(), comm);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(comm->getRank() * g_localDim));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(comm->getSize() * g_localDim));
  ECHO(const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs));
  TEST_ASSERT(vs->isCompatible(*v->space()));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(DefaultSpmdVectorSpace,
                                                  deprecatedParallelConstruct)

//#endif // THYRA_HIDE_DEPRECATED_CODE

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DefaultSpmdVectorSpace, parallelFullExtract,
                                  Scalar) {
  const RCP<const VectorSpaceBase<Scalar> > vs =
      createSpmdVectorSpace<Scalar>(g_localDim);

  const RCP<VectorBase<Scalar> > v = createMember(vs);
#ifdef THYRA_DEFAULT_SPMD_VECTOR_SPACE_UNIT_TESTS_DUMP
  RTOpPack::show_spmd_apply_op_dump = true;
#endif
  {
    out << "\nSetting up v[i] = i, i=0...n-1 ...\n";
    const DetachedSpmdVectorView<Scalar> dv(v);
    const Ordinal localOffset = dv.spmdSpace()->localOffset();
    const Ordinal localSubDim = dv.spmdSpace()->localSubDim();
    for (Ordinal i = 0; i < localSubDim; ++i) {
      dv[i] = as<Scalar>(localOffset + i);
    }
  }

  out << "\nv = " << *v;

  {
    const ConstDetachedVectorView<Scalar> dv(v);

    TEST_EQUALITY(dv.subDim(), vs->dim());

    out << "\nTest that dv[i] == i, i=0...n-1 ... ";
    bool local_success = true;
    for (Ordinal i = 0; i < dv.subDim(); ++i) {
      TEST_ARRAY_ELE_EQUALITY(dv, i, as<Scalar>(i));
    }
    if (local_success)
      out << "passed\n";
    else
      success = false;
  }

#ifdef THYRA_DEFAULT_SPMD_VECTOR_SPACE_UNIT_TESTS_DUMP
  RTOpPack::show_spmd_apply_op_dump = false;
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(DefaultSpmdVectorSpace,
                                                  parallelFullExtract)

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DefaultSpmdVectorSpace, dangling_vs,
                                  Scalar) {
  RCP<const Thyra::VectorSpaceBase<Scalar> > vs =
      defaultSpmdVectorSpace<Scalar>(g_localDim);
  const int globalDim         = vs->dim();
  RCP<VectorBase<Scalar> > x1 = createMember(vs);
  {
    // x1 owns a false RCP to vs
    TEST_EQUALITY_CONST(x1->space().has_ownership(), true);
    // RCP<> for x1 owns a true RCP to vs
    const std::string label = "VectorSpaceBase";
    RCP<const VectorSpaceBase<Scalar> > extra_data_x1 =
        get_extra_data<RCP<const VectorSpaceBase<Scalar> >, VectorBase<Scalar> >(x1, label);
    TEST_EQUALITY_CONST(extra_data_x1.has_ownership(), true);
  }
  RCP<Thyra::VectorBase<Scalar> > x0 = x1->clone_v();
  {
    // x0 owns a false RCP to vs
    TEST_EQUALITY_CONST(x0->space().has_ownership(), true);
    // RCP<> for x0 owns a true RCP to a _DIFFERENT_ VectorSpaceBase
    // object because the one used to clone x1 is a false RCP, so the
    // VectorSpaceBase was cloned and that is the one that was set on the RCP.
    std::string label = "VectorSpaceBase";
    RCP<const VectorSpaceBase<Scalar> > extra_data_x0 =
        get_extra_data<RCP<const VectorSpaceBase<Scalar> >, VectorBase<Scalar> >(x0, label);
    TEST_EQUALITY_CONST(extra_data_x0.has_ownership(), true);
    TEST_EQUALITY(extra_data_x0.ptr(), vs.ptr());
  }
  vs = null;  // vs still around because x1's RCP owns it
  x1 = null;  // vs deleted
  {
    RCP<const VectorSpaceBase<Scalar> > vs_old = x0->space();
    TEST_EQUALITY_CONST(vs_old->dim(), globalDim);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(DefaultSpmdVectorSpace,
                                                  dangling_vs)

}  // namespace
