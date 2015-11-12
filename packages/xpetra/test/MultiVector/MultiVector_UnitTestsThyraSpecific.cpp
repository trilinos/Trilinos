// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Xpetra_UnitTestHelpers.hpp>
#include <Teuchos_Comm.hpp>

#include "RTOpPack_ROpNorm1.hpp"

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"

#include "Xpetra_MultiVector.hpp"
#include "Xpetra_Vector.hpp"

#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

#include "Xpetra_ThyraUtils.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMultiVector.hpp"
#include "Xpetra_TpetraVector.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMultiVector.hpp"
#include "Xpetra_EpetraVector.hpp"
#endif

#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"

namespace {

  bool testMpi = true;

  TEUCHOS_STATIC_SETUP()
  {
     Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
     clp.addOutputSetupOptions(true);
     clp.setOption(
                   "test-mpi", "test-serial", &testMpi,
                   "Test MPI (if available) or force test of serial.  In a serial build,"
                   " this option is ignored and a serial comm is always used." );
  }

  Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
  {
    Teuchos::RCP<const Teuchos::Comm<int> > ret;
    if (testMpi) {
      ret = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    }
    else {
      ret = rcp(new Teuchos::SerialComm<int>());
    }
    return ret;
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL( MultiVector, Create,        MV, V, Ordinal, Scalar , Node )
  {
    typedef Ordinal LO;
    typedef Ordinal GO;
    typedef Scalar scalar_type;
    typedef Xpetra::Map<LO, GO, Node> map_type;
    typedef Xpetra::MapFactory<LO, GO, Node> map_factory_type;
    typedef Xpetra::MultiVector<Scalar, LO, GO, Node> mv_type;
    typedef Xpetra::ThyraUtils<Scalar, LO, GO, Node> th_utils_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm ();

#ifdef HAVE_XPETRA_TPETRA
    Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;
#else
#  ifdef HAVE_XPETRA_EPETRA
    Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;
#  else
#    error "Should never get here!"
#  endif // HAVE_XPETRA_EPETRA
#endif // HAVE_XPETRA_TPETRA

    // create an Xpetra map
    const LO numInd = 63;
    Teuchos::RCP<const map_type> map = map_factory_type::Build (lib, numInd, 0, comm);

    // create Thyra vector space out of Xpetra Map
    Teuchos::RCP<const Thyra::VectorSpaceBase<scalar_type> > thMap = th_utils_type::toThyra(map);
    TEUCHOS_TEST_FOR_EXCEPTION(map->getGlobalNumElements()!=thMap->dim(), std::logic_error, "Global dimension of Xpetra map and Thyra VectorSpaceBase are different.");
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<scalar_type> > thSpmdMap = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<scalar_type> >(thMap);
    TEUCHOS_TEST_FOR_EXCEPTION(thSpmdMap == Teuchos::null, std::logic_error, "Cannot cast VectorSpaceBase to SpmdVectorSpaceBase.");
    TEUCHOS_TEST_FOR_EXCEPTION(map->getNodeNumElements()!=thSpmdMap->localSubDim(), std::logic_error, "Local dimension of Xpetra map and Thyra VectorSpaceBase on one (or more) processor(s) are different.");

    // create Thyra MultiVector
    Teuchos::RCP< Thyra::MultiVectorBase<scalar_type> > thMVec = Thyra::createMembers(thMap, 2);
    Teuchos::RCP< Thyra::SpmdMultiVectorBase<scalar_type> > thSpmdMVec = Teuchos::rcp_dynamic_cast<Thyra::SpmdMultiVectorBase<scalar_type> >(thMVec);
    TEUCHOS_TEST_FOR_EXCEPTION(thSpmdMVec == Teuchos::null, std::logic_error, "Cannot cast MultiVectorBase to SpmdMultiVectorBase.");

    // fill multivector with some data
    const Ordinal localOffset = ( thSpmdMap != Teuchos::null ? thSpmdMap->localOffset() : 0 );
    const Ordinal localSubDim = ( thSpmdMap != Teuchos::null ? thSpmdMap->localSubDim() : thMap->dim() );
    Teuchos::RCP<Thyra::DetachedMultiVectorView<scalar_type> > thyData =
        Teuchos::rcp(new Thyra::DetachedMultiVectorView<scalar_type>(*thSpmdMVec,Teuchos::Range1D(localOffset,localOffset+localSubDim-1)));

    // loop over all vectors in multivector
    for(Ordinal j = 0; j < thSpmdMVec->domain()->dim(); ++j) {
      // loop over all local rows
      for(Ordinal i = 0; i < localSubDim; ++i) {
        (*thyData)(i,0) = 1;
        (*thyData)(i,1) = 2;
      }
    }

    // calculate and check 1-norm of Thyra MultiVector
    RTOpPack::ROpNorm1<scalar_type> op;
    const Ordinal numVec = thSpmdMVec->domain()->dim();
    Teuchos::Array<Teuchos::RCP<RTOpPack::ReductTarget> > rcp_op_targs(numVec);
    Teuchos::Array<Teuchos::Ptr<RTOpPack::ReductTarget> > op_targs(numVec);
    for( Ordinal kc = 0; kc < numVec; ++kc ) {
      rcp_op_targs[kc] = op.reduct_obj_create();
      op_targs[kc] = rcp_op_targs[kc].ptr();
    }
    ::Thyra::applyOp<scalar_type>(op, Teuchos::tuple(Teuchos::ptrInArg(*(thMVec.ptr()))),
      Teuchos::ArrayView<Teuchos::Ptr<Thyra::MultiVectorBase<scalar_type> > >(Teuchos::null),
      op_targs );
    TEST_EQUALITY( op(*op_targs[0]), numInd );
    TEST_EQUALITY( op(*op_targs[1]), 2*numInd );

    // create Xpetra multivector from Thyra multi vector
    Teuchos::RCP<mv_type> xpMVec = th_utils_type::toXpetra(thMVec,comm);
    TEUCHOS_TEST_FOR_EXCEPTION(xpMVec == Teuchos::null, std::logic_error, "Failed to convert Thyra::MultiVector to Xpetra::MultiVector.");
    TEST_EQUALITY( xpMVec->getNumVectors(), numVec );
    TEST_EQUALITY( xpMVec->getLocalLength(), localSubDim );
    TEST_EQUALITY( xpMVec->getGlobalLength(), numInd );

    std::vector<scalar_type> norms (numVec, STS::zero() );
    Teuchos::ArrayView<scalar_type> normsView(norms);
    xpMVec->norm1(normsView);
    TEST_EQUALITY( normsView[0], numInd );
    TEST_EQUALITY( normsView[1], 2*numInd );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL( MultiVector, CreateProductMV,        MV, V, Ordinal, Scalar , Node )
  {
    typedef Ordinal LO;
    typedef Ordinal GO;
    typedef Scalar scalar_type;
    typedef Xpetra::Map<LO, GO, Node> map_type;
    typedef Xpetra::MapFactory<LO, GO, Node> map_factory_type;
    typedef Xpetra::MultiVector<Scalar, LO, GO, Node> mv_type;
    typedef Xpetra::ThyraUtils<Scalar, LO, GO, Node> th_utils_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm ();

#ifdef HAVE_XPETRA_TPETRA
    Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;
#else
#  ifdef HAVE_XPETRA_EPETRA
    Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;
#  else
#    error "Should never get here!"
#  endif // HAVE_XPETRA_EPETRA
#endif // HAVE_XPETRA_TPETRA

    // create an Xpetra map
    const LO numA = 63;
    const LO numB = 24;
    Teuchos::RCP<const map_type> mapA  = map_factory_type::Build (lib, numA, 0, comm);
    Teuchos::RCP<const map_type> mapB  = map_factory_type::Build (lib, numB, 0, comm);

    // create Thyra vector space out of Xpetra Map
    Teuchos::RCP<const Thyra::VectorSpaceBase<scalar_type> > thMapA = th_utils_type::toThyra(mapA);
    Teuchos::RCP<const Thyra::VectorSpaceBase<scalar_type> > thMapB = th_utils_type::toThyra(mapB);

    TEUCHOS_TEST_FOR_EXCEPTION(mapA->getGlobalNumElements()!=thMapA->dim(), std::logic_error, "Global dimension of Xpetra map and Thyra VectorSpaceBase are different.");
    TEUCHOS_TEST_FOR_EXCEPTION(mapB->getGlobalNumElements()!=thMapB->dim(), std::logic_error, "Global dimension of Xpetra map and Thyra VectorSpaceBase are different.");

    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<scalar_type> > thSpmdMapA = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<scalar_type> >(thMapA);
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<scalar_type> > thSpmdMapB = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<scalar_type> >(thMapB);
    TEUCHOS_TEST_FOR_EXCEPTION(thSpmdMapA == Teuchos::null, std::logic_error, "Cannot cast VectorSpaceBase to SpmdVectorSpaceBase.");
    TEUCHOS_TEST_FOR_EXCEPTION(thSpmdMapB == Teuchos::null, std::logic_error, "Cannot cast VectorSpaceBase to SpmdVectorSpaceBase.");

    // create Thyra MultiVector
    Teuchos::RCP< Thyra::MultiVectorBase<scalar_type> > thMVecA = Thyra::createMembers(thMapA, 1);
    Teuchos::RCP< Thyra::MultiVectorBase<scalar_type> > thMVecB = Thyra::createMembers(thMapB, 1);
    Teuchos::RCP< Thyra::SpmdMultiVectorBase<scalar_type> > thSpmdMVecA = Teuchos::rcp_dynamic_cast<Thyra::SpmdMultiVectorBase<scalar_type> >(thMVecA);
    Teuchos::RCP< Thyra::SpmdMultiVectorBase<scalar_type> > thSpmdMVecB = Teuchos::rcp_dynamic_cast<Thyra::SpmdMultiVectorBase<scalar_type> >(thMVecB);
    TEUCHOS_TEST_FOR_EXCEPTION(thSpmdMVecA == Teuchos::null, std::logic_error, "Cannot cast MultiVectorBase to SpmdMultiVectorBase.");
    TEUCHOS_TEST_FOR_EXCEPTION(thSpmdMVecB == Teuchos::null, std::logic_error, "Cannot cast MultiVectorBase to SpmdMultiVectorBase.");

    // fill multivector A with some data
    const Ordinal localOffsetA = ( thSpmdMapA != Teuchos::null ? thSpmdMapA->localOffset() : 0 );
    const Ordinal localSubDimA = ( thSpmdMapA != Teuchos::null ? thSpmdMapA->localSubDim() : thMapA->dim() );
    Teuchos::RCP<Thyra::DetachedMultiVectorView<scalar_type> > thyDataA =
        Teuchos::rcp(new Thyra::DetachedMultiVectorView<scalar_type>(*thSpmdMVecA,Teuchos::Range1D(localOffsetA,localOffsetA+localSubDimA-1)));

    // loop over all vectors in multivector
    for(Ordinal j = 0; j < thSpmdMVecA->domain()->dim(); ++j) {
      // loop over all local rows
      for(Ordinal i = 0; i < localSubDimA; ++i) {
        (*thyDataA)(i,j) = 1;
      }
    }

    // fill multivector B with some data
    const Ordinal localOffsetB = ( thSpmdMapB != Teuchos::null ? thSpmdMapB->localOffset() : 0 );
    const Ordinal localSubDimB = ( thSpmdMapB != Teuchos::null ? thSpmdMapB->localSubDim() : thMapB->dim() );
    Teuchos::RCP<Thyra::DetachedMultiVectorView<scalar_type> > thyDataB =
        Teuchos::rcp(new Thyra::DetachedMultiVectorView<scalar_type>(*thSpmdMVecB,Teuchos::Range1D(localOffsetB,localOffsetB+localSubDimB-1)));

    // loop over all vectors in multivector
    for(Ordinal j = 0; j < thSpmdMVecB->domain()->dim(); ++j) {
      // loop over all local rows
      for(Ordinal i = 0; i < localSubDimB; ++i) {
        (*thyDataB)(i,j) = 2;
      }
    }

    Teuchos::RCP<Thyra::DefaultProductVectorSpace<scalar_type> > thyProdVecSpace = Thyra::productVectorSpace(Teuchos::tuple(thMapA,thMapB));
    TEUCHOS_TEST_FOR_EXCEPTION(thyProdVecSpace == Teuchos::null, std::logic_error, "Failed to create product vector space.");

    // create product multi vector from multivectors A and B
    Teuchos::RCP<Thyra::DefaultProductMultiVector<scalar_type> > thyProdAB = Thyra::defaultProductMultiVector<scalar_type>(thyProdVecSpace, Teuchos::tuple(thMVecA,thMVecB));
    TEUCHOS_TEST_FOR_EXCEPTION(thyProdAB == Teuchos::null, std::logic_error, "Failed to create product multivector.");
    Teuchos::RCP<Thyra::MultiVectorBase<scalar_type> > thyProdMVec = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<scalar_type> >(thyProdAB);
    TEUCHOS_TEST_FOR_EXCEPTION(thyProdMVec == Teuchos::null, std::logic_error, "Downcast of product multivector to multivector failed.");

    // calculate and check 1-norm of Thyra MultiVector
    RTOpPack::ROpNorm1<scalar_type> op;
    const Ordinal numVec = thSpmdMVecA->domain()->dim();
    Teuchos::Array<Teuchos::RCP<RTOpPack::ReductTarget> > rcp_op_targs(numVec);
    Teuchos::Array<Teuchos::Ptr<RTOpPack::ReductTarget> > op_targs(numVec);
    for( Ordinal kc = 0; kc < numVec; ++kc ) {
      rcp_op_targs[kc] = op.reduct_obj_create();
      op_targs[kc] = rcp_op_targs[kc].ptr();
    }
    ::Thyra::applyOp<scalar_type>(op, Teuchos::tuple(Teuchos::ptrInArg(*(thyProdMVec.ptr()))),
      Teuchos::ArrayView<Teuchos::Ptr<Thyra::MultiVectorBase<scalar_type> > >(Teuchos::null),
      op_targs );
    TEST_EQUALITY( op(*op_targs[0]), numA + 2*numB );

    // create Xpetra multivector from Thyra multi vector
    Teuchos::RCP<mv_type> xpMVec = th_utils_type::toXpetra(thyProdMVec,comm);
    TEUCHOS_TEST_FOR_EXCEPTION(xpMVec == Teuchos::null, std::logic_error, "Downcast of product multivector to multivector failed.");

    std::vector<scalar_type> norms (1, STS::zero() );
    Teuchos::ArrayView<scalar_type> normsView(norms);
    xpMVec->norm1(normsView);

    std::cout << normsView[0] << std::endl;
    /*TEST_EQUALITY( op(*op_targs[1]), 2*numInd );

    // create Xpetra multivector from Thyra multi vector
    Teuchos::RCP<mv_type> xpMVec = th_utils_type::toXpetra(thMVec,comm);
    TEUCHOS_TEST_FOR_EXCEPTION(xpMVec == Teuchos::null, std::logic_error, "Failed to convert Thyra::MultiVector to Xpetra::MultiVector.");
    TEST_EQUALITY( xpMVec->getNumVectors(), numVec );
    TEST_EQUALITY( xpMVec->getLocalLength(), localSubDim );
    TEST_EQUALITY( xpMVec->getGlobalLength(), numInd );

    std::vector<scalar_type> norms (numVec, STS::zero() );
    Teuchos::ArrayView<scalar_type> normsView(norms);
    xpMVec->norm1(normsView);
    TEST_EQUALITY( normsView[0], numInd );
    TEST_EQUALITY( normsView[1], 2*numInd );*/
  }



#ifdef HAVE_XPETRA_TPETRA
#define UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( MV, V, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT( MultiVector, Create, MV, V, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT( MultiVector, CreateProductMV, MV, V, ORDINAL, SCALAR, NODE) \

#else
#define UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( MV, V, ORDINAL, SCALAR, NODE )
#endif // HAVE_XPETRA_TPETRA


//#if defined(HAVE_XPETRA_TPETRA) && defined(HAVE_XPETRA_INT_INT) && defined(HAVE_XPETRA_SERIAL)
#if defined(HAVE_XPETRA_INT_INT) && defined(HAVE_XPETRA_SERIAL)
  typedef Xpetra::TpetraMultiVector<double,int,int> MMultiVector;
  typedef Xpetra::TpetraVector<double,int,int> MVector;
  typedef Kokkos::Compat::KokkosSerialWrapperNode MNode;
  UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( MMultiVector, MVector, int, double, MNode )
#endif

}
