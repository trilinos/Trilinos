// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

//#include "Stokhos_ThyraMPVectorUnitTest.hpp"

// Tpetra
#include "Stokhos_Tpetra_MP_Vector.hpp"
#include "Stokhos_Tpetra_Utilities_MP_Vector.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Map.hpp"
//#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
//#include "Tpetra_CrsGraph.hpp"
//#include "Tpetra_CrsMatrix.hpp"

// Thyra
#include "Thyra_TpetraThyraWrappers.hpp"

//#include "Thyra_ScaledModelEvaluator.hpp"
#include "Thyra_Simple2DModelEvaluator.hpp"

#include "Kokkos_Core.hpp"
#include "KokkosCompat_ClassicNodeAPI_Wrapper.hpp"

template <typename scalar, typename ordinal>
inline
scalar generate_vector_coefficient( const ordinal nFEM,
                                    const ordinal nStoch,
                                    const ordinal iColFEM,
                                    const ordinal iStoch )
{
  const scalar X_fem = 100.0 + scalar(iColFEM) / scalar(nFEM);
  const scalar X_stoch =  1.0 + scalar(iStoch) / scalar(nStoch);
  return X_fem + X_stoch;
  //return 1.0;
}

const int VectorSize = 16;
typedef Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Serial> SerialWrapperNode;
typedef Stokhos::DeviceForNode<SerialWrapperNode>::type Device;
typedef Stokhos::StaticFixedStorage<int,double,VectorSize,Device::execution_space> SFS;

typedef SFS Storage;
typedef int LocalOrdinal;
typedef int GlobalOrdinal;
typedef SerialWrapperNode Node;


int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Kokkos::global_sacado_mp_vector_size = VectorSize;

  // Initialize serial
  Kokkos::Serial::initialize();

  // Run tests
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef typename Storage::execution_space Device;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::ScalarTraits<BaseScalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;

  // Ensure device is initialized
  if (!Kokkos::HostSpace::execution_space::is_initialized())
    Kokkos::HostSpace::execution_space::initialize();
  if (!Device::is_initialized())
    Device::initialize();

  // Comm
  RCP<const Tpetra_Comm> comm =
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  // Map
  GlobalOrdinal nrow = 10;
  RCP<Node> node = KokkosClassic::Details::getNode<Node>();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal>(
      nrow, comm, node);
  ArrayView<const GlobalOrdinal> myGIDs = map->getNodeElementList();
  const size_t num_my_row = myGIDs.size();

  // Fill vectors
  RCP<Tpetra_Vector> x1 = Tpetra::createVector<Scalar>(map);
  RCP<Tpetra_Vector> x2 = Tpetra::createVector<Scalar>(map);
  ArrayRCP<Scalar> x1_view = x1->get1dViewNonConst();
  ArrayRCP<Scalar> x2_view = x2->get1dViewNonConst();
  Scalar val1(VectorSize, BaseScalar(0.0)), val2(VectorSize, BaseScalar(0.0));
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      val1.fastAccessCoeff(j) = generate_vector_coefficient<BaseScalar,size_t>(nrow, VectorSize, row, j);
      val2.fastAccessCoeff(j) = 0.12345 * generate_vector_coefficient<BaseScalar,size_t>(nrow, VectorSize, row, j);
    }
    x1_view[i] = val1;
    x2_view[i] = val2;
  }
  x1_view = Teuchos::null;
  x2_view = Teuchos::null;

  //x1->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))), Teuchos::VERB_EXTREME);
  //x2->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))), Teuchos::VERB_EXTREME);

  // Try to create Thyra vectors
  RCP<Thyra::VectorBase<Scalar> > xx1 = Thyra::createVector(x1);
  RCP<Thyra::VectorBase<Scalar> > xx2 = Thyra::createVector(x2);

  //xx1->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))), Teuchos::VERB_EXTREME);
  //xx2->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))), Teuchos::VERB_EXTREME);

  // Add
  Scalar alpha = 2.1;
  Scalar beta = 3.7;
  RCP<Tpetra_Vector> y = Tpetra::createVector<Scalar>(map);
  RCP<Thyra::VectorBase<Scalar> > yy = Thyra::createVector(y);
  //y->update(alpha, *x1, beta, *x2, Scalar(0.0));
  //yy->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))), Teuchos::VERB_EXTREME);
  Thyra::V_StVpStV(yy.ptr(), alpha, *xx1, beta, *xx2);

  //y->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))), Teuchos::VERB_EXTREME);
  //yy->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))), Teuchos::VERB_EXTREME);

  // Check
  ArrayRCP<Scalar> y_view = y->get1dViewNonConst();
  //Scalar yy_view;
  Scalar val(VectorSize, BaseScalar(0.0));
  BaseScalar tol = 1.0e-14;
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      BaseScalar v = generate_vector_coefficient<BaseScalar,size_t>(
        nrow, VectorSize, row, j);
      val.fastAccessCoeff(j) = alpha.coeff(j)*v + 0.12345*beta.coeff(j)*v;
    }
    //yy_view = Thyra::get_ele(*yy, i); 
    //yy_view = Thyra::get_ele(*yy, row); 
    //TEST_EQUALITY( y_view[i].size(), VectorSize );
    //TEST_EQUALITY( yy_view.size(), VectorSize );
    std::cout << y_view[i] << std::endl;
    std::cout << val << std::endl;
    for (LocalOrdinal j=0; j<VectorSize; ++j)
    {
      //TEST_FLOATING_EQUALITY( y_view[i].fastAccessCoeff(j), val.fastAccessCoeff(j), tol );
      //TEST_FLOATING_EQUALITY( yy_view.fastAccessCoeff(j), val.fastAccessCoeff(j), tol );
    }
  }

  Scalar val1b(VectorSize, BaseScalar(0.0)), val2b(VectorSize, BaseScalar(0.0));
  for (LocalOrdinal j=0; j<VectorSize; ++j)
  {
      val1b.fastAccessCoeff(j) = 3.1415*j;
      val2b.fastAccessCoeff(j) = 0.1234+j;
  }
  //std::cout << "val1b = " << val1b << std::endl;
  //std::cout << "val2b = " << val2b << std::endl;

  // Set up and run Model Evaluator using multipoint type
  RCP<Thyra::ModelEvaluator<Scalar> > model_mp = Thyra::simple2DModelEvaluator<Scalar>();
  Thyra::ModelEvaluatorBase::InArgs<Scalar> in_args_mp = model_mp->getNominalValues();
  RCP<Thyra::VectorBase<Scalar> > x_mp = Thyra::createMember(model_mp->get_x_space());
  Thyra::set_ele(0, val1b, x_mp.ptr());
  Thyra::set_ele(1, val2b, x_mp.ptr());
  in_args_mp.set_x(x_mp);
  //x_mp->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))), Teuchos::VERB_EXTREME);

  Thyra::ModelEvaluatorBase::OutArgs<Scalar> out_args_mp = model_mp->createOutArgs();
  RCP<const Thyra::VectorSpaceBase<Scalar> > f_space_mp = model_mp->get_f_space();
  RCP<Thyra::VectorBase<Scalar> > f_mp = Thyra::createMember(f_space_mp);
  RCP<Thyra::LinearOpBase<Scalar> > W_op_mp = model_mp->create_W_op() ;
  out_args_mp.set_f(f_mp);
  out_args_mp.set_W_op(W_op_mp);

  model_mp->evalModel(in_args_mp, out_args_mp);

  //f_mp->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))), Teuchos::VERB_EXTREME);
  //W_op_mp->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))), Teuchos::VERB_EXTREME);

  const RCP<Thyra::SimpleDenseLinearOp<Scalar> > W_sdlo_mp = Teuchos::rcp_dynamic_cast<Thyra::SimpleDenseLinearOp<Scalar> >(W_op_mp, true);
  const RCP<Thyra::MultiVectorBase<Scalar> > W_mv_mp = W_sdlo_mp->getNonconstMultiVector();
  const Thyra::ConstDetachedVectorView<Scalar> f_dv_mp(f_mp);
  const Thyra::ConstDetachedMultiVectorView<Scalar> W_dv_mp(*W_mv_mp);


  // Set up Model Evaluator using BaseScalar type
  RCP<Thyra::ModelEvaluator<BaseScalar> > model = Thyra::simple2DModelEvaluator<BaseScalar>();
  Thyra::ModelEvaluatorBase::InArgs<BaseScalar> in_args = model->getNominalValues();
  RCP<Thyra::VectorBase<BaseScalar> > x = Thyra::createMember(model->get_x_space());
  in_args.set_x(x);

  Thyra::ModelEvaluatorBase::OutArgs<BaseScalar> out_args = model->createOutArgs();
  RCP<const Thyra::VectorSpaceBase<BaseScalar> > f_space = model->get_f_space();
  RCP<Thyra::VectorBase<BaseScalar> > f = Thyra::createMember(f_space);
  RCP<Thyra::LinearOpBase<BaseScalar> > W_op = model->create_W_op() ;
  out_args.set_f(f);
  out_args.set_W_op(W_op);

  const RCP<Thyra::SimpleDenseLinearOp<BaseScalar> > W_sdlo = Teuchos::rcp_dynamic_cast<Thyra::SimpleDenseLinearOp<BaseScalar> >(W_op, true);
  const RCP<Thyra::MultiVectorBase<BaseScalar> > W_mv = W_sdlo->getNonconstMultiVector();
  const Thyra::ConstDetachedVectorView<BaseScalar> f_dv(f);
  const Thyra::ConstDetachedMultiVectorView<BaseScalar> W_dv(*W_mv);

  // Evaluate ModelEvaluator<BaseScalar> at multiple points and compare with
  // results of ModelEvaluator<Scalar> run
  ScalarMag tolb = Teuchos::as<ScalarMag>(10.0) * ST::eps();
  for (LocalOrdinal j=0; j<VectorSize; ++j)
  {
    Thyra::set_ele(0, val1b.fastAccessCoeff(j), x.ptr());
    Thyra::set_ele(1, val2b.fastAccessCoeff(j), x.ptr());

    model->evalModel(in_args, out_args);

    //x->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))), Teuchos::VERB_EXTREME);
    //f->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))), Teuchos::VERB_EXTREME);
    //W_op->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))), Teuchos::VERB_EXTREME);

    //TEST_FLOATING_EQUALITY( f_dv_mp(0).fastAccessCoeff(j), f_dv(0), tolb );
    //TEST_FLOATING_EQUALITY( f_dv_mp(1).fastAccessCoeff(j), f_dv(1), tolb );
    //TEST_FLOATING_EQUALITY(W_dv_mp(0,0).fastAccessCoeff(j), W_dv(0,0), tolb);
    //TEST_FLOATING_EQUALITY(W_dv_mp(0,1).fastAccessCoeff(j), W_dv(0,1), tolb);
    //TEST_FLOATING_EQUALITY(W_dv_mp(1,0).fastAccessCoeff(j), W_dv(1,0), tolb);
    //TEST_FLOATING_EQUALITY(W_dv_mp(1,1).fastAccessCoeff(j), W_dv(1,1), tolb);

    std::cout << f_dv_mp(0).fastAccessCoeff(j) << " == " << f_dv(0) << std::endl;
    std::cout << f_dv_mp(1).fastAccessCoeff(j) << " == " <<  f_dv(1) << std::endl;
    std::cout << W_dv_mp(0,0).fastAccessCoeff(j) << " == " <<  W_dv(0,0) << std::endl;
    std::cout << W_dv_mp(0,1).fastAccessCoeff(j) << " == " <<  W_dv(0,1) << std::endl;
    std::cout << W_dv_mp(1,0).fastAccessCoeff(j) << " == " <<  W_dv(1,0) << std::endl;
    std::cout << W_dv_mp(1,1).fastAccessCoeff(j) << " == " <<  W_dv(1,1) << std::endl;
  }

  // Finish up
  Kokkos::Serial::finalize();

  return 0;
}
