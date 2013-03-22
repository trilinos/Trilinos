/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <cmath>
#include <test_utils/test_Vector.hpp>
#include <test_utils/test_VectorSpace.hpp>
#include <test_utils/test_MatrixGraph.hpp>

#include <fei_Factory.hpp>
#include <snl_fei_Factory.hpp>
#include <fei_Vector_Impl.hpp>

#include <test_utils/LibraryFactory.hpp>

#ifdef HAVE_FEI_AZTECOO
#include <fei_Aztec_LinSysCore.hpp>
#endif
#include <fei_Factory_Trilinos.hpp>

#undef fei_file
#define fei_file "test_Vector.cpp"
#include <fei_ErrMacros.hpp>


test_Vector::test_Vector(MPI_Comm comm)
  : tester(comm)
{
}

test_Vector::~test_Vector()
{
}

int test_Vector::runtests()
{
  //-------------------------------
  // We'll test the vector produced by Factory_Trilinos 
  fei::SharedPtr<fei::Factory> factory_trilinos(new Factory_Trilinos(comm_));

  if (localProc_==0) FEI_COUT << "getting fei::Vector from Factory_Trilinos..."
			      << FEI_ENDL;
  fei::SharedPtr<fei::Vector> fei_vec = create_vector(factory_trilinos);

  vector_test1(fei_vec);

  if (localProc_==0) FEI_COUT << FEI_ENDL;

  return(0);
}

fei::SharedPtr<fei::Vector>
test_Vector::create_vector(fei::SharedPtr<fei::Factory> factory)
{
  testData test_data(localProc_, numProcs_);

  fei::SharedPtr<fei::VectorSpace> vspace =
    test_VectorSpace::create_VectorSpace(comm_, &test_data, localProc_, numProcs_,
					 false, false, (const char*)0, factory);
  int err = vspace->initComplete();
  if (err != 0) {
    FEI_COUT << "ERROR, failed to create valid fei::VectorSpace." << FEI_ENDL;
    throw std::runtime_error("test_Vector::create_vector: ERROR, failed to create valid fei::VectorSpace.");
  }

  if (localProc_==0) FEI_COUT << "  creating fei::Vector instance... ";

  fei::SharedPtr<fei::Vector> vec = factory->createVector(vspace);

  if (localProc_==0) FEI_COUT << "ok" << FEI_ENDL;

  return(vec);
}

void test_Vector::vector_test1(fei::SharedPtr<fei::Vector> fei_vec)
{
  if (localProc_==0)
  FEI_COUT << "  vector_test1: testing fei::Vector with type '"
	   << fei_vec->typeName() << "':"<<FEI_ENDL;

  fei::SharedPtr<fei::VectorSpace> vspace = fei_vec->getVectorSpace();

  std::vector<int> global_offsets;
  vspace->getGlobalIndexOffsets(global_offsets);

  int i, my_first_offset = global_offsets[localProc_];
  int my_last_offset = global_offsets[localProc_+1]-1;
  int num_local_indices = my_last_offset - my_first_offset + 1;

  std::vector<double> coefs(num_local_indices, 1.0);
  std::vector<double> check_coefs(num_local_indices);
  std::vector<int> indices(num_local_indices);
  for(i=0; i<num_local_indices; ++i) {
    indices[i] = my_first_offset + i;
  }

  if (localProc_==0)
    FEI_COUT << "   testing fei::Vector::copyIn/copyOut...";

  int errcode = fei_vec->copyIn(num_local_indices, &indices[0], &coefs[0]);
  if (errcode != 0) {
    throw std::runtime_error("nonzero errcode from fei_vec->copyIn");
  }

  errcode = fei_vec->copyOut(num_local_indices, &indices[0], &check_coefs[0]);
  if (errcode != 0) {
    throw std::runtime_error("nonzero errcode from fei_vec->copyOut");
  }

  if (coefs != check_coefs) {
    throw std::runtime_error("fei_vec->copyOut didn't produce the right coefs");
  }

  if (localProc_==0)
    FEI_COUT << "ok"<<FEI_ENDL << "   testing fei::Vector::putScalar...";

  errcode = fei_vec->putScalar(0.0);

  if (errcode != 0) {
    throw std::runtime_error("nonzero errcode from fei_vec->putScalar");
  }

  errcode = fei_vec->copyOut(num_local_indices, &indices[0], &check_coefs[0]);
  if (errcode != 0) {
    throw std::runtime_error("nonzero errcode from fei_vec->copyOut");
  }

  for(i=0; i<num_local_indices; ++i) {
    if (std::abs(check_coefs[i]) > 1.e-38) {
      throw std::runtime_error("fei_vec->putScalar(0.0) didn't zero the vector");
    }
  }

  if (localProc_==0) FEI_COUT << "ok"<<FEI_ENDL;
}

