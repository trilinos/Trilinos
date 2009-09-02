/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_iostream.hpp>

#include <fei_utils.hpp>
#include <test_utils/fei_test_utils.hpp>

#include <test_utils/test_Factory.hpp>
#include <test_utils/test_Vector.hpp>
#include <test_utils/test_Matrix.hpp>

#include <fei_Factory_Trilinos.hpp>

#undef fei_file
#define fei_file "fei_test_plugins.cpp"
#include <fei_ErrMacros.hpp>


//--------------------------------------------------
// main
//--------------------------------------------------
int main(int argc, char** argv) {

  int numProcs, localProc;
  CHK_ERR( fei_test_utils::initialize_mpi(argc, argv, localProc, numProcs) );

  if (localProc == 0) {
    FEI_COUT << "FEI version: " << fei::utils::version() << FEI_ENDL;
  }

  fei::SharedPtr<fei::Factory> factory(new Factory_Trilinos(MPI_COMM_WORLD));

  fei::ParameterSet pset;
  pset.add(fei::Param("LPM_EpetraBasic", true));
  factory->parameters(pset);

  test_Factory ftester(MPI_COMM_WORLD);

  ftester.factory_test1(factory);

  test_Vector vtester(MPI_COMM_WORLD);

  fei::SharedPtr<fei::Vector> fei_vec = vtester.create_vector(factory);

  vtester.vector_test1(fei_vec);

  test_Matrix mtester(MPI_COMM_WORLD);

  fei::SharedPtr<fei::Matrix> fei_mat = mtester.create_matrix(factory);

  mtester.matrix_test1(fei_mat);

#ifndef FEI_SER
  if (MPI_Finalize() != MPI_SUCCESS) ERReturn(-1);
#endif

  return(0);
}

