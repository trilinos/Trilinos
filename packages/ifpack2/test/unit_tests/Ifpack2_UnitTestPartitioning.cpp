// @HEADER
// ***********************************************************************
// 
//                IFPACK
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#ifdef HAVE_IFPACK2_QD
#include <qd/dd_real.h>
#endif

#include <Ifpack2_UnitTestHelpers.hpp>
#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#else
#include <Teuchos_DefaultSerialComm.hpp>
#endif
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Ifpack2_OverlappingPartitioner.hpp>
#include <Ifpack2_LinearPartitioner.hpp>





using Tpetra::global_size_t;
typedef tif_utest::Node Node;
using namespace std;
using Teuchos::rcp;
using Teuchos::RCP;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Partitioning, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{  
  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CRS;
  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> ROW;
  
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> CRSG;

  // Useful matrices and such (tridiagonal test)
  global_size_t num_rows_per_proc = 5;
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc); 
  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Matrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  // ====================================== //
  //            point blocking              //
  // ====================================== //
  {
    Teuchos::ParameterList List;
    List.set("partitioner: local parts",(int) num_rows_per_proc);
    Ifpack2::LinearPartitioner<CRSG > MyPart(Matrix->getGraph());
    MyPart.setParameters(List);
    MyPart.compute();
    const Teuchos::ArrayView<const LocalOrdinal>  & myview = MyPart.nonOverlappingPartition();
    
    Teuchos::Array<LocalOrdinal> correct_solution((int)num_rows_per_proc);
    for(int i=0;i<(int)num_rows_per_proc;i++)
      correct_solution[i]=i;

    TEST_COMPARE_ARRAYS(myview,correct_solution);
  }

  // ====================================== //
  //            full blocking               //
  // ====================================== //
  // Point blocking
  {
    Teuchos::ParameterList List;
    List.set("partitioner: local parts",1);
    Ifpack2::LinearPartitioner<CRSG > MyPart(Matrix->getGraph());
    MyPart.setParameters(List);
    MyPart.compute();
    const Teuchos::ArrayView<const LocalOrdinal>  & myview = MyPart.nonOverlappingPartition();
    
    Teuchos::Array<LocalOrdinal> correct_solution((int)num_rows_per_proc);
    for(int i=0;i<(int)num_rows_per_proc;i++)
      correct_solution[i]=0;

    TEST_COMPARE_ARRAYS(myview,correct_solution);
  }


}



#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Partitioning, Test0, Scalar, LocalOrdinal,GlobalOrdinal)

UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)
#ifndef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
UNIT_TEST_GROUP_SCALAR_ORDINAL(float, short, int)
#endif
