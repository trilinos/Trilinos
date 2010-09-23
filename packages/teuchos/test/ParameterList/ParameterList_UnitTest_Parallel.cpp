// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
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

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace Teuchos {


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, xmlUpdateAndBroadcast ) {
  const RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  // Test the broadcast functionality to avoid unscalable I/O collisions
  std::string inputFile="input.xml";
  ParameterList A;
  ParameterList B;
  updateParametersFromXmlFile(inputFile, &A);
  updateParametersFromXmlFileAndBroadcast(inputFile, &B, *comm);
  out << "B = " << B;
  TEST_ASSERT( B.begin() != B.end() ); // Avoid false positive from empty lists

  // See if any process returned a failed (i.e. a non-zero local_failed)
  int local_failed = !(A == B);
  int global_failed = -1;
  reduceAll( *comm, Teuchos::REDUCE_SUM, local_failed, outArg(global_failed) );
  TEST_EQUALITY_CONST( global_failed, 0 );
}


} // namespace Teuchos



