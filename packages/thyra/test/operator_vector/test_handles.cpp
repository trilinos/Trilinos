// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#include "Teuchos_MPISession.hpp"
#include "Thyra_VectorImpl.hpp"
#include "Thyra_VectorSpaceImpl.hpp"
#include "Thyra_VectorOpTester.hpp"
#include "Thyra_DefaultSerialVectorSpace.hpp"
#include "Teuchos_ScalarTraits.hpp"

using namespace Teuchos;
using namespace Thyra;

int main(int argc, void *argv[]) 
{
  try
    {
      MPISession::init(&argc, &argv);
      typedef Teuchos::ScalarTraits<double> ST;

      int n = 4;

      RefCountPtr<const VectorSpaceBase<double> > sp 
        = rcp(new DefaultSerialVectorSpace<double>(n));
      VectorSpace<double> space = sp;

      Vector<double> x = space.createMember();

      Thyra::randomize((double)(-ST::one()), (double)(+ST::one()), 
                       x.ptr().get());

      for (int i=0; i<n; i++)
        {
          double x_i = x[i];
          cerr << "i=" << i << " x=" << x_i << endl;
        }

      //  VectorOpTester<double> tester(space, TestSpecifier<double>(true, 1.0e-13, 1.0e-10));

      //      tester.runAllTests();
      
    }
  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
    }
  MPISession::finalize();

}

