//@HEADER
// ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
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
//@HEADER



#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Thyra_VectorImpl.hpp"
#include "Thyra_VectorSpaceImpl.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_LinearOperatorImpl.hpp"
#include "Thyra_LinearCombinationImpl.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"


using namespace Teuchos;
using namespace Thyra;


LinearOperator<double> makeRandomDenseOperator(int nc, const VectorSpace<double>& rowSp)
{
  RefCountPtr<Thyra::MultiVectorBase<double> > mv = rowSp.createMembers(nc);
  Thyra::randomize(-1.0, 1.0, &*mv);
  RefCountPtr<Thyra::LinearOpBase<double> > rtn = mv;
  return rtn;
}

int main(int argc, char *argv[]) 
{
  bool success = false;
  
  GlobalMPISession mpiSession(&argc, &argv);
  typedef Teuchos::ScalarTraits<double> ST;
  
  // Get stream that can print to just root or all streams!
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  
  try
    {
      Array<int> rangeSpaceSizes = tuple(5, 10);
      Array<int> domainSpaceSizes = tuple(5, 8);

      Array<VectorSpace<double> > rangeBlocks(rangeSpaceSizes.size());

      Array<Array<LinearOperator<double> > > blocks(rangeSpaceSizes.size());

      for (unsigned int br=0; br<rangeSpaceSizes.size(); br++)
        {
          int n = rangeSpaceSizes[br];
          rangeBlocks[br] 
            = new DefaultSpmdVectorSpace<double>(DefaultComm<Index>::getComm(),n,-1);

          blocks[br].resize(domainSpaceSizes.size());

          for (unsigned int bc=0; bc<domainSpaceSizes.size(); bc++)
            {
              blocks[br][bc] = makeRandomDenseOperator(domainSpaceSizes[bc],
                                                       rangeBlocks[br]);
            }
        }
      
      
      LinearOperator<double> A = block2x2(blocks[0][0], blocks[0][1],
                                          blocks[1][0], blocks[1][1]);

      VectorSpace<double> domain = A.domain();
      VectorSpace<double> range = A.range();

      *out << "A num block rows = " << A.numBlockRows() << endl;
      *out << "A num block cols = " << A.numBlockCols() << endl;

      *out << "A domain size = " << domain.dim() << endl;
      *out << "A range size = " << range.dim() << endl;

      Vector<double> x = domain.createMember();
      cerr << "randomizing trial vector" << endl;
      Thyra::randomize(-1.0, 1.0, x.ptr().get());

      Array<Vector<double> > xBlock(domain.numBlocks());
      for (unsigned int i=0; i<xBlock.size(); i++)
        {
          xBlock[i] = x.getBlock(i);
        }

      *out << "x size = " << space(x).dim() << endl;

      cerr << "------------------------------------------------------------" << endl;
      cerr << "computing A*x..." << endl;
      Vector<double> y0;
      y0 = A * x;
      cout << "y0 = " << y0.ptr().get() << endl;
      for (int i=0; i<space(y0).numBlocks(); i++)
        {
          cerr << "y0[" << i << "] = " << endl << y0.getBlock(i) << endl;
        }
      

      Vector<double> y1 = range.createMember();
      cerr << "------------------------------------------------------------" << endl;
      cerr << "computing A*x block-by-block..." << endl;
      Array<Vector<double> > yBlock(range.numBlocks());
      for (unsigned int i=0; i<yBlock.size(); i++)
        {
          yBlock[i] = range.getBlock(i).createMember();
          zeroOut(yBlock[i]);
          for (unsigned int j=0; j<xBlock.size(); j++)
            {
              LinearOperator<double> Aij = A.getBlock(i,j);
              if (Aij.ptr().get() != 0)
                {
                  cerr << "A(" << i << ", " << j << ") = " << endl 
                       << Aij << endl;
                }
              else
                {
                  cerr << "A(" << i << ", " << j << ") = 0 " << endl;
                }
              cerr << "x[" << j << "] = " << endl << xBlock[j] << endl;
              if (Aij.ptr().get()==0) continue;
              yBlock[i] = yBlock[i] + Aij * xBlock[j];
            }
          y1.setBlock(i, yBlock[i]);
        }

      for (int i=0; i<space(y1).numBlocks(); i++)
        {
          cerr << "y1[" << i << "] = " << endl << y1.getBlock(i) << endl;
        }
      double err = norm2(y1 - y0);
      cerr << "error = " << err << endl;

      double tol = 1.0e-13;
      if (err < tol)
        {
          cerr << "block op test PASSED" << endl;
        }
      else
        {
          cerr << "block op test FAILED" << endl;
        }
    }
  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
    }
}



