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


template <class Scalar> inline
LinearOperator<Scalar> makeRandomDenseOperator(int nc, const VectorSpace<Scalar>& rowSp)
{
  typedef typename Teuchos::ScalarTraits<Scalar> ST;
  RCP<Thyra::MultiVectorBase<Scalar> > mv = rowSp.createMembers(nc);
  Thyra::randomize(-ST::one(), ST::one(), &*mv);
  RCP<Thyra::LinearOpBase<Scalar> > rtn = mv;
  return rtn;
}

template <class Scalar>
bool runTest(Teuchos::RCP<Teuchos::FancyOStream>& out);

int main(int argc, char *argv[]) 
{
  bool success = false;
  
  GlobalMPISession mpiSession(&argc, &argv);
  typedef Teuchos::ScalarTraits<double> ST;
  
  // Get stream that can print to just root or all streams!
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  
  try
    {
      CommandLineProcessor  clp;
      clp.throwExceptions(false);
      clp.addOutputSetupOptions(true);
      bool verbose = false;
      clp.setOption( "verbose", "quiet", &verbose, 
                     "Determines if any output is printed or not." );

      
      CommandLineProcessor::EParseCommandLineReturn parse_return 
        = clp.parse(argc,argv);
      if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

      if (!verbose) out = rcp(new FancyOStream(rcp(new oblackholestream())));


      success = runTest<double>(out);

      success = runTest<float>(out) && success;

#if defined(HAVE_TEUCHOS_COMPLEX)
      success = runTest<std::complex<double> >(out) && success;
#endif

    }

  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,out.get()?*out:std::cerr,success)

    if (success)
      {
        *out << "all tests PASSED!" << std::endl;
        return 0;
      }
    else
      {
        *out << "at least one test FAILED!" << std::endl;
        return 1;
      }
}






template <class Scalar> inline
bool runTest(Teuchos::RCP<Teuchos::FancyOStream>& out)
{
  typedef typename Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  *out << "==========================================================================="
       << std::endl;
  *out << "running block operator test for " << ST::name() << std::endl;
  *out << "==========================================================================="
       << std::endl;
  Array<int> rangeSpaceSizes = tuple(5, 10);
  Array<int> domainSpaceSizes = tuple(5, 8);
  
  Array<VectorSpace<Scalar> > rangeBlocks(rangeSpaceSizes.size());
  
  Array<Array<LinearOperator<Scalar> > > blocks(rangeSpaceSizes.size());
  
  for (unsigned int br=0; br<rangeSpaceSizes.size(); br++)
    {
      int n = rangeSpaceSizes[br];
      rangeBlocks[br] 
        = new DefaultSpmdVectorSpace<Scalar>(DefaultComm<Index>::getComm(),n,-1);
      
      blocks[br].resize(domainSpaceSizes.size());
      
      for (unsigned int bc=0; bc<domainSpaceSizes.size(); bc++)
        {
          blocks[br][bc] = makeRandomDenseOperator<Scalar>(domainSpaceSizes[bc],
                                                           rangeBlocks[br]);
        }
    }
  
  
  LinearOperator<Scalar> A = block2x2(blocks[0][0], blocks[0][1],
                                      blocks[1][0], blocks[1][1]);
  
  VectorSpace<Scalar> domain = A.domain();
  VectorSpace<Scalar> range = A.range();
  
  *out << "A num block rows = " << A.numBlockRows() << std::endl;
  *out << "A num block cols = " << A.numBlockCols() << std::endl;
  
  *out << "A domain size = " << domain.dim() << std::endl;
  *out << "A range size = " << range.dim() << std::endl;
  
  Vector<Scalar> x = domain.createMember();
  *out << "randomizing trial std::vector" << std::endl;
  Thyra::randomize(-ST::one(), ST::one(), x.ptr().get());
  
  Array<Vector<Scalar> > xBlock(domain.numBlocks());
  for (unsigned int i=0; i<xBlock.size(); i++)
    {
      xBlock[i] = x.getBlock(i);
    }
  
  *out << "x size = " << space(x).dim() << std::endl;
  
  *out << "------------------------------------------------------------" << std::endl;
  *out << "computing A*x..." << std::endl;
  Vector<Scalar> y0;
  y0 = A * x;
  
  
  Vector<Scalar> y1 = range.createMember();
  *out << "------------------------------------------------------------" << std::endl;
  *out << "computing A*x block-by-block..." << std::endl;
  Array<Vector<Scalar> > yBlock(range.numBlocks());
  for (unsigned int i=0; i<yBlock.size(); i++)
    {
      yBlock[i] = range.getBlock(i).createMember();
      zeroOut(yBlock[i]);
      for (unsigned int j=0; j<xBlock.size(); j++)
        {
          LinearOperator<Scalar> Aij = A.getBlock(i,j);
          if (Aij.ptr().get()==0) continue;
          yBlock[i] = yBlock[i] + Aij * xBlock[j];
        }
      y1.setBlock(i, yBlock[i]);
    }
  
  ScalarMag err = norm2(y1 - y0);
  *out << "error = " << err << std::endl;
  
  ScalarMag tol = 1.0e2 * ST::prec();

  bool ok = err < tol;
  if (ok)
    {
      *out << "err=" << err << " tol=" << tol << ", test PASSED!" << std::endl;
    }
  else
    {
      *out << "err=" << err << " tol=" << tol << ", test FAILED!" << std::endl;
    }

  return err < tol;
}

