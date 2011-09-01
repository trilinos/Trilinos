// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "RTOpPack_version.hpp"
#include "RTOpPack_LapackWrappers.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ScalarTraits.hpp"


//
// Local testing support
//

#include "Teuchos_LocalTestingHelpers.hpp"


// Overridde some of the macros

#undef TEST_ARRAY_ELE_EQUALITY
#define TEST_ARRAY_ELE_EQUALITY( a, i, val ) \
   TEUCHOS_TEST_ARRAY_ELE_EQUALITY( a, i, val, true, out, success )

#define TEST_MATRIX_ELE_EQUALITY( a, i, j, val ) \
  TEUCHOS_TEST_MATRIX_ELE_EQUALITY( a, i, j, val, true, out, success )

#define TEST_MATRIX_ELE_FLOATING_EQUALITY( a, i, j, val, tol ) \
  TEUCHOS_TEST_MATRIX_ELE_FLOATING_EQUALITY( a, i, j, val, tol, true, out, success )


//
// Testing functions
//


template<class Scalar>
bool testLapackWrappers(
  const int n,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType &tol,
  Teuchos::FancyOStream &out 
  )
{

  using Teuchos::outArg;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::OSTab;
  using Teuchos::Array;
  using RTOpPack::SubMultiVectorView;

  typedef Teuchos_Index Ordinal;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  bool success = true;

  out << "\n" << RTOpPack::version() << "\n";

  out
    << "\n***"
    << "\n*** Testing RTOp wrappers to LAPACK " << ST::name()
    << "\n***\n";

  //
  out << "\nA) Construct a simple "<<n<<"-by-"<<n<<" random matrix M ...\n\n";
  //

  SubMultiVectorView<Scalar> M(n, n);
  ST::seedrandom(0);
  for ( int i = 0; i < n; ++i ) {
    for ( int j = 0; j < n; ++j ) {
      M(i,j) = ST::random();
    }
  }

  //
  out << "\nB) Create the factorization M = P * L * U ...\n\n";
  //
  
  SubMultiVectorView<Scalar> LU(n, n);
  Array<int> ipiv(n);
  
  {

    OSTab tab(out);

    RTOpPack::assign_entries<Scalar>( outArg(LU), M );

    int rank;
    RTOpPack::getrf<Scalar>( LU, ipiv, outArg(rank) );

    // Note: The contents of ipiv and LU is technically an implementation
    // detail of LAPACK and you should not try to directly test its contents!
    
    TEST_EQUALITY( rank, n );

  }

  //
  out << "\nC) Perform backsolves with the L and U factors and check result explicitly ...\n\n";
  //

  {

    OSTab tab(out);

    SubMultiVectorView<Scalar> X_known(n,1);
    ST::seedrandom(0);
    for ( int i = 0; i < n; ++i ) {
      X_known(i,0) = ST::random();
    }

    SubMultiVectorView<Scalar> B(n,1);
    SubMultiVectorView<Scalar> X(n,1);

    //
    out << "\nC.1) Performing non-transposed solve ...\n\n";
    //

    {

      OSTab tab2(out);
    
      for ( int i = 0; i < n; ++i ) {
        B(i,0) = 0.0;
        for ( int j = 0; j < n; ++j ) {
          B(i,0) += M(i,j) * X_known(j,0);
        }
      }
    
      RTOpPack::assign_entries<Scalar>( outArg(X), B );
      RTOpPack::getrs<Scalar>( LU, ipiv, RTOpPack::NOTRANS, outArg(X) );
      
      for ( int i = 0; i < n; ++i ) {
        TEST_MATRIX_ELE_FLOATING_EQUALITY( X, i, 0, X_known(i,0), tol );
      }

    }

    //
    out << "\nC.2) Performing non-conjugate transposed solve ...\n\n";
    //

    {

      OSTab tab2(out);
    
      for ( int i = 0; i < n; ++i ) {
        B(i,0) = 0.0;
        for ( int j = 0; j < n; ++j ) {
          B(i,0) += M(j,i) * X_known(j,0);
        }
      }
    
      RTOpPack::assign_entries<Scalar>( outArg(X), B );
      RTOpPack::getrs<Scalar>( LU, ipiv, RTOpPack::TRANS, outArg(X) );
      
      for ( int i = 0; i < n; ++i ) {
        TEST_MATRIX_ELE_FLOATING_EQUALITY( X, i, 0, X_known(i,0), tol );
      }

    }

    //
    out << "\nC.3) Performing conjugate transposed solve ...\n\n";
    //

    {

      OSTab tab2(out);
    
      for ( int i = 0; i < n; ++i ) {
        B(i,0) = 0.0;
        for ( int j = 0; j < n; ++j ) {
          B(i,0) += ST::conjugate(M(j,i)) * X_known(j,0);
        }
      }
    
      RTOpPack::assign_entries<Scalar>( outArg(X), B );
      RTOpPack::getrs<Scalar>( LU, ipiv, RTOpPack::CONJTRANS, outArg(X) );
      
      for ( int i = 0; i < n; ++i ) {
        TEST_MATRIX_ELE_FLOATING_EQUALITY( X, i, 0, X_known(i,0), tol );
      }

    }

  }
    
  return success;

}

int main(int argc, char* argv[])
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::OSTab;
  using Teuchos::CommandLineProcessor;
  using Teuchos::ScalarTraits;

  bool success = true, result;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  try {

		CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    int n = 4;
    clp.setOption( "n", &n, "Dimension of the system." );

    double epsScale = 200.0;
    clp.setOption( "eps-scale", &epsScale, "Constant (greater than 1) to scale eps by in error tests." );
    
		CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
      return parse_return;

    RCP<FancyOStream>
      out = VerboseObjectBase::getDefaultOStream();

#if defined(HAVE_TEUCHOS_BLASFLOAT) && defined(HAVE_TEUCHOS_FLOAT)
    result = testLapackWrappers<float>(n, epsScale*ScalarTraits<float>::eps(), *out);
    if(!result) success = false;
#endif

    result = testLapackWrappers<double>(n, epsScale*ScalarTraits<double>::eps(), *out);
    if(!result) success = false;

#if defined(HAVE_COMPLEX) && defined(HAVE_TEUCHOS_COMPLEX)
    result = testLapackWrappers<std::complex<double> >(
      n, epsScale*ScalarTraits<std::complex<double> >::eps(), *out);
    if(!result) success = false;
#endif
    
    if(success)
      *out << "\nEnd Result: TEST PASSED\n";
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);
    
  return ( success ? 0 : 1 );
  
}
