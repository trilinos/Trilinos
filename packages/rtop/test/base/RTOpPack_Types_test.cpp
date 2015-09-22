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

#include "RTOpPack_Types.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"


//
// Define local macros to make defining tests easier for this particular test
// code.
//
// Note, macros with these types of names should only exist in a *.cpp file
// after all #includes are done!
//


#define TEST_EQUALITY_CONST( v1, v2 ) \
  TEUCHOS_TEST_EQUALITY_CONST( v1, v2, out, success )

#define TEST_EQUALITY( v1, v2 ) \
  TEUCHOS_TEST_EQUALITY( v1, v2, out, success )

#define TEST_ARRAY_ELE_EQUALITY( a, i, val ) \
   TEUCHOS_TEST_ARRAY_ELE_EQUALITY( a, i, val, false, out, local_success )

#define TEST_THROW( code, ExceptType  ) \
  TEUCHOS_TEST_THROW( code, ExceptType, out, success  )

#define TEST_NOTHROW( code  ) \
  TEUCHOS_TEST_NOTHROW( code, out, success  )


//
// Templated unit testing functions
//


template<class Scalar>
bool testSubVectorView(
  const int subDim, const int stride, Teuchos::FancyOStream &out
  )
{
  
  using Teuchos::TypeNameTraits;
  using RTOpPack::ConstSubVectorView;
  using RTOpPack::SubVectorView;
  using Teuchos::as;
  using Teuchos::null;
  using Teuchos::ArrayRCP;
  using Teuchos::OSTab;

  bool success = true;
 
  out
    << "\n***"
    << "\n*** Testing "<<TypeNameTraits<SubVectorView<Scalar> >::name()
    <<" with subDim="<<subDim <<", stride="<<stride
    << "\n***\n";
  
  OSTab tab(out);

  {
    out << "\nTesting default construction for SubVectorView ...\n";
    SubVectorView<Scalar> sv;
    TEST_EQUALITY_CONST(sv.globalOffset(),0);
    TEST_EQUALITY_CONST(sv.subDim(),0);
    TEST_EQUALITY_CONST(sv.values(),null);
    TEST_EQUALITY_CONST(sv.stride(),0);
#ifdef TEUCHOS_DEBUG
      TEST_THROW(sv[0],std::out_of_range);
      TEST_THROW(sv(0),std::out_of_range);
#endif
  }

  {
    out << "\nTesting default construction for ConstSubVectorView ...\n";
    ConstSubVectorView<Scalar> sv;
    TEST_EQUALITY_CONST(sv.globalOffset(),0);
    TEST_EQUALITY_CONST(sv.subDim(),0);
    TEST_EQUALITY_CONST(sv.values(),null);
    TEST_EQUALITY_CONST(sv.stride(),0);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(sv[0],std::out_of_range);
    TEST_THROW(sv(0),std::out_of_range);
#endif
  }

  out << "\nCreating value data for subview ...\n";
  const ArrayRCP<Scalar>
    values = Teuchos::arcp<Scalar>(subDim*std::abs(stride));
  std::fill(values.begin(), values.end(), as<Scalar>(-1));
  {
    typename ArrayRCP<Scalar>::iterator itr = values.begin();
    if (stride < 0)
      itr += ( (-subDim*stride) - 1 );
    for ( int i = 0; i < subDim; ++i ) {
      *itr = as<Scalar>(i);
      itr+=stride;
    }
  }

  {
    out << "\nCreating SubVectorView for serial ...\n";
    SubVectorView<Scalar> sv(values);
    TEST_EQUALITY_CONST(sv.globalOffset(), 0);
    TEST_EQUALITY(sv.subDim(), as<int>(values.size()));
    TEST_EQUALITY(sv.values(), values);
    TEST_EQUALITY_CONST(sv.stride(), 1);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(sv[-1],std::out_of_range);
    TEST_THROW(sv(-1),std::out_of_range);
    TEST_THROW(sv[values.size()],std::out_of_range);
    TEST_THROW(sv(values.size()),std::out_of_range);
#endif
  }

  out << "\nCreating SubVectorView ...\n";
  const Teuchos_Ordinal globalOffset = 20;
  SubVectorView<Scalar> sv(globalOffset, subDim, values, stride);

  {
    out << "\nTesting SubVectorView access functions ...\n";
    TEST_EQUALITY(sv.globalOffset(),globalOffset);
    TEST_EQUALITY(sv.subDim(),subDim);
    TEST_EQUALITY(sv.values(),values);
    TEST_EQUALITY(sv.stride(),stride);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(sv[-1],std::out_of_range);
    TEST_THROW(sv(-1),std::out_of_range);
    TEST_THROW(sv[subDim],std::out_of_range);
    TEST_THROW(sv(subDim),std::out_of_range);
#endif
    out << "\nTest that sv[i] == i ... ";
    bool local_success = true;
    for ( int i = 0; i < subDim; ++i )
      TEST_ARRAY_ELE_EQUALITY(sv,i,as<Scalar>(i));
    if (local_success) out << "passed\n";
    else success = false;
  }
  
  out << "\nCreating ConstSubVectorView ...\n";
  SubVectorView<Scalar> csv(globalOffset,subDim,values,stride);

  {
    out << "\nTesting ConstSubVectorView access functions ...\n";
    TEST_EQUALITY(csv.globalOffset(),globalOffset);
    TEST_EQUALITY(csv.subDim(),subDim);
    TEST_EQUALITY(csv.values(),values);
    TEST_EQUALITY(csv.stride(),stride);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(csv[-1],std::out_of_range);
    TEST_THROW(csv(-1),std::out_of_range);
    TEST_THROW(csv[subDim],std::out_of_range);
    TEST_THROW(csv(subDim),std::out_of_range);
#endif
    out << "\nTest that csv[i] == i ... ";
    bool local_success = true;
    for ( int i = 0; i < subDim; ++i )
      TEST_ARRAY_ELE_EQUALITY(csv,i,as<Scalar>(i));
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    out << "\nUninitialize SubVectorView ...\n";
    sv.uninitialize();
    TEST_EQUALITY_CONST(sv.globalOffset(),0);
    TEST_EQUALITY_CONST(sv.subDim(),0);
    TEST_EQUALITY_CONST(sv.values(),null);
    TEST_EQUALITY_CONST(sv.stride(),0);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(sv[0],std::out_of_range);
    TEST_THROW(sv(0),std::out_of_range);
#endif
  }
  
  {
    out << "\nUninitialize ConstSubVectorView ...\n";
    csv.uninitialize();
    TEST_EQUALITY_CONST(csv.globalOffset(),0);
    TEST_EQUALITY_CONST(csv.subDim(),0);
    TEST_EQUALITY_CONST(csv.values(),null);
    TEST_EQUALITY_CONST(csv.stride(),0);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(csv[0],std::out_of_range);
    TEST_THROW(csv(0),std::out_of_range);
#endif
  }

  {
    out << "\nInitialize SubVectorView ...\n";
    sv.initialize(globalOffset,subDim,values,stride);
    TEST_EQUALITY(sv.globalOffset(),globalOffset);
    TEST_EQUALITY(sv.subDim(),subDim);
    TEST_EQUALITY(sv.values(),values);
    TEST_EQUALITY(sv.stride(),stride);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(sv[-1],std::out_of_range);
    TEST_THROW(sv(-1),std::out_of_range);
    TEST_THROW(sv[subDim],std::out_of_range);
    TEST_THROW(sv(subDim),std::out_of_range);
#endif
  }

  {
    out << "\nInitialize ConstSubVectorView ...\n";
    csv.initialize(globalOffset, subDim, values, stride);
    TEST_EQUALITY(csv.globalOffset(), globalOffset);
    TEST_EQUALITY(csv.subDim(), subDim);
    TEST_EQUALITY(csv.values(), values);
    TEST_EQUALITY(csv.stride(), stride);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(csv[-1], std::out_of_range);
    TEST_THROW(csv(-1), std::out_of_range);
    TEST_THROW(csv[subDim], std::out_of_range);
    TEST_THROW(csv(subDim), std::out_of_range);
#endif
  }

  {
    out << "\nCopy construct SubVectorView ...\n";
    const SubVectorView<Scalar> sv2(sv);
    TEST_EQUALITY(sv2.globalOffset(),globalOffset);
    TEST_EQUALITY(sv2.subDim(),subDim);
    TEST_EQUALITY(sv2.values(),values);
    TEST_EQUALITY(sv2.stride(),stride);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(sv2[-1],std::out_of_range);
    TEST_THROW(sv2(-1),std::out_of_range);
    TEST_THROW(sv2[subDim],std::out_of_range);
    TEST_THROW(sv2(subDim),std::out_of_range);
#endif
  }

  {
    out << "\nCopy construct ConstSubVectorView ...\n";
    const ConstSubVectorView<Scalar> csv2(csv);
    TEST_EQUALITY(csv2.globalOffset(),globalOffset);
    TEST_EQUALITY(csv2.subDim(),subDim);
    TEST_EQUALITY(csv2.values(),values);
    TEST_EQUALITY(csv2.stride(),stride);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(csv2[-1],std::out_of_range);
    TEST_THROW(csv2(-1),std::out_of_range);
    TEST_THROW(csv2[subDim],std::out_of_range);
    TEST_THROW(csv2(subDim),std::out_of_range);
#endif
  }

  {
    out << "\nTesting Change of SubVectorView through element access  ...\n";
    for ( int i = 0; i < subDim; ++i )
      sv[i] = as<Scalar>(subDim-i-1);
    out << "\nTest that sv[i] == subDim-i-1 ... ";
    bool local_success = true;
    for ( int i = 0; i < subDim; ++i )
      TEST_ARRAY_ELE_EQUALITY(sv,i,as<Scalar>(subDim-i-1));
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    out << "\nTest copy of elements ...\n";
    SubVectorView<Scalar> sv2(subDim);
    RTOpPack::assign_entries<Scalar>(Teuchos::outArg(sv2), csv);
    out << "\nTest that sv2[i] == subDim-i-1 ... ";
    bool local_success = true;
    for ( int i = 0; i < subDim; ++i )
      TEST_ARRAY_ELE_EQUALITY(sv2,i,as<Scalar>(subDim-i-1));
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    out << "\nTest change of globalOffset ...\n";
    const Teuchos_Ordinal newGlobalOffset = globalOffset + 1;
    csv.setGlobalOffset(newGlobalOffset);
    TEST_EQUALITY(csv.globalOffset(),newGlobalOffset);
  }

  out << "\nCreating empty SubVectorView ...\n";
  SubVectorView<Scalar> sv_0(globalOffset, 0, null, stride);

  {
    out << "\nTesting SubVectorView access functions ...\n";
    TEST_EQUALITY_CONST(sv_0.globalOffset(), globalOffset);
    TEST_EQUALITY_CONST(sv_0.subDim(), 0);
    TEST_EQUALITY_CONST(sv_0.values(), null);
    TEST_EQUALITY_CONST(sv_0.stride(), stride);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(sv_0[0], std::out_of_range);
    TEST_THROW(sv_0(0), std::out_of_range);
#endif
  }
  
  out << "\nCreating empty ConstSubVectorView ...\n";
  SubVectorView<Scalar> csv_0(globalOffset, 0, null, stride);

  {
    out << "\nTesting ConstSubVectorView access functions ...\n";
    TEST_EQUALITY_CONST(csv_0.globalOffset(), globalOffset);
    TEST_EQUALITY_CONST(csv_0.subDim(), 0);
    TEST_EQUALITY_CONST(csv_0.values(), null);
    TEST_EQUALITY_CONST(csv_0.stride(), stride);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(csv_0[0], std::out_of_range);
    TEST_THROW(csv_0(0), std::out_of_range);
#endif
  }
  
  return success;

}


template<class Scalar>
bool testSubMultiVectorView(
  const int subDim, const int numSubCols, const int leadingDim,
  Teuchos::FancyOStream &out
  )
{
  
  using RTOpPack::ConstSubMultiVectorView;
  using RTOpPack::SubMultiVectorView;
  using RTOpPack::ConstSubVectorView;
  using RTOpPack::SubVectorView;
  using Teuchos::TypeNameTraits;
  using Teuchos::ArrayRCP;
  using Teuchos::Array;
  using Teuchos::OSTab;
  using Teuchos::null;
  using Teuchos::as;

  bool success = true;
 
  out
    << "\n***"
    << "\n*** Testing "<<TypeNameTraits<SubMultiVectorView<Scalar> >::name()
    <<" with subDim="<<subDim <<", numSubCols="<<numSubCols<<", leadingDim="<<leadingDim
    << "\n***\n";
  
  OSTab tab(out);

  {
    out << "\nTesting default construction for SubMultiVectorView ...\n";
    SubMultiVectorView<Scalar> smv;
    TEST_EQUALITY_CONST(smv.globalOffset(),0);
    TEST_EQUALITY_CONST(smv.subDim(),0);
    TEST_EQUALITY_CONST(smv.colOffset(),0);
    TEST_EQUALITY_CONST(smv.numSubCols(),0);
    TEST_EQUALITY_CONST(smv.values(),null);
    TEST_EQUALITY_CONST(smv.leadingDim(),0);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(smv(0,0),std::out_of_range);
    TEST_THROW(smv.col(0),std::out_of_range);
#endif
  }

  {
    out << "\nTesting default construction for ConstSubMultiVectorView ...\n";
    ConstSubMultiVectorView<Scalar> csmv;
    TEST_EQUALITY_CONST(csmv.globalOffset(),0);
    TEST_EQUALITY_CONST(csmv.subDim(),0);
    TEST_EQUALITY_CONST(csmv.colOffset(),0);
    TEST_EQUALITY_CONST(csmv.numSubCols(),0);
    TEST_EQUALITY_CONST(csmv.values(),null);
    TEST_EQUALITY_CONST(csmv.leadingDim(),0);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(csmv(0,0),std::out_of_range);
    TEST_THROW(csmv.col(0),std::out_of_range);
#endif
  }

  out << "\nCreating value data for subview ...\n";

  const ArrayRCP<Scalar>
    values = Teuchos::arcp<Scalar>(numSubCols*leadingDim);
  std::fill(values.begin(),values.end(),as<Scalar>(-1));

  // 2008/03/04: rabartl: Here I am creating a simple 2D array to hold the
  // expected values for the SubMultiVector objects.  This is to avoid strange
  // roundoff problems that seem to be happening on some platforms
  // (e.g. beowulf) for some reason.  By storing the expected values instead
  // of recomputing them, this should compare correctly on all platforms.
  Array<Array<Scalar> > contigValues(subDim);
  for ( int i = 0; i < subDim; ++i ) contigValues[i].resize(numSubCols); 

  typename ArrayRCP<Scalar>::iterator itr = values.begin();
  for ( int j = 0; j < numSubCols; ++j ) {
    typename ArrayRCP<Scalar>::iterator col_itr = itr + j*leadingDim;
    for ( int i = 0; i < subDim; ++i ) {
      contigValues[i][j] = as<Scalar>(i) + as<Scalar>(j)/as<Scalar>(1000);
      *col_itr++ = contigValues[i][j];
    }
  }

  out << "\nCreating SubMultiVectorView ...\n";
  const Teuchos_Ordinal globalOffset = 20;
  const Teuchos_Ordinal colOffset = 15;
  SubMultiVectorView<Scalar>
    smv(globalOffset, subDim, colOffset, numSubCols, values, leadingDim);

  {
    out << "\nTesting SubMultiVectorView access functions ...\n";
    TEST_EQUALITY(smv.globalOffset(),globalOffset);
    TEST_EQUALITY(smv.subDim(),subDim);
    TEST_EQUALITY(smv.values(),values);
    TEST_EQUALITY(smv.colOffset(),colOffset);
    TEST_EQUALITY(smv.numSubCols(),numSubCols);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(smv(-1,0),std::out_of_range);
    TEST_THROW(smv(subDim,0),std::out_of_range);
    TEST_THROW(smv(0,-1),std::out_of_range);
    TEST_THROW(smv(0,numSubCols),std::out_of_range);
#endif
    out << "\nTest that smv(i,j) == i + j/1000 ...\n";
    for ( int j = 0; j < numSubCols; ++j ) {
      for ( int i = 0; i < subDim; ++i ) {
        TEST_EQUALITY(smv(i,j), contigValues[i][j]);
      }
    }
    out << "\nTest that smv.col(j)(i) == i + j/1000 ...\n";
    for ( int j = 0; j < numSubCols; ++j ) {
      const SubVectorView<Scalar> smv_col_j = smv.col(j);
      for ( int i = 0; i < subDim; ++i ) {
        TEST_EQUALITY(smv_col_j[i], contigValues[i][j]);
      }
    }
  }

  out << "\nCreating ConstSubMultiVectorView ...\n";
  ConstSubMultiVectorView<Scalar>
    csmv(globalOffset, subDim, colOffset, numSubCols, values, leadingDim);

  {
    out << "\nTesting ConstSubMultiVectorView access functions ...\n";
    TEST_EQUALITY(csmv.globalOffset(),globalOffset);
    TEST_EQUALITY(csmv.subDim(),subDim);
    TEST_EQUALITY(csmv.values(),values);
    TEST_EQUALITY(csmv.colOffset(),colOffset);
    TEST_EQUALITY(csmv.numSubCols(),numSubCols);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(csmv(-1,0),std::out_of_range);
    TEST_THROW(csmv(subDim,0),std::out_of_range);
    TEST_THROW(csmv(0,-1),std::out_of_range);
    TEST_THROW(csmv(0,numSubCols),std::out_of_range);
#endif
    out << "\nTest that csmv(i,j) == i + j/1000 ...\n";
    for ( int j = 0; j < numSubCols; ++j ) {
      for ( int i = 0; i < subDim; ++i ) {
        TEST_EQUALITY(csmv(i,j), contigValues[i][j]);
      }
    }
    out << "\nTest that csmv.col(j)(i) == i + j/1000 ...\n";
    for ( int j = 0; j < numSubCols; ++j ) {
      const ConstSubVectorView<Scalar> csmv_col_j = csmv.col(j);
      for ( int i = 0; i < subDim; ++i ) {
        TEST_EQUALITY(csmv_col_j[i], contigValues[i][j]);
      }
    }
  }

  {
    out << "\nUninitialize SubMultiVectorView ...\n";
    smv.uninitialize();
    TEST_EQUALITY_CONST(smv.globalOffset(),0);
    TEST_EQUALITY_CONST(smv.subDim(),0);
    TEST_EQUALITY_CONST(smv.colOffset(),0);
    TEST_EQUALITY_CONST(smv.numSubCols(),0);
    TEST_EQUALITY_CONST(smv.values(),null);
    TEST_EQUALITY_CONST(smv.leadingDim(),0);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(smv(0,0),std::out_of_range);
    TEST_THROW(smv.col(0),std::out_of_range);
#endif
  }

  {
    out << "\nUninitialize ConstSubMultiVectorView ...\n";
    csmv.uninitialize();
    TEST_EQUALITY_CONST(csmv.globalOffset(),0);
    TEST_EQUALITY_CONST(csmv.subDim(),0);
    TEST_EQUALITY_CONST(csmv.colOffset(),0);
    TEST_EQUALITY_CONST(csmv.numSubCols(),0);
    TEST_EQUALITY_CONST(csmv.values(),null);
    TEST_EQUALITY_CONST(csmv.leadingDim(),0);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(csmv(0,0),std::out_of_range);
    TEST_THROW(csmv.col(0),std::out_of_range);
#endif
  }

  {
    out << "\nInitialize SubMultiVectorView ...\n";
    smv.initialize(globalOffset, subDim, colOffset, numSubCols, values, leadingDim);
    TEST_EQUALITY(smv.globalOffset(),globalOffset);
    TEST_EQUALITY(smv.subDim(),subDim);
    TEST_EQUALITY(smv.values(),values);
    TEST_EQUALITY(smv.colOffset(),colOffset);
    TEST_EQUALITY(smv.numSubCols(),numSubCols);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(smv(-1,0),std::out_of_range);
    TEST_THROW(smv(subDim,0),std::out_of_range);
    TEST_THROW(smv(0,-1),std::out_of_range);
    TEST_THROW(smv(0,numSubCols),std::out_of_range);
#endif
  }

  {
    out << "\nInitialize ConstSubMultiVectorView ...\n";
    csmv.initialize(globalOffset, subDim, colOffset, numSubCols, values, leadingDim);
    TEST_EQUALITY(csmv.globalOffset(),globalOffset);
    TEST_EQUALITY(csmv.subDim(),subDim);
    TEST_EQUALITY(csmv.values(),values);
    TEST_EQUALITY(csmv.colOffset(),colOffset);
    TEST_EQUALITY(csmv.numSubCols(),numSubCols);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(csmv(-1,0),std::out_of_range);
    TEST_THROW(csmv(subDim,0),std::out_of_range);
    TEST_THROW(csmv(0,-1),std::out_of_range);
    TEST_THROW(csmv(0,numSubCols),std::out_of_range);
#endif
  }

  {
    out << "\nCopy constructor for SubMultiVectorView ...\n";
    const SubMultiVectorView<Scalar> smv2(smv);
    TEST_EQUALITY(smv2.globalOffset(),globalOffset);
    TEST_EQUALITY(smv2.subDim(),subDim);
    TEST_EQUALITY(smv2.values(),values);
    TEST_EQUALITY(smv2.colOffset(),colOffset);
    TEST_EQUALITY(smv2.numSubCols(),numSubCols);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(smv2(-1,0),std::out_of_range);
    TEST_THROW(smv2(subDim,0),std::out_of_range);
    TEST_THROW(smv2(0,-1),std::out_of_range);
    TEST_THROW(smv2(0,numSubCols),std::out_of_range);
#endif
  }

  {
    out << "\nCopy constructor for ConstSubMultiVectorView ...\n";
    const ConstSubMultiVectorView<Scalar> csmv2(smv);
    TEST_EQUALITY(csmv2.globalOffset(),globalOffset);
    TEST_EQUALITY(csmv2.subDim(),subDim);
    TEST_EQUALITY(csmv2.values(),values);
    TEST_EQUALITY(csmv2.colOffset(),colOffset);
    TEST_EQUALITY(csmv2.numSubCols(),numSubCols);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(csmv2(-1,0),std::out_of_range);
    TEST_THROW(csmv2(subDim,0),std::out_of_range);
    TEST_THROW(csmv2(0,-1),std::out_of_range);
    TEST_THROW(csmv2(0,numSubCols),std::out_of_range);
#endif
  }

  {
    out << "\nChange of SubMultiVectorView to smv(i,j) = j + i/1000  ...\n";
    for ( int j = 0; j < numSubCols; ++j ) {
      for ( int i = 0; i < subDim; ++i ) {
        contigValues[i][j] = as<Scalar>(j) + as<Scalar>(i)/as<Scalar>(1000);
        smv(i,j) = contigValues[i][j];
      }
    }
    out << "\nTest that smv(i,j) == j + i/1000 ...\n";
    for ( int j = 0; j < numSubCols; ++j ) {
      for ( int i = 0; i < subDim; ++i ) {
        TEST_EQUALITY(smv(i,j), contigValues[i][j]);
      }
    }
  }

  {
    out << "\nTest copy of elements ...\n";
    SubMultiVectorView<Scalar> smv2(subDim, numSubCols);
    RTOpPack::assign_entries<Scalar>(Teuchos::outArg(smv2), csmv);
    out << "\nTest that smv2(i,j) == j + i/1000 ...\n";
    for ( int j = 0; j < numSubCols; ++j ) {
      for ( int i = 0; i < subDim; ++i ) {
        TEST_EQUALITY(smv2(i,j), contigValues[i][j]);
      }
    }
  }

  {
    out << "\nTest change of globalOffset ...\n";
    const Teuchos_Ordinal newGlobalOffset = globalOffset + 1;
    csmv.setGlobalOffset(newGlobalOffset);
    TEST_EQUALITY(csmv.globalOffset(),newGlobalOffset);
  }

  out << "\nCreating empty SubMultiVectorView ...\n";
  SubMultiVectorView<Scalar>
    smv_0(globalOffset, 0, colOffset, numSubCols, null, leadingDim);

  {
    out << "\nTesting SubMultiVectorView access functions ...\n";
    TEST_EQUALITY(smv_0.globalOffset(),globalOffset);
    TEST_EQUALITY(smv_0.subDim(), 0);
    TEST_EQUALITY(smv_0.values(), null);
    TEST_EQUALITY(smv_0.colOffset(), colOffset);
    TEST_EQUALITY(smv_0.numSubCols(), numSubCols);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(smv_0(0,0),std::out_of_range);
#endif
  }

  out << "\nCreating empty ConstSubMultiVectorView ...\n";
  ConstSubMultiVectorView<Scalar>
    csmv_0(globalOffset, 0, colOffset, numSubCols, null, leadingDim);

  {
    out << "\nTesting ConstSubMultiVectorView access functions ...\n";
    TEST_EQUALITY(csmv_0.globalOffset(), globalOffset);
    TEST_EQUALITY(csmv_0.subDim(), 0);
    TEST_EQUALITY(csmv_0.values(), null);
    TEST_EQUALITY(csmv_0.colOffset(), colOffset);
    TEST_EQUALITY(csmv_0.numSubCols(), numSubCols);
#ifdef TEUCHOS_DEBUG
    TEST_THROW(csmv_0(0,0),std::out_of_range);
#endif
  }
  
  return success;

}


//
// Main testing program
//

int main( int argc, char* argv[] ) {

  using Teuchos::CommandLineProcessor;
	
	bool success = true;
 
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  //const int procRank = Teuchos::GlobalMPISession::getRank();
 
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
 
	try {
    
    //
		// Read options from the commandline
    //

    CommandLineProcessor clp(false); // Don't throw exceptions

    int subDim = 4;
    clp.setOption( "sub-dim", &subDim,
      "Number of elements in the (multi)vector subview" );

    int stride = 1;
    clp.setOption( "stride", &stride,
      "Stride between elements" );

    int numCols = 3;
    clp.setOption( "num-cols", &numCols,
      "Number of columns in the multivector subview" );

    int leadingDim = 4;
    clp.setOption( "leading-dim", &leadingDim,
      "Leading dimension of the sub-multivector subviews" );

		CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv);

		if ( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) {
			*out << "\nEnd Result: TEST FAILED" << std::endl;
			return parse_return;
		}

    *out << std::endl << Teuchos::Teuchos_Version() << std::endl;

    {
      bool result;
      typedef double Scalar;
      result = testSubVectorView<Scalar>(subDim, stride, *out);
      if (!result) success = false;
      result = testSubMultiVectorView<Scalar>(subDim, numCols, leadingDim, *out);
      if (!result) success = false;
    }

#ifdef HAVE_TEUCHOS_FLOAT
    {
      bool result;
      typedef float Scalar;
      result = testSubVectorView<Scalar>(subDim, stride, *out);
      if (!result) success = false;
      result = testSubMultiVectorView<Scalar>(subDim, numCols, leadingDim, *out);
      if (!result) success = false;
    }
#endif
 
	}
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);
 
  if (success)
    *out << "\nEnd Result: TEST PASSED" << std::endl;
  else
    *out << "\nEnd Result: TEST FAILED" << std::endl;
 
  return ( success ? 0 : 1 );
 
}
