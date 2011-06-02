/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/


#ifndef ARRAY_CONVERSIONS_UNIT_TEST_HELPERS
#define ARRAY_CONVERSIONS_UNIT_TEST_HELPERS


namespace ArrayConversionsUnitTestHelpers {


using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Ptr;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::as;

extern Teuchos_Ordinal n;


template<class T>
Array<RCP<T> > generateArrayRcp(const Teuchos_Ordinal n_in)
{
  Array<RCP<T> > a(n_in);
  for (Teuchos_Ordinal i=0 ; i<n_in ; ++i) {
    RCP<T> data = rcp(new T(as<T>(i)));
    a[i] = data;
  }
  return a;
}


template<class T>
Array<RCP<T> > generateArrayRcpGen(const Teuchos_Ordinal n_in)
{
  Array<RCP<T> > a;
  for (Teuchos_Ordinal i=0 ; i<n_in ; ++i) {
    a.push_back(rcp(new T));
  }
  return a;
}


template<class T>
T testArrayViewInput(const ArrayView<const Ptr<const T> >& a_in)
{
  typedef Teuchos::ScalarTraits<T> ST;
  T a = ST::zero();
  for (Teuchos_Ordinal i=0; i<a_in.size(); ++i) {
    a += *a_in[i];
  }
  return a;
}


template<class T>
void testArrayViewOutput(const ArrayView<const Ptr<T> >& a_out)
{
  typedef Teuchos::ScalarTraits<T> ST;
  for (Teuchos_Ordinal i=0 ; i<a_out.size() ; ++i) {
    *a_out[i] = as<T>(i);
  }
}


} // namespace ArrayConversionsUnitTestHelpers 


#endif // ARRAY_CONVERSIONS_UNIT_TEST_HELPERS
