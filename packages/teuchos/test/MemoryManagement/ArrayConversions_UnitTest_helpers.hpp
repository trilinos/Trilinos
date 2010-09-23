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
