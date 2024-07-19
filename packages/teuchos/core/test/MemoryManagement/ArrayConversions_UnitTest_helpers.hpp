// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


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
  for (Teuchos_Ordinal i=0 ; i<a_out.size() ; ++i) {
    *a_out[i] = as<T>(i);
  }
}


} // namespace ArrayConversionsUnitTestHelpers


#endif // ARRAY_CONVERSIONS_UNIT_TEST_HELPERS
