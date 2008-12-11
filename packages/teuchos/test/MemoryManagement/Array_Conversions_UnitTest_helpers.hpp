
namespace ArrayConversionsUnitTestHelpers {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Ptr;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::as;

extern Teuchos_Ordinal n;

// Declarations:
template<class T>
T TestArrayViewInput(const ArrayView<const Ptr<const T> >& a_in);

template<class T>
void TestArrayViewOutput(const ArrayView<const Ptr<T> >& a_out);

template<class T>
Array<RCP<T> > generateArray(const Teuchos_Ordinal n);


// Definitions:
template<class T>
T TestArrayViewInput(const ArrayView<const Ptr<const T> >& a_in)
{
  typedef Teuchos::ScalarTraits<T> ST;
  T a = ST::zero();
  for (Teuchos_Ordinal i=0; i<a_in.size(); ++i) {
    a += *a_in[i];
  }
  return a;
}

template<class T>
void TestArrayViewOutput(const ArrayView<const Ptr<T> >& a_out)
{
  typedef Teuchos::ScalarTraits<T> ST;
  for (Teuchos_Ordinal i=0 ; i<a_out.size() ; ++i) {
    *a_out[i] = as<T>(i);
  }
}

template<class T>
Array<RCP<T> > generateArray(const Teuchos_Ordinal n)
{
  Array<RCP<T> > a(n);
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    RCP<T> data = rcp(new T(as<T>(i)));
    a[i] = data;
  }
  return a;
}


} // namespace ArrayConversionsUnitTestHelpers 

