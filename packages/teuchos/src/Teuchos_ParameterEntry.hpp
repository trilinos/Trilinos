// Kris
// 07.08.03 -- Move into Teuchos package/namespace



namespace Teuchos {

class ParameterList; // forward declaration so an Entry can be of type ParameterList itself

template<typename EntryType>
class Entry
{
public:
  Entry();
  Entry(EntryType newData, bool isCreatedByGet = 0);
  Entry(const Entry &Source);
  Entry & operator= (const Entry &Source);
  EntryType GetData();
  bool isUsed();

private:
  void reset();

  EntryType data;
  bool isGotten;
  bool isSetByGet;
}; // class Entry

template<typename EntryType>
Entry<EntryType>::Entry()
{
  // nothing
}

template<typename EntryType>
Entry<EntryType>::Entry(EntryType newData, bool isCreatedByGet)
{
  data = newData;
  isSetByGet = isCreatedByGet;
}

template<typename EntryType>
Entry<EntryType>::Entry(const Entry &Source)
{
  data = Source.data;
  isGotten = Source.isGotten;
  isSetByGet = Source.isSetByGet;
}

template<typename EntryType>
Entry<EntryType> & Entry<EntryType>::operator= (const Entry &Source)
{
  if(this != &Source)
    {
      data = Source.data;
      isGotten = Source.isGotten;
      isSetByGet = Source.isSetByGet;
    }
  return *this;
}

template<typename EntryType>
EntryType Entry<EntryType>::GetData()
{
  return data;
}

template<typename EntryType>
bool Entry<EntryType>::isUsed()
{
  return isGotten;
}

template<typename EntryType>
void Entry<EntryType>::reset()
{
  // What should reset do?
}

} // namespace Teuchos
