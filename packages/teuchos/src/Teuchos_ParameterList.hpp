// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#include <complex>
#include <iostream>
#include <map>
#include <string>
// #include "Teuchos_ParameterArbitrary.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Teuchos {

class ParameterList
{
public:
  ParameterList();

  // templated constructor function
  template<typename CPType>
  ParameterList(std::string, CPType);

  ParameterList(const ParameterList &);
  ParameterList & operator= (const ParameterList &);

  template<typename SPType>
  void SetParameter(std::string, SPType);

  template<typename GPType>
  GPType GetParameter(std::string, GPType);

  void Print(int);

private:

  void PrintTabs(int);

  std::map<std::string, Entry<char> >                      CharMap;
  std::map<std::string, Entry<complex<double> > > ComplexDoubleMap;
  std::map<std::string, Entry<complex<float> > >   ComplexFloatMap;
  std::map<std::string, Entry<double> >                  DoubleMap;
  std::map<std::string, Entry<float> >                    FloatMap;
  std::map<std::string, Entry<int> >                        IntMap;
  std::map<std::string, Entry<std::string> >             StringMap;
  // std::map<std::string, Entry<Arbitrary*> >         ArbitraryMap;
  // std::map<std::string, Entry<ParameterList*> > ParameterListMap;
  
}; // class ParameterList

ParameterList::ParameterList()
{
  // nothing
}

template<typename CPType>
ParameterList::ParameterList(std::string name, CPType newData)
{
  cout << "Unsupported Parameter Type! Unable to add Parameter " << name << std::endl;
}

// Need a specialized constructor template function for each supported type
template<>
ParameterList::ParameterList(std::string name, char newData)
{
  Entry<char> newEntry(newData);
  CharMap[name] = newEntry;
}

template<>
ParameterList::ParameterList(std::string name, complex<double> newData)
{
  Entry<complex<double> > newEntry(newData);
  ComplexDoubleMap[name] = newEntry;
}

template<>
ParameterList::ParameterList(std::string name, complex<float> newData)
{
  Entry<complex<float> > newEntry(newData);
  ComplexFloatMap[name] = newEntry;
}

template<>
ParameterList::ParameterList(std::string name, double newData)
{
  Entry<double> newEntry(newData);
  DoubleMap[name] = newEntry;
}

template<>
ParameterList::ParameterList(std::string name, float newData)
{
  Entry<float> newEntry(newData);
  FloatMap[name] = newEntry;
}

template<>
ParameterList::ParameterList(std::string name, int newData)
{
  Entry<int> newEntry(newData);
  IntMap[name] = newEntry;
}

template<>
ParameterList::ParameterList(std::string name, std::string newData)
{
  Entry<std::string> newEntry(newData);
  StringMap[name] = newEntry;
}

// template<>
// ParameterList::ParameterList(std::string name, Arbitrary* newData)
// {
//   Entry<Arbitrary*> newEntry(newData);
//   ArbitraryMap[name] = newEntry;
// }

// template<>
// ParameterList::ParameterList(std::string name, ParameterList* newData)
// {
//   Entry<ParameterList*> newEntry(newData);
//   ParameterListMap[name] = newEntry;
// }

ParameterList::ParameterList(const ParameterList &Source)
{
  CharMap = Source.CharMap;
  ComplexDoubleMap = Source.ComplexDoubleMap;
  ComplexFloatMap = Source.ComplexFloatMap;
  DoubleMap = Source.DoubleMap;
  FloatMap = Source.FloatMap;
  IntMap = Source.IntMap;
  StringMap = Source.StringMap;
//   ArbitraryMap = Source.ArbitraryMap;
//   ParameterListMap = Source.ParameterListMap;

  std::cout << "copy constructor!!" << std::endl;
}

ParameterList & ParameterList::operator= (const ParameterList &Source)
{
  std::cout << "operator=!!!" << std::endl;
  if(this != &Source)
    {
      CharMap = Source.CharMap;
      ComplexDoubleMap = Source.ComplexDoubleMap;
      ComplexFloatMap = Source.ComplexFloatMap;
      DoubleMap = Source.DoubleMap;
      FloatMap = Source.FloatMap;
      IntMap = Source.IntMap;
      StringMap = Source.StringMap;
//       ArbitraryMap = Source.ArbitraryMap;
//       ParameterListMap = Source.ParameterListMap;
    }
  return *this;
}

template<typename SPType>
void ParameterList::SetParameter(std::string name, SPType newData)
{
  std::cout << "Unsupported Parameter Type! Unable to add Parameter " << name << std::endl;
}

// Need a specialized SetParameter function for each supported type
template<>
void ParameterList::SetParameter<char>(std::string name, char newData)
{
  Entry<char> newEntry(newData);
  CharMap[name] = newEntry;
}

template<>
void ParameterList::SetParameter<complex<double> >(std::string name, complex<double> newData)
{
  Entry<complex<double> > newEntry(newData);
  ComplexDoubleMap[name] = newEntry;
}

template<>
void ParameterList::SetParameter<complex<float> >(std::string name, complex<float> newData)
{
  Entry<complex<float> > newEntry(newData);
  ComplexFloatMap[name] = newEntry;
}

template<>
void ParameterList::SetParameter<double>(std::string name, double newData)
{
  Entry<double> newEntry(newData);
  DoubleMap[name] = newEntry;
}

template<>
void ParameterList::SetParameter<float>(std::string name, float newData)
{
  Entry<float> newEntry(newData);
  FloatMap[name] = newEntry;
}

template<>
void ParameterList::SetParameter<int>(std::string name, int newData)
{
  Entry<int> newEntry(newData);
  IntMap[name] = newEntry;
}

template<>
void ParameterList::SetParameter<std::string>(std::string name, std::string newData)
{
  Entry<std::string> newEntry(newData);
  StringMap[name] = newEntry;
}

// template<>
// void ParameterList::SetParameter<Arbitrary*>(std::string name, Arbitrary* newData)
// {
//   Entry<Arbitrary*> newEntry(newData);
//   ArbitraryMap[name] = newEntry;
// }

// template<>
// void ParameterList::SetParameter<ParameterList*>(std::string name, ParameterList* newData)
// {
//   Entry<ParameterList*> newEntry(newData);
//   ParameterListMap[name] = newEntry;
// }

template<typename GPType>
GPType ParameterList::GetParameter(std::string name, GPType nominal)
{
  std::cout << "Unsupported Parameter Type! Unable to retrieve Parameter " << name << std::endl;
  return nominal;
}

// Need a specialized GetParameter function for each supported type
template<>
char ParameterList::GetParameter(std::string name, char nominal)
{
  char result = nominal;
  std::map<std::string, Entry<char> >::iterator findIter = CharMap.find(name);
  if(findIter != CharMap.end())
    {
      result = (*findIter).second.GetData();
    }
  else
    {
      Entry<char> newEntry(nominal, 1);
      CharMap[name] = newEntry;
    }
  return result;
}

template<>
complex<double> ParameterList::GetParameter(std::string name, complex<double> nominal)
{
  complex<double> result = nominal;
  std::map<std::string, Entry<complex<double> > >::iterator findIter = ComplexDoubleMap.find(name);
  if(findIter != ComplexDoubleMap.end())
    {
      result = (*findIter).second.GetData();
    }
  else
    {
      Entry<complex<double> > newEntry(nominal, 1);
      ComplexDoubleMap[name] = newEntry;
    }
  return result;
}

template<>
complex<float> ParameterList::GetParameter(std::string name, complex<float> nominal)
{
  complex<float> result = nominal;
  std::map<std::string, Entry<complex<float> > >::iterator findIter = ComplexFloatMap.find(name);
  if(findIter != ComplexFloatMap.end())
    {
      result = (*findIter).second.GetData();
    }
  else
    {
      Entry<complex<float> > newEntry(nominal, 1);
      ComplexFloatMap[name] = newEntry;
    }
  return result;
}

template<>
double ParameterList::GetParameter(std::string name, double nominal)
{
  double result = nominal;
  std::map<std::string, Entry<double> >::iterator findIter = DoubleMap.find(name);
  if(findIter != DoubleMap.end())
    {
      result = (*findIter).second.GetData();
    }
  else
    {
      Entry<double> newEntry(nominal, 1);
      DoubleMap[name] = newEntry;
    }
  return result;
}

template<>
float ParameterList::GetParameter(std::string name, float nominal)
{
  float result = nominal;
  std::map<std::string, Entry<float> >::iterator findIter = FloatMap.find(name);
  if(findIter != FloatMap.end())
    {
      result = (*findIter).second.GetData();
    }
  else
    {
      Entry<float> newEntry(nominal, 1);
      FloatMap[name] = newEntry;
    }
  return result;
}

template<>
int ParameterList::GetParameter(std::string name, int nominal)
{
  cout << "GetParameter!!" << endl;
  int result = nominal;
  std::map<std::string, Entry<int> >::iterator findIter = IntMap.find(name);
  if(findIter != IntMap.end())
    {
      result = (*findIter).second.GetData();
    }
  else
    {
      Entry<int> newEntry(nominal, 1);
      IntMap[name] = newEntry;
    }
  return result;
}

template<>
std::string ParameterList::GetParameter(std::string name, std::string nominal)
{
  std::string result = nominal;
  std::map<std::string, Entry<std::string> >::iterator findIter = StringMap.find(name);
  if(findIter != StringMap.end())
    {
      result = (*findIter).second.GetData();
    }
  else
    {
      Entry<std::string> newEntry(nominal, 1);
      StringMap[name] = newEntry;
    }
  return result;
}

// template<>
// Arbitrary* ParameterList::GetParameter(std::string name, Arbitrary* nominal)
// {
//   Arbitrary* result = nominal;
//   std::map<std::string, Entry<Arbitrary*> >::iterator findIter = ArbitraryMap.find(name);
//   if(findIter != ArbitraryMap.end())
//     {
//       result = (*findIter).second.GetData();
//     }
//   else
//     {
//       Entry<Arbitrary*> newEntry(nominal, 1);
//       ArbitraryMap[name] = newEntry;
//     }
//   return result;
// }

// template<>
// ParameterList* ParameterList::GetParameter(std::string name, ParameterList* nominal)
// {
//   ParameterList* result = nominal;
//   std::map<std::string, Entry<ParameterList*> >::iterator findIter = ParameterListMap.find(name);
//   if(findIter != ParameterListMap.end())
//     {
//       result = (*findIter).second.GetData();
//     }
//   // How should GetParameter() behave if we are trying to get a list parameter?
//   return result;
// }

void ParameterList::Print(int indent)
{
  if(!CharMap.empty())
    {
      std::map<std::string, Entry<char> >::iterator charIter = CharMap.begin();
      while(charIter != CharMap.end())
	{
	  PrintTabs(indent);
	  std::cout << "[(char) " << (*charIter).first << " = '" << (*charIter).second.GetData() << "']" << std::endl;
	  charIter++;
	}
    }
  if(!ComplexDoubleMap.empty())
    {
      std::map<std::string, Entry<complex<double> > >::iterator complexdoubleIter = ComplexDoubleMap.begin();
      while(complexdoubleIter != ComplexDoubleMap.end())
	{
	  PrintTabs(indent);
	  std::cout << "[(complex<double>) " << (*complexdoubleIter).first << " = " << (*complexdoubleIter).second.GetData() << "]" << std::endl;
	  complexdoubleIter++;
	}
    }
  if(!ComplexFloatMap.empty())
    {
      std::map<std::string, Entry<complex<float> > >::iterator complexfloatIter = ComplexFloatMap.begin();
      while(complexfloatIter != ComplexFloatMap.end())
	{
	  PrintTabs(indent);
	  std::cout << "[(complex<float>) " << (*complexfloatIter).first << " = " << (*complexfloatIter).second.GetData() << "]" << std::endl;
	  complexfloatIter++;
	}
    }
  if(!DoubleMap.empty())
    {
      std::map<std::string, Entry<double> >::iterator doubleIter = DoubleMap.begin();
      while(doubleIter != DoubleMap.end())
	{
	  PrintTabs(indent);
	  std::cout << "[(double) " << (*doubleIter).first << " = " << (*doubleIter).second.GetData() << "]" << std::endl;
	  doubleIter++;
	}
    }
  if(!FloatMap.empty())
    {
      std::map<std::string, Entry<float> >::iterator floatIter = FloatMap.begin();
      while(floatIter != FloatMap.end())
	{
	  PrintTabs(indent);
	  std::cout << "[(float) " << (*floatIter).first << " = " << (*floatIter).second.GetData() << "]" << std::endl;
	  floatIter++;
	}
    }
  if(!IntMap.empty())
    {
      std::map<std::string, Entry<int> >::iterator intIter = IntMap.begin();
      while(intIter != IntMap.end())
	{
	  PrintTabs(indent);
	  std::cout << "[(int) " << (*intIter).first << " = " << (*intIter).second.GetData() << "]" << std::endl;
	  intIter++;
	}
    }
  if(!StringMap.empty())
    {
      std::map<std::string, Entry<std::string> >::iterator stringIter = StringMap.begin();
      while(stringIter != StringMap.end())
	{
	  PrintTabs(indent);
	  std::cout << "[(string) " << (*stringIter).first << " = \"" << (*stringIter).second.GetData() << "\"]" << std::endl;
	  stringIter++;
	}
    }
//   if(!ArbitraryMap.empty())
//     {
//       std::map<std::string, Entry<Arbitrary> >::iterator ArbitraryIter = ArbitraryMap.begin();
//       while(ArbitraryIter != ArbitraryMap.end())
// 	{
// 	  PrintTabs(indent);
// 	  std::cout << "[(Arbitrary) " << (*ArbitraryIter).first << " = " << (*ArbitraryIter).second.GetData() << "]" << std::endl;
// 	  ArbitraryIter++;
// 	}
//     }
//   if(!ParameterListMap.empty())
//     {
//       std::map<std::string, Entry<ParameterList*> >::iterator ParameterListIter = ParameterListMap.begin();
//       while(ParameterListIter != ParameterListMap.end())
// 	{
// 	  PrintTabs(indent);
// 	  std::cout << "[(ParameterList) " << (*ParameterListIter).first << ":]" << std::endl;
// 	  (*(*ParameterListIter).second.GetData()).Print(indent + 1);
// 	  ParameterListIter++;
// 	}
//     }
}

void ParameterList::PrintTabs(int n)
{
  int i;
  for(i = 0; i < n; i++)
    {
      std::cout << "\t";
    }
}

} // namespace Teuchos
