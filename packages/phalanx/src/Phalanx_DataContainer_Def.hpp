// @HEADER
// @HEADER

#ifndef PHX_DATA_CONTAINER_DEF_HPP
#define PHX_DATA_CONTAINER_DEF_HPP

#include "Teuchos_TestForException.hpp"
#include <iostream>
#include <typeinfo>
#include <sstream>
#include <iterator>
#include "boost/mpl/at.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_DebugStrings.hpp"

// ************************************************************************
template <typename DataT, typename Traits>
Teuchos::ArrayRCP<DataT> PHX::DataContainer<DataT, Traits>::
getFieldData(const PHX::FieldTag& t)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  typename map< RCP<const FieldTag>, ArrayRCP<DataT>, FTComp >::iterator it;
  it = m_data.find(Teuchos::rcp(&t, false));

  if (it == m_data.end()) {
    std::string type = PHX::getTypeString<DataT, Traits>();
    std::ostringstream msg;
    msg << "The field:\n\n" << t
	<< "\n\ndoes not exist in DataContainer of type: " 
	<< type << std::endl;
    TEST_FOR_EXCEPTION(it == m_data.end(), std::logic_error, msg.str());
  }

  return it->second;
}

// ************************************************************************
template <typename DataT, typename Traits>
void PHX::DataContainer<DataT, Traits>::
allocateField(const Teuchos::RCP<PHX::FieldTag>& t, 
	      std::size_t max_num_cells,
	      typename Traits::Allocator& a)
{
  std::size_t num_elements = t->dataLayout().size() * max_num_cells;
  m_data[t] = a.template allocate<DataT>(num_elements);
}

// ************************************************************************
template <typename DataT, typename Traits>
const std::type_info& PHX::DataContainer<DataT, Traits>::
dataTypeInfo() const
{
  return typeid(DataT);
}

// ************************************************************************
template <typename DataT, typename Traits>
std::size_t PHX::DataContainer<DataT, Traits>::
getSizeOfDataType() const
{
  std::size_t size = sizeof(DataT);
  return size;
}

// ************************************************************************
template <typename DataT, typename Traits>
void PHX::DataContainer<DataT, Traits>::
print(std::ostream& os) const
{
  using namespace std;
  using namespace Teuchos;

  std::string type = PHX::getTypeString<DataT, Traits>();
  
  os << "********************************************" << endl;
  os << "PHX::DataContainer Output" << endl;
  os << "********************************************" << endl;
  os << "  Data Type = " << type << endl;
  os << "  My FieldTags:";

  if (m_data.size() == 0)
    os << " None!" << endl;
  else {
    os << endl;
    typename map< RCP<const PHX::FieldTag>, ArrayRCP<DataT> >::const_iterator it = m_data.begin();
    for (; it != m_data.end(); ++it)
      os << "    " << *(it->first) << endl;
  }

  os << "********************************************" << endl;
}

// ************************************************************************
// ************************************************************************
#endif 
