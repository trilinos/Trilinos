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
getFieldData(const PHX::FieldTag& v)
{
  typename std::map<PHX::FieldTag, Teuchos::ArrayRCP<DataT> >::iterator it;

  it = data_.find(v);
  if (it == data_.end()) {
    DataT data;
    std::ostringstream msg;
    msg << "The field:\n\n" << v
	<< "\n\ndoes not exist in DataContainer of type: " 
	<< typeid(data).name() << std::endl;
    TEST_FOR_EXCEPTION(it == data_.end(), std::logic_error, msg.str());
  }

  return it->second;
}

// ************************************************************************
template <typename DataT, typename Traits>
void PHX::DataContainer<DataT, Traits>::
allocateField(const PHX::FieldTag& v, std::size_t max_num_cells,
	      typename Traits::Allocator& a)
{
  std::size_t num_elements = v.dataLayout()->size() * max_num_cells;
  data_[v] = a.template allocate<DataT>(num_elements);
}

// ************************************************************************
template <typename DataT, typename Traits>
const std::type_info& PHX::DataContainer<DataT, Traits>::
getAlgebraicTypeInfo() const
{
  return typeid(typename boost::mpl::at<typename Traits::DataToAlgebraicMap, 
		DataT>::type);
}

// ************************************************************************
template <typename DataT, typename Traits>
void PHX::DataContainer<DataT, Traits>::
print(std::ostream& os) const
{
  std::string type = PHX::getTypeString<DataT, Traits>();
  
  using namespace std;
  os << "********************************************" << endl;
  os << "PHX::DataContainer Output" << endl;
  os << "********************************************" << endl;
  os << "  Data Type = " << type << endl;
  os << "  My FieldTags:";

  if (data_.size() == 0)
    os << " None!" << endl;
  else {
    os << endl;
    typename std::map< PHX::FieldTag, Teuchos::ArrayRCP<DataT> >::const_iterator it =
      data_.begin();
    for (; it != data_.end(); ++it)
      os << "    " << it->first << endl;
  }

  os << "********************************************" << endl;
}

// ************************************************************************
// ************************************************************************
#endif 
