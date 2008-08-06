#ifndef PHX_DATA_CONTAINER_HPP
#define PHX_DATA_CONTAINER_HPP

#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_FieldTag_Comparison.hpp"
#include "Phalanx_DataContainer_Base.hpp"

namespace PHX {

  /*! \brief Container that holds all fields associated with a specific DataT.
    
      One DataContainer is instantiated for each data type in each
      evaluation type.

  */
  template <typename DataT, typename Traits>
  class DataContainer : public PHX::DataContainerBase<Traits> {
    
  public:

    DataContainer() {}
    
    ~DataContainer() {}
    
    Teuchos::ArrayRCP<DataT> getFieldData(const PHX::FieldTag& t);
    
    void allocateField(const Teuchos::RCP<PHX::FieldTag>& t,
		       std::size_t max_num_cells,
		       typename Traits::Allocator& a);

    const std::type_info& dataTypeInfo() const;

    std::size_t getSizeOfDataType() const;

    void print(std::ostream& os) const;
    
  private:
    
    std::map< Teuchos::RCP<const PHX::FieldTag>, 
	      Teuchos::ArrayRCP<DataT>, 
	      FTComp > 
    m_data;
    
  };
  
  template <typename DataT, typename Traits>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::DataContainer<DataT, Traits>& dc)
  {
    dc.print(os);
    return os;
  }

} 

#include "Phalanx_DataContainer_Def.hpp"

#endif 
