#ifndef PHX_DATA_CONTAINER_HPP
#define PHX_DATA_CONTAINER_HPP

#include <iostream>
#include <map>
#include <vector>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_DataContainer_Base.hpp"

namespace PHX {
  
  /*! \brief Contain that holds fields associated with a specific data type.
    
  */
  template <typename DataT, typename Traits>
  class DataContainer : public PHX::DataContainerBase<Traits> {
    
  public:
    
    DataContainer() {}
    
    ~DataContainer() {}
    
    Teuchos::ArrayRCP<DataT> getFieldData(const PHX::FieldTag& v);
    
    void allocateField(const PHX::FieldTag& v,
		       std::size_t max_num_cells,
		       typename Traits::Allocator& a);

    const std::type_info& getAlgebraicTypeInfo() const;

    void print(std::ostream& os) const;
    
  private:
    
    std::map< PHX::FieldTag, Teuchos::ArrayRCP<DataT> > data_;
    
  };
  
  template <typename DataT, typename Traits>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::DataContainer<DataT, Traits>& v)
  {
    v.print(os);
    return os;
  }

} 

#include "Phalanx_DataContainer_Def.hpp"

#endif 
