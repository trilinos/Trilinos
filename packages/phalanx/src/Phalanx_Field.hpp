#ifndef PHX_FIELD_H
#define PHX_FIELD_H

#include <iostream>
#include <string>
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_FieldTag.hpp"

namespace PHX {

  template<typename DataT>
  class Field {
    
  public:

    typedef DataT value_type;
    
    Field(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& t);
    
    Field(const PHX::FieldTag& v);
    
    Field();
    
    ~Field();
    
    const PHX::FieldTag& fieldTag() const;

    DataT& operator[](int index);

    typename Teuchos::ArrayRCP<DataT>::Ordinal size() const;

    void setFieldTag(const PHX::FieldTag& v);
    
    void setFieldData(const Teuchos::ArrayRCP<DataT>& d);
    
    void print(std::ostream& os) const;

  private:
    
    PHX::FieldTag tag;
    
    Teuchos::ArrayRCP<DataT> field_data;

#ifdef PHX_DEBUG
    bool tag_set;
    bool data_set;
    static const std::string field_tag_error_msg;
    static const std::string field_data_error_msg;
#endif

  };
  
  template<typename DataT>
  std::ostream& operator<<(std::ostream& os, const PHX::Field<DataT>& h);
  
} 

#include "Phalanx_Field_Def.hpp"

#endif 
