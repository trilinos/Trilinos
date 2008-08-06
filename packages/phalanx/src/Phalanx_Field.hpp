// @HEADER
// @HEADER

#ifndef PHX_FIELD_H
#define PHX_FIELD_H

#include <iostream>
#include <string>
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_FieldTag_Tag.hpp"

namespace PHX {

  template<typename DataT>
  class Field {
    
  public:

    typedef DataT value_type;
    
    Field(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& t);
    
    Field(const PHX::Tag<DataT>& v);
    
    Field();
    
    ~Field();
    
    const PHX::FieldTag& fieldTag() const;

    DataT& operator[](int index);

    typename Teuchos::ArrayRCP<DataT>::Ordinal size() const;

    void setFieldTag(const PHX::Tag<DataT>& t);
    
    void setFieldData(const Teuchos::ArrayRCP<DataT>& d);
    
    void print(std::ostream& os) const;

  private:
    
    PHX::Tag<DataT> m_tag;
    
    Teuchos::ArrayRCP<DataT> m_field_data;

#ifdef PHX_DEBUG
    bool m_tag_set;
    bool m_data_set;
    static const std::string m_field_tag_error_msg;
    static const std::string m_field_data_error_msg;
#endif

  };
  
  template<typename DataT>
  std::ostream& operator<<(std::ostream& os, const PHX::Field<DataT>& h);
  
} 

#include "Phalanx_Field_Def.hpp"

#endif 
