// @HEADER
// @HEADER

#ifndef PHX_FIELDTAG_TAG_HPP
#define PHX_FIELDTAG_TAG_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_DataLayout.hpp"

namespace PHX {

  /*! \brief Typed Field Tag

      This class is a concrete implementation of the FieldTag base
      class that is templated on the data type to determine type
      information.

  */
  template<typename DataT>
  class Tag : public PHX::FieldTag {

  public:

    typedef DataT value_type;

    Tag(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& dl);
    
    ~Tag();

    Teuchos::RCP<FieldTag> clone() const;

    void operator=(const PHX::Tag<DataT>& t);
    
    bool operator==(const FieldTag& t) const;
    
    const std::string& name() const;

    const PHX::DataLayout& dataLayout() const;

    const std::type_info& dataTypeInfo() const;

    const std::string identifier() const;

    void print(std::ostream& os) const;

  protected:

    std::string m_name;
    
    Teuchos::RCP<PHX::DataLayout> m_data_layout;

  };

  template<typename DataT>
  std::ostream& operator<<(std::ostream& os, const PHX::Tag<DataT>& t);
  
} 

#include "Phalanx_FieldTag_Tag_Def.hpp"

#endif 
