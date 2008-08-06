// @HEADER
// @HEADER

#ifndef PHX_FIELDTAG_HPP
#define PHX_FIELDTAG_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"

namespace PHX {

  class DataLayout;

  class FieldTag {

  public:

    FieldTag() {}
    
    virtual ~FieldTag() {}

    virtual Teuchos::RCP<FieldTag> clone() const = 0;

    virtual bool operator==(const FieldTag& t) const = 0;
    
    virtual bool operator!=(const FieldTag& t) const
    { return !(*this == t); };
    
    virtual const std::string& name() const = 0;

    virtual const PHX::DataLayout& dataLayout() const = 0;

    virtual const std::type_info& dataTypeInfo() const = 0;
    
    //! Unique name identifier that can be used for strict weak ordering in stl std::map keys.
    virtual const std::string identifier() const = 0;

    virtual void print(std::ostream& os) const = 0;

  };

  std::ostream& operator<<(std::ostream& os, const PHX::FieldTag& t)
  { 
    t.print(os); 
    return os; 
  }
  
} 

#endif 
