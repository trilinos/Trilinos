#ifndef PHX_FIELDTAG_H
#define PHX_FIELDTAG_H

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Teuchos_RCP.hpp"

namespace PHX {

  class FieldTag {

  public:

    FieldTag(const std::string& name,
	     const Teuchos::RCP<PHX::DataLayout>& dl);
    
    ~FieldTag();

    FieldTag& operator=(const FieldTag& a);
    
    bool operator==(const FieldTag& a) const;
    
    bool operator<(const FieldTag& a) const;
    
    const std::string& name() const;

    const Teuchos::RCP<PHX::DataLayout> dataLayout() const;

    void print(std::ostream& os, int indent = 0) const;

  protected:

    std::string m_name;
    
    Teuchos::RCP<PHX::DataLayout> m_data_layout;

  };

  std::ostream& operator<<(std::ostream& os, const PHX::FieldTag& v);
  
} 

#endif 
