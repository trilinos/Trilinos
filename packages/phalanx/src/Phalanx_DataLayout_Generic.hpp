// @HEADER
// @HEADER

#ifndef PHX_DATA_LAYOUT_GENERIC
#define PHX_DATA_LAYOUT_GENERIC

#include "Phalanx_DataLayout.hpp"
#include <string>
#include <iostream>

namespace PHX{


  /*! \brief A concrete implementation of the DataLayout class that should cover most user requirements.

      This concrete class should be used for DataLayouts unless the
      user must pass specific external information via Data Layouts.

  */
  class Generic : public DataLayout {

  public:

    Generic(const std::string& unique_identifier, std::size_t i);

    virtual ~Generic();

    virtual bool operator==(const DataLayout& right) const;

    virtual const std::string& name() const;

    virtual std::size_t size() const;

    virtual const std::string identifier() const;

    virtual void print(std::ostream& os, int indent = 0) const;

  private:

    const std::string m_name;

    const std::size_t m_size;

  };

  std::ostream& operator<<(std::ostream& os, const PHX::Generic& t);

}

#endif
