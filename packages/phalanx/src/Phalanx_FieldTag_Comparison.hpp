// @HEADER
// @HEADER

#ifndef PHX_FIELDTAG_COMPARISON_HPP
#define PHX_FIELDTAG_COMPARISON_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "Phalanx_FieldTag.hpp"

namespace PHX {

  //! Functor for RCP<FieldTag> comparisons in stl std::map.  This must follow strict weak ordering rules.
  struct FTComp
  {
    bool operator()(const Teuchos::RCP<const PHX::FieldTag> a, 
		    const Teuchos::RCP<const PHX::FieldTag> b) const
    {
      return strcmp(a->identifier().c_str(), b->identifier().c_str()) < 0;
    }
  };

  //! Predicate for searches of stl std::vector< RCP<FieldTag> > using find_if algorithm.
  struct FTPred {
    const Teuchos::RCP<FieldTag> m_a;
    FTPred(const Teuchos::RCP<FieldTag> a) : m_a(a) {}
    bool operator()(const Teuchos::RCP<FieldTag> b) 
    {
      if (*m_a == *b)
	return true;
      else 
	return false;
    }
  };

  //! Predicate for searches of stl std::vector< RCP<FieldTag> > using find_if algorithm.
  struct FTPredRef {
    const FieldTag& m_a;
    FTPredRef(const FieldTag& a) : m_a(a) {}
    bool operator()(const Teuchos::RCP<FieldTag> b) 
    {
      if (m_a == *b)
	return true;
      else 
	return false;
    }
  };

} 

#endif 
