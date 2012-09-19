#include <fstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include <stk_sddm/Taxonomy.hpp>
#include <stk_sddm/Property.hpp>
#include <stk_sddm/Type.hpp>

namespace sierra {
namespace Prsr {

std::ostream &printCommandSpecForTaxonomy(std::ostream &os, const std::string &key_id);

} // namespace Prsr
} // namespace sierra


namespace stk {
namespace sddm {

namespace {

std::string
serialize(
  const Taxon &                   taxon)
{
  std::string t;
  
  if (taxon.getName().find(':') != std::string::npos) {
    const std::string &s = taxon.getName();
    t += "\"";
    for (std::string::size_type i = 0; i < s.size(); ++i)
      if (s[i] == '"')
        t += "\\\"";
      else
        t += s[i];
    t += "\"";
  }
  else
    t = taxon.getName();

  // if (taxon.getOccurs().min == 0) {
  //   if (taxon.getOccurs().max == 1)
  //     t += '?';
  //   else if (taxon.getOccurs().max > 1)
  //     t += '*';
  // }
  // else if (taxon.getOccurs().min == 1) {
  //   if (taxon.getOccurs().max == 1)
  //     ;
  //   else if (taxon.getOccurs().max > 1)
  //     t += '+';
  // }

  return t;
}


std::ostream &
dump_children(
  std::ostream &                os,
  const Taxon &                 taxon,
  std::vector<std::string> &    active_parents)
{
  os << std::setw(active_parents.size()*2) << "";
  os << serialize(taxon);
  if (taxon.getType())
    os  << ": " << *taxon.getType();
  os << std::endl;
  
  if (active_parents.size() < 4) {
    active_parents.push_back(taxon.getName());
    
    for (TaxonVector::const_iterator it = taxon.getChildren().begin(); it != taxon.getChildren().end(); ++it) 
      dump_children(os, *(*it), active_parents);  

    active_parents.pop_back();
  }
  
  return os;
}


std::ostream &
print_children(
  std::ostream &	        os,
  const Taxon &                 taxon,
  std::vector<const Taxon *> &  active_parents)
{
  static size_t depth = 0;
  static std::vector<std::vector<const Taxon *> > visited;

  if (depth == 0)
    visited.clear();

  active_parents.push_back(&taxon);

  for (std::vector<const Taxon *>::const_iterator it = active_parents.begin(); it != active_parents.end(); ++it)
    os << (*it)->getId() << " ";
  os << ", ";

  for (std::vector<const Taxon *>::const_iterator it = active_parents.begin(); it != active_parents.end(); ++it)
    os << (*it)->getName() << " ";
  os << ", ";
  
  // if (getType())
  //   os << "[" << getType()->name() << "]";

// #ifdef HACK
  sierra::Prsr::printCommandSpecForTaxonomy(os, taxon.getId());
  os << std::endl;
// #endif
  ++depth;

  // recurse depth == 1, no recursion; 2 or more, is recursive length
  size_t recurse_depth = active_parents.end() - std::find(active_parents.begin(), active_parents.end(), &taxon);

  bool recurse = false;
  for (stk::sddm::TaxonVector::const_iterator it = taxon.getChildren().begin(); it != taxon.getChildren().end() && !recurse; ++it) {
    if (recurse_depth > 1) 
      for (std::vector<std::vector<const Taxon *> >::const_iterator it2 = visited.begin(); it2 != visited.end() && !recurse; ++it2) 
        if ((*it2).size() > recurse_depth && std::equal((*it2).end() - recurse_depth + 1, (*it2).end(), active_parents.end() - recurse_depth))
          recurse = true;

    if (!recurse) {
      visited.push_back(active_parents);
      print_children(os, *(*it), active_parents);
    }
  }
  
  --depth;

  active_parents.pop_back();  
  
  return os;
}

} // namespace <unnamed>

std::ostream &
dump(
  std::ostream &        os,
  const Taxonomy &      taxonomy)
{
  // os << "Full Taxonomy: " << std::endl;
  // for (TaxonVector::const_iterator it = m_taxa.begin(); it != m_taxa.end(); ++it) {
  //   os << *(*it) << std::endl;
  //   for (TaxonVector::const_iterator it2 = (*it)->m_children.begin(); it2 != (*it)->m_children.end(); ++it2) {
  //     os << "  " << *(*it2) << std::endl;
  //   }
  // }
  
  // os << std::endl
  //    << "Taxonomy groups: " << std::endl;

  std::vector<std::string> active_parents;
  
  for (TaxonVector::const_iterator it = taxonomy.getRoot().getChildren().begin(); it != taxonomy.getRoot().getChildren().end(); ++it)
    dump_children(os, *(*it), active_parents);

  return os;
}


std::ostream &
xml(
  std::ostream &	os,
  const Taxon &         taxon)
{
  os << "  <Taxon id=\"" << taxon.getId() << "\" taxon=\"" << taxon.getName() << "\"";
  if (taxon.getType())  
    os << " type=\"" << taxon.getType()->name() << "\"";
  if (taxon.getChildren().empty() && taxon.getAnnotations().first == taxon.getAnnotations().second) 
    os << "/>" << std::endl;
  else {
    os << ">" << std::endl;
    for (stk::sddm::AnnotationMap::const_iterator it = taxon.getAnnotations().first; it != taxon.getAnnotations().second; ++it) {
      (*it)->xml(os);
    }
    for (stk::sddm::TaxonVector::const_iterator it = taxon.getChildren().begin(); it != taxon.getChildren().end(); ++it)
      os << "    <Child id=\"" << (*it)->getId() << "\"/>" << std::endl;
    os << "  </Taxon>" << std::endl;
  }
  
  return os;
}


std::ostream &
xml(
  std::ostream &        os,
  const Taxonomy &      taxonomy)
{
  os << "<Taxonomy name=\"" << taxonomy.getName() << "\">" << std::endl;

  xml(os, taxonomy.getRoot());
  
  for (TaxonVector::const_iterator it = taxonomy.getTaxa().begin(); it != taxonomy.getTaxa().end(); ++it)
    xml(os, *(*it));
  
  // for (TaxonVector::const_iterator it = m_taxa.begin(); it != m_taxa.end(); ++it) {
  //   (*it)->xml(os);
  //   os << "  <Group id=\"" << (*it)->getId() << "\" taxon=\"" << (*it)->getTaxon() << "\">" << std::endl; 
 
  //   for (TaxonVector::const_iterator it2 = (*it)->m_children.begin(); it2 != (*it)->m_children.end(); ++it2) {
  //     (*it2)->xml(os);
  //   }
  // }

  os << "</Taxonomy>" << std::endl;

  return os;
}


std::ostream &
printGrammar(
  std::ostream &        os,
  const Taxonomy &      taxonomy)
{
  std::vector<const Taxon *>    active_parents;
  
  for (TaxonVector::const_iterator it = taxonomy.getRoot().getChildren().begin(); it != taxonomy.getRoot().getChildren().end(); ++it) {
    os << "Group" << std::endl;
    print_children(os, *(*it), active_parents);
  }
  
  return os;
}


void
load(
  std::istream &                        is,
  Taxonomy &                            taxonomy) 
{
  std::string s;

  std::vector<Taxon *> taxon_levels;
  
//  getline(is, s);
  
  while (getline(is, s)) {
    std::string::size_type i = 0;
    while (i < s.size() && s[i] == ' ')
      ++i;

    if (i == s.size() || s[i] == '@')
      continue;
    
    size_t level = i/2;
    while (taxon_levels.size() > level)
      taxon_levels.pop_back();

    if (level > taxon_levels.size())
      throw std::runtime_error("level trouble");
    
    std::string taxon_string;

    bool quote = false;
    while (i < s.size()) {
      if (!quote && s[i] == ':')
        break;
      
      if (s[i] == '\"')
        quote = !quote;
      else if (s[i] == '\\')
        taxon_string += s[++i];
      else
        taxon_string += s[i];
      ++i;
    }

    // Taxon::Occurs occurs(1, 1);
    // const char last_char = taxon_string[taxon_string.size() - 1];
    // if (last_char == '+') {
    //   occurs.max = Taxon::Occurs::UNLIMITED;
    //   taxon_string = taxon_string.substr(0, taxon_string.size() - 1);
    // }
    // else if (last_char == '*') {
    //   occurs.min = 0;
    //   occurs.max = Taxon::Occurs::UNLIMITED;
    //   taxon_string = taxon_string.substr(0, taxon_string.size() - 1);
    // }
    // else if (last_char == '?') {
    //   occurs.min = 0;
    //   occurs.max = 1;
    //   taxon_string = taxon_string.substr(0, taxon_string.size() - 1);
    // }

    Taxon *taxon = 0;
    if (level == 0) {
      taxon = taxonomy.makeGroup("x", taxon_string);
      taxon_levels.push_back(taxon);
    }
    else {
      taxon = taxonomy.createTaxon("x", taxon_string);
      taxon_levels.push_back(taxon);
      taxon_levels[level - 1]->addChild(taxon_levels.back());
    }

    // taxon->setOccurs(occurs);
      
    if (i < s.size() && s[i] == ':') {
      ++i;

      while (i < s.size() && s[i] == ' ')
        ++i;

      std::string type_string;
      while (i < s.size() && s[i] != ' ')
        type_string += s[i++];

      taxon->setType(taxonomy.findType(type_string));
    }

    // if (i < size() && s[i] == '(')
  }
}


std::ostream &
operator<<(
  std::ostream &        os,
  const Taxonomy &      taxonomy)
{
  dump(os, taxonomy);

  return os;
}


std::ostream &
operator<<(
  std::ostream &        os,
  const Taxon &           taxon)
{
  os << serialize(taxon);
  if (taxon.getType())
    os  << ": " << *taxon.getType();
  
  return os;
}

} // namespace sddm
} // namespace stk
