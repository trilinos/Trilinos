#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include <stk_sddm/Property.hpp>
#include <stk_sddm/Taxonomy.hpp>
#include <stk_sddm/Type.hpp>

namespace stk {
namespace sddm {

namespace {

struct taxon_name_finder
{
  taxon_name_finder(const std::string &name)
    : m_name(name)
  {}
  
  bool operator()(const Taxon *taxon) {
    return taxon->getName() == m_name;
  }
  
private:
  std::string m_name;
};


struct taxon_id_finder
{
  taxon_id_finder(const TaxonId &id)
    : m_id(id)
  {}
  
  bool operator()(const Taxon *taxon) {
    return taxon->getId() == m_id;
  }
  
private:
  TaxonId         m_id;
};

} // namespace <unnamed>


Taxon::Taxon(
  Taxonomy &            taxonomy,
  const TaxonId &       id,
  const std::string &   name)
  : m_id(id),
    m_name(name),
    m_type(0),
    m_taxonomy(taxonomy),
    m_children()
{}


Taxon::~Taxon()
{}


void
Taxon::addChild(
  Taxon *                 taxon)
{
  // TaxonVector::const_iterator it = std::find_if(m_children.begin(), m_children.end(), taxon_id_finder(taxon->getId()));
  // if (it != m_children.end())
  //   std::cout << "duplicate id " << taxon->getId() << " " << taxon->getName() << " already in " << getId() << " " << getName() << std::endl;
  
  TaxonVector::const_iterator it = std::find_if(m_children.begin(), m_children.end(), taxon_name_finder(taxon->getName()));
  if (it != m_children.end())
    ; // std::cout << "duplicate taxon " << taxon->getId() << " " << taxon->getName() << " already in " << getId() << " " << getName() << std::endl;
  else
    m_children.push_back(taxon);
}


Taxon &
Taxon::setType(
  const AnyType *       type)
{
  if (!m_taxonomy.findType(type)) {
    std::ostringstream oss;
    oss << "Cannot set the type for taxon " << getName() << " because the type " << type->name() << " has not been registered with the taxonomy " << m_taxonomy.getName();
    
    throw std::runtime_error(oss.str());
  }
  
  m_type = type;
  
  return *this;
}


const AnyType *
Taxonomy::findType(
  const std::string &   type_name) const
{
  for (std::vector<const AnyType *>::const_iterator it = m_dataTypes.begin(); it != m_dataTypes.end(); ++it)
    if (type_name == (*it)->name())
      return (*it);
  return 0;
}
  

const AnyType *
Taxonomy::findType(
  const AnyType *       type) const
{
  std::vector<const AnyType *>::const_iterator it = std::find(m_dataTypes.begin(), m_dataTypes.end(), type);
  if (it != m_dataTypes.end())
    return (*it);
  else
    return 0;
}
  

const Taxon *
Taxon::match(
  const std::string &   name) const
{
  for (TaxonVector::const_iterator it = m_children.begin(); it != m_children.end(); ++it) {
    if ((*it)->getName() == name)
      return (*it);
  }

  return 0;
}


Taxon *
Taxon::match(
  const std::string &   name)
{
  for (TaxonVector::const_iterator it = m_children.begin(); it != m_children.end(); ++it) {
    if ((*it)->getName() == name)
      return (*it);
  }

  return 0;
}


bool
Taxon::validate(
  const Property &      property) const
{
  bool result = true;
  
  for (AnnotationMap::const_iterator it = m_annotationMap.begin(); it != m_annotationMap.end(); ++it)
    result &= (*it)->validate(property);

  return result;
}


Taxonomy::~Taxonomy()
{
  destroy_and_clear();
}

void Taxonomy::destroy_and_clear()
{
  while (!m_children.empty()) {
    delete m_children.back();
    m_children.pop_back();
  }

  //vector swap trick to actually delete memory for m_children
  TaxonVector tmp;
  tmp.swap(m_children);

  delete m_root;
  m_root = NULL;
}

Taxon *
Taxonomy::createTaxon(
  const TaxonId &       id,
  const std::string &   name)
{
  TaxonId taxon_id = id;
  Taxon *new_taxon = new Taxon(*this, taxon_id, name);
  m_children.push_back(new_taxon);
  return new_taxon;
}


Taxon *
Taxonomy::makeGroup(
  const TaxonId &         id,
  const std::string &   name)
{
  Taxon *group = match(name);
  if (!group) {
    group = createTaxon(id, name);
    m_root->addChild(group);
  }
  
  return group;
}


Taxon *
Taxonomy::match(
  const std::string &   name)
{
  return m_root->match(name);
}


const Taxon *
Taxonomy::match(
  const std::string &   name) const
{
  return m_root->match(name);
}

} // namespace sddm
} // namespace stk
