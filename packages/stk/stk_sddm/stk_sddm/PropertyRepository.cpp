#include <iostream>

#include <stk_sddm/PropertyRepository.hpp>
#include <stk_sddm/Taxonomy.hpp>

namespace stk {
namespace sddm {

namespace {

void
invalid_taxon(
  const Property &      property,
  const std::string &   name)
{
  std::cerr << "Invalid taxon '" << name << "' does exist as a child of property '" << property.getName() << "' valid taxons are:" << std::endl;
}


void
validation_error(
  const Property &      property,
  const std::string &   message)
{
  std::cerr << "Validation error for property " << property.getName() << ": " << message << std::endl;
}

} // namespace <empty>


PropertyRepository::PropertyRepository(
  const std::string &   name,
  AnyValue *            value)
  : Property(name, value, *this),
    m_invalidTaxonFunc(invalid_taxon),
    m_validationErrorFunc(validation_error),
    m_taxonomies()
{}


PropertyRepository::~PropertyRepository()
{}


const Taxon *
PropertyRepository::match(
  const std::string &   message) const
{
  if (!m_taxonomies.empty()) {
    for (TaxonomyVector::const_iterator it = m_taxonomies.begin(); it != m_taxonomies.end(); ++it) {
      const Taxon *taxon = (*it)->match(message);
      if (taxon)
        return taxon;
    }
  }

  (*m_invalidTaxonFunc)(*this, message);

  return 0;
}

} // namespace sddm
} // namespace stk
