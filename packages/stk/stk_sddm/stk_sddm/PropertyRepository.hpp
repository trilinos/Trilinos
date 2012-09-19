#ifndef STK_SDDM_PROPERTY_REPOSITORY_HPP
#define STK_SDDM_PROPERTY_REPOSITORY_HPP

#include <stk_sddm/Property.hpp>

namespace stk {
namespace sddm {

class Taxonomy;

typedef std::vector<const Taxonomy *> TaxonomyVector;

/** 
 * @brief Class <b>PropertyRepository</b> 
 *
 */
class PropertyRepository : public Property
{
public:
  /** 
   * @brief Typedef <b>InvalidTaxonFunc</b> specifies the signature for
   *        the invalid taxon callback function.
   * 
   */
  typedef void (*ReportFunc)(const Property &property, const std::string &message);
  
  /** 
   * @brief <b>PropertyRepository</b> constructs a new property
   *        repository.  The property repository serves as the root to
   *        a property tree and contains a vector of taxonomies used
   *        to validate the insertion of new properties into the tree.
   * 
   * @param name        Specifies the name of the root property of
   *                    the property repository.
   *
   * @param value       Specifies an initial value for the root
   *                    property. 
   */
  PropertyRepository(const std::string &name, AnyValue *value = 0);

  /** 
   * @brief <b>~PropertyRepository</b> 
   * 
   */
  ~PropertyRepository();

  void addTaxonomy(const Taxonomy &taxonomy) {
    m_taxonomies.push_back(&taxonomy);
  }

  /** 
   * @brief <b>match</b> 
   * 
   * @param property 
   * @param name 
   * 
   * @return 
   */
  virtual const Taxon *match(const std::string &name) const;

  /** 
   * @brief <b>setInvalidTaxonFunction</b> 
   * 
   * @param invalid_taxon_func 
   */
  void setInvalidTaxonFunction(ReportFunc invalid_taxon_func) {
    m_invalidTaxonFunc = invalid_taxon_func;
  }

  void reportInvalidTaxon(const Property &property, const std::string &name) const {
    return (*m_invalidTaxonFunc)(property, name);
  }
  
  /** 
   * @brief <b>setInvalidTaxonFunction</b> 
   * 
   * @param invalid_taxon_func 
   */
  void setValidationErrorFunction(ReportFunc validation_error_func) {
    m_validationErrorFunc = validation_error_func;
  }

  void reportValidationError(const Property &property, const std::string &message) const {
    return (*m_validationErrorFunc)(property, message);
  }
  
private:
  ReportFunc            m_invalidTaxonFunc;             ///< Function to call on invalid taxon
  ReportFunc            m_validationErrorFunc;          ///< Function to call on validation error value
  TaxonomyVector        m_taxonomies;                   ///< Active taxonomies
};

} // namespace sddm
} // namespace stk

#endif // STK_SDDM_PROPERTY_REPOSITORY_HPP
