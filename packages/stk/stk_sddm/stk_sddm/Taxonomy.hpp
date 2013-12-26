#ifndef SDDM_TAXONOMY_HPP
#define SDDM_TAXONOMY_HPP

#include <vector>
#include <map>
#include <set>
#include <string>
#include <iosfwd>

#include <stk_sddm/AnnotationMap.hpp>
#include <stk_sddm/AnnotationInterface.hpp>
#include <stk_sddm/DataTypeTraits.hpp>

namespace stk {
namespace sddm {

class AnyType;
class Taxonomy;
class Property;
class Taxon;

typedef std::vector<Taxon *> TaxonVector;
typedef std::string TaxonId;

class Taxon
{
  friend class Taxonomy;
  
private:
/** 
 * @brief <b>Taxon</b> constructs a new taxon node with specified id
 * and name.  Use Taxonomy::createTaxon() to create new nodes in a
 * taxonomy.
 * 
 * @param taxonomy      owning taxonomy
 * @param id            unique identification for this taxon
 * @param name          
 */
  Taxon(Taxonomy &taxonomy, const TaxonId &id, const std::string &name);

  /** 
   * @brief <b>~Taxon</b> 
   * 
   */
  ~Taxon();

  /** 
   * @brief <b>Taxon</b> 
   * 
   */
  Taxon(const Taxon &);

  /** 
   * @brief <b>operator=</b> 
   * 
   * 
   * @return 
   */
  Taxon &operator=(const Taxon &);
  
public:
  /** 
   * @brief <b>addChild</b> 
   * 
   * @param taxon 
   */
  void addChild(Taxon *taxon);

  /** 
   * @brief <b>match</b> returns the child node with the specified
   * name.  The search is case sensitive.
   * 
   * @param name        child's name to find
   * 
   * @return pointer to the child taxon with specified name or 0 if not found.
   */
  Taxon *match(const std::string &name);
  
  /** 
   * @brief <b>match</b> returns the child node with the specified
   * name.  The search is case sensitive.
   * 
   * @param name        child's name to find
   * 
   * @return pointer to the child taxon with specified name or 0 if not found.
   */
  const Taxon *match(const std::string &name) const;

  /** 
   * @brief <b>getId</b> returns the id of the taxon.
   * 
   * 
   * @return const reference to the id of the taxon.
   */
  const TaxonId &getId() const {
    return m_id;
  }
  
  /** 
   * @brief <b>getName</b> 
   * 
   * 
   * @return const reference to the name of the taxon.
   */
  const std::string &getName() const {
    return m_name;
  }
  
  /** 
   * @brief <b>getType</b> returns the AnyType description of the
   * taxon.  This description defines the type of data that may be
   * stored in the property associated with this taxon.
   * 
   * 
   * @return pointer to the AnyType object which describes the type of
   * property value that may be stored.
   */
  const AnyType *getType() const {
    return m_type;
  }
  
  /** 
   * @brief <b>setType</b> sets the AnyType description of the taxon.
   * This description defines the type of data that may be stored in
   * the property associated with this taxon.
   * 
   * @param any_type pointer to the AnyType object which describes the
   * type of property value that may be stored.  The taxonomy does NOT take responsibility to destroy the type.  Gener
   * 
   * @return reference to this taxon.
   */
  Taxon &setType(const AnyType *any_type);

  /** 
   * @brief <b>addAnnotation</b> adds an annotation to the taxon.  The
   * annotation provides additional information about the property
   * associated with the taxon.
   * 
   * @param annotation pointer to an annotation.  The taxonomy takes
   * ownership of the Annotation object and will destroy it when the
   * taxonomy is destroyed.
   * 
   * @return 
   */
  template <class T>
  Taxon &addAnnotation(AnnotationInterface<T> *annotation) {
    m_annotationMap.addAnnotation(annotation);
    return *this;
  }
  
  /** 
   * @brief <b>findAnnotation</b> returns a pair of iterators which
   * specify a range of annotations of the specified annotation
   * interface.
   */
  template <class T, class Out>
  void findAnnotations(Out out) const {
    m_annotationMap.findAnnotation<T>(out);
  }
  
  /** 
   * @brief <b>getAnnotations</b> 
   */
  std::pair<AnnotationMap::const_iterator, AnnotationMap::const_iterator> getAnnotations() const {
    return std::make_pair(m_annotationMap.begin(), m_annotationMap.end());
  }
  
  /** 
   * @brief <b>validate</b> performs annotation validation on
   * property.  All failing annotations will result in a
   * validation_error() call on the property's repository.
   * 
   * @param property    property to be validated
   * 
   * @return true if property passes all validation tests.
   */
  bool validate(const Property &property) const;
  
  /** 
   * @brief <b>getChildren</b> 
   */
  const TaxonVector &getChildren() const {
    return m_children;
  }

private:
  const TaxonId         m_id;                           ///< Taxon's unique identifier
  const std::string     m_name;                         ///< Name that may be placed in property tree
  const AnyType *       m_type;                         ///< Type of value associated with property
  Taxonomy &            m_taxonomy;                     ///< Owner taxonomy
  AnnotationMap         m_annotationMap;                ///< Annotation's applied to assocated property
  TaxonVector           m_children;                     ///< Children of this taxon
};


class Taxonomy 
{
public:
  /** 
   * @brief <b>Taxonomy</b> 
   * 
   * @param name 
   */
  Taxonomy(const char *name)
    : m_name(name),
      m_root(new Taxon(*this, name, name)),
      m_children(),
      m_dataTypes()
  {}

  /** 
   * @brief <b>~Taxonomy</b> 
   * 
   */
  ~Taxonomy();

private:
  Taxonomy(const Taxonomy &);
  Taxonomy &operator=(const Taxonomy &);

public:
  /** 
   * @brief <b>createTaxon</b> 
   * 
   * @param id 
   * @param name 
   * 
   * @return 
   */
  Taxon *createTaxon(const TaxonId &id, const std::string &name);

  /** 
   * @brief <b>makeGroup</b> 
   * 
   * @param id 
   * @param name 
   * 
   * @return 
   */
  Taxon *makeGroup(const TaxonId &id, const std::string &name);

  /** 
   * @brief <b>match</b> 
   * 
   * @param name 
   * 
   * @return 
   */
  Taxon *match(const std::string &name);

  /** 
   * @brief <b>match</b> 
   * 
   * @param name 
   * 
   * @return 
   */
  const Taxon *match(const std::string &name) const;

  /** 
   * @brief <b>getName</b> 
   * 
   * 
   * @return 
   */
  const std::string &getName() const {
    return m_name;
  }
  
  /** 
   * @brief <b>getRoot</b> 
   * 
   * 
   * @return 
   */
  const Taxon &getRoot() const {
    return *m_root;
  }
  
  /** 
   * @brief <b>getTaxa</b> 
   * 
   * 
   * @return 
   */
  const TaxonVector &getTaxa() const {
    return m_children;
  }
  
  /** 
   * @brief <b>findType</b> 
   * 
   * @param type_name 
   * 
   * @return 
   */
  const AnyType *findType(const std::string &type_name) const;

  /** 
   * @brief <b>findType</b> 
   * 
   * @param type
   * 
   * @return 
   */
  const AnyType *findType(const AnyType *type) const;

  /** 
   * @brief <b>findType</b> 
   * 
   * @return a const AnyType* 
   */
  template <class T>
  const AnyType *findType() const {
    return findType(DataTypeTraits<T>::name());
  }

  /** 
   * @brief <b>registerType</b> 
   * 
   * @param any_type 
   */
  void registerType(const AnyType *any_type) {
    m_dataTypes.push_back(any_type);
  }
  
  /** 
   * @brief <b>registerType</b> 
   */
  template <class It>
  void registerTypes(It first, It last) {
    std::copy(first, last, std::back_inserter(m_dataTypes));
  }
  
  void destroy_and_clear();

private:
  const std::string             m_name;                 ///< 
  Taxon *                       m_root;                 ///< 
  TaxonVector                   m_children;             ///< 
  std::vector<const AnyType *>  m_dataTypes;            ///< 
};

typedef std::vector<const Taxonomy *> TaxonomyVector;

} // namespace sddm
} // namespace stk

#endif // SDDM_TAXONOMY_HPP
