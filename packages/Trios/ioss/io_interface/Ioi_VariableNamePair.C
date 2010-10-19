/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioi_VariableNamePair.h>
#include <Ioss_EntityBlock.h>
#include <assert.h>

namespace Ioi {
  VariableNamePair::VariableNamePair()
    : isValid(false), wasWritten(false), isExternal(false), subsetInclude(false), subsetExclude(false),
externalNameSpecified(false)
  {}

  VariableNamePair::VariableNamePair( const std::string& i, const std::string& e )
    : intName(i.c_str()), extName(e.c_str()), isValid(false), wasWritten(false), isExternal(false),
subsetInclude(false), subsetExclude(false), externalNameSpecified(false)

  {}

  VariableNamePair::VariableNamePair( const char *i, const char *e )
    : intName(i), extName(e), isValid(false), wasWritten(false), isExternal(false),
subsetInclude(false), subsetExclude(false), externalNameSpecified(false)
  {}

  VariableNamePair::VariableNamePair( const VariableNamePair& n )
    : intName(n.intName), extName(n.extName), entities_(n.entities_),
isValid(n.isValid), wasWritten(n.wasWritten), isExternal(n.isExternal),
subsetInclude(n.subsetInclude), subsetExclude(n.subsetExclude),
externalNameSpecified(n.externalNameSpecified)
  {}

  VariableNamePair&
  VariableNamePair::operator=( const VariableNamePair& n )
  {
    intName       = n.intName;
    extName       = n.extName;
    isValid       = n.isValid;
    wasWritten    = n.wasWritten;
    isExternal    = n.isExternal;
    externalNameSpecified = n.externalNameSpecified;

    subsetInclude = n.subsetInclude;
    subsetExclude = n.subsetExclude;

    if (!n.entities_.empty())
entities_     = n.entities_;
    return *this;
  }

  bool VariableNamePair::operator==( const VariableNamePair& n ) const
  {
    return intName == n.intName && extName == n.extName && isValid == n.isValid
&& wasWritten == n.wasWritten 
      && isExternal == n.isExternal
&& externalNameSpecified == n.externalNameSpecified
      && subsetInclude == n.subsetInclude
&& subsetExclude == n.subsetExclude;

  }

  bool VariableNamePair::operator!=( const VariableNamePair& n ) const
  {
    return !(*this == n);
  }

  void VariableNamePair::internal_name(const std::string& new_name) const
  {
    intName = new_name;
  }

  const std::string& VariableNamePair::internal_name() const { return intName; }

  const std::string& VariableNamePair::external_name() const { return extName; }

  bool VariableNamePair::is_valid() const
  {
    return isValid;
  }

  const VariableNamePair &VariableNamePair::set_valid(bool validity) const
  {
    isValid = validity;
    return *this;
  }

  bool VariableNamePair::was_written() const
  {
    return wasWritten;
  }

  const VariableNamePair& VariableNamePair::set_written(bool was_written) const
  {
    wasWritten = was_written;
    return *this;
  }


  bool VariableNamePair::external() const
  {
    return isExternal;
  }

  const VariableNamePair& VariableNamePair::set_external(bool is_external) const
  {
    isExternal = is_external;
    return *this;
  }

  void VariableNamePair::set_subset_include()
  {
    subsetInclude = true;
    assert(!subsetExclude);
  }

  void VariableNamePair::set_subset_exclude()
  {
    subsetExclude = true;
    assert(!subsetInclude);
  }

  void VariableNamePair::add_entity_name(const std::string& entity_name)
  {
    assert(subsetInclude || subsetExclude);

    // Could do a check to see if entity is already in list, but there should
    // be no problems if the same entity is in the list multiple times...
    entities_.push_back(const_cast<std::string&>(entity_name));
  }

  bool VariableNamePair::apply_to_entity(const Ioss::GroupingEntity *entity) const
  {
    assert(!(subsetInclude && subsetExclude));

    // If neither subset option specified, then variable is valid for all entities.
    if (!subsetInclude && !subsetExclude)
return true;

    // NOTE: Need to handle aliases here.  For now brute force the check.
    bool in_list = false;
    std::vector<std::string>::const_iterator I = entities_.begin();
    while (I != entities_.end() && !in_list) {
if (entity->is_alias(*I)) {
  in_list = true;
}
++I;
    }

    if (subsetInclude) {
// List specifies the entities that are to be included...
return in_list;
    } else {
assert(subsetExclude);
// List specifies the entities that are to be excluded...
return !in_list;
    }
  }

  void VariableNamePair::remove_subsetting()
  {
    // If neither subset option specified, just return; nothing to do.
    if (!subsetInclude && !subsetExclude)
return;
    
    subsetInclude = false;
    subsetExclude = false;
    std::vector<std::string>().swap(entities_);
    assert(entities_.empty());
  }
}
