/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef Ioi_VariableNamePair_h
#define Ioi_VariableNamePair_h

#include <vector>
#include <string>

namespace Ioss {
  class GroupingEntity;
}

namespace Ioi {
  /*! the Ioi::VariableNamePair class is used to associate an external
   *  (i.e. for a results file) name with an internal sierra name for
   *  a variable
   */
  class VariableNamePair {
  public:
    VariableNamePair();
    VariableNamePair( const std::string& int_name, const std::string& ext_name );
    VariableNamePair( const char *int_name, const char *ext_name );
    VariableNamePair( const VariableNamePair& );
    VariableNamePair& operator=( const VariableNamePair& );
    bool operator==( const VariableNamePair& ) const;
    bool operator!=( const VariableNamePair& ) const;
    void internal_name(const std::string& new_name) const;
    const std::string& internal_name() const;
    const std::string& external_name() const;
    bool is_valid() const;
    const VariableNamePair& set_valid(bool validity) const;
    bool was_written() const;
    const VariableNamePair& set_written(bool was_written) const;
    bool was_external_name_specified() const {return externalNameSpecified;}
    VariableNamePair& set_external_name_specified(bool was_specified)
    { externalNameSpecified = was_specified; return *this; }

    //
    //  Flag to record if the variable is being managed by the application rather than io
    // 
    bool external() const;
    const VariableNamePair& set_external(bool external_) const;


    bool subset_variable() const {return subsetInclude || subsetExclude;}
    bool subset_exclude() const  {return subsetExclude;}
    bool subset_include() const  {return subsetInclude;}
    void remove_subsetting();

    /*!
     * Determine whether this variable should be written for the
     * entity named 'entity_name'
     */
    bool apply_to_entity(const Ioss::GroupingEntity *entity) const;

    void set_subset_include();
    void set_subset_exclude();
    void add_entity_name(const std::string& entity_name);

  private:
    mutable std::string intName;
    std::string extName;

    /*!
     * An optional list of entities to subset this variable to the
     * specifed subset of entities specified (or not specified) in
     * the list. The inclusion/exclusion is specified by the
     * boolean variable includeEntities.
     */
    std::vector<std::string> entities_;

    mutable bool isValid;
    mutable bool wasWritten;
    mutable bool isExternal;
    /*!
     * True if the entities list specifies a subset of entities
     * applied to this variable; Only applies if the enttities list
     * is non-empty.
     */
    bool subsetInclude;
    /*!
     * True if it specifies the subset of entities to be excluded
     * from this variable.  Only applies if the enttities list is
     * non-empty.
     */
    bool subsetExclude;
    bool externalNameSpecified;  //! True if user specified the external name with "as external_name"
  };
}
#endif
