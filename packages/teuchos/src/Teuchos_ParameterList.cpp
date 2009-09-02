// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

//#define TEUCHOS_PARAMETER_LIST_SHOW_TRACE

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_StrUtils.hpp"

#ifdef TEUCHOS_PARAMETER_LIST_SHOW_TRACE
#include "Teuchos_VerboseObject.hpp"
#endif

/* NOTE: ASCI Red (TFLOP) does not support the i-> function for iterators 
 * in the STL.  Therefore when compiling for the TFLOP we must redefine the 
 * iterator from i-> to (*i). This slows things down on other platforms 
 * so we switch between the two when necessary.
 */

namespace {

std::string filterValueToString(const Teuchos::ParameterEntry& entry )
{
  return ( entry.isList() ? std::string("...") : toString(entry.getAny()) );
}

struct ListPlusValidList {
  Teuchos::ParameterList   *list;
  Teuchos::ParameterList   *validList;
  ListPlusValidList(
    Teuchos::ParameterList   *_list
    ,Teuchos::ParameterList  *_validList
    )
    :list(_list),validList(_validList)
    {}
};

} // namespace 

namespace Teuchos {

ParameterList::ParameterList()
  :name_("ANONYMOUS"), disableRecursiveValidation_(false)
{}

ParameterList::ParameterList(const std::string &name_in)
  :name_(name_in), disableRecursiveValidation_(false)
{}

ParameterList::ParameterList(const ParameterList& source)
{
  name_ = source.name_;
  params_ = source.params_;
  disableRecursiveValidation_ = source.disableRecursiveValidation_;
}

ParameterList& ParameterList::operator=(const ParameterList& source) 
{
  if (&source == this)
    return *this;
  name_ = source.name_;
  params_ = source.params_;
  disableRecursiveValidation_ = source.disableRecursiveValidation_;
  return *this;
}

ParameterList& ParameterList::setParameters(const ParameterList& source)
{
  for( ConstIterator i = source.begin(); i != source.end(); ++i ) {
    const std::string     &name_i  = this->name(i);
    const ParameterEntry  &entry_i = this->entry(i);
    if(entry_i.isList()) {
      this->sublist(name_i,false,entry_i.docString()).setParameters(
        getValue<ParameterList>(entry_i) );
    }
    else {
      this->setEntry(name_i,entry_i);
    }
  }
  this->updateSubListNames();
  return *this;
}

ParameterList& ParameterList::setParametersNotAlreadySet(
  const ParameterList& source
  ) 
{
  for( ConstIterator i = source.begin(); i != source.end(); ++i ) {
    const std::string     &name_i  = this->name(i);
    const ParameterEntry  &entry_i = this->entry(i);
    if(entry_i.isList()) {
      this->sublist(name_i,false,entry_i.docString()).setParametersNotAlreadySet(
        getValue<ParameterList>(entry_i) );
    }
    else {
      const ParameterEntry
        *thisEntryPtr = this->getEntryPtr(name_i);
      // If the entry does not already exist, then set it.  Otherwise, leave the
      // existing intery allow
      if(!thisEntryPtr)
        this->setEntry(name_i,entry_i);
    }
  }
  this->updateSubListNames();
  return *this;
}

ParameterList& ParameterList::disableRecursiveValidation()
{
  disableRecursiveValidation_ = true;
  return *this;
}

ParameterList::~ParameterList() 
{}

void ParameterList::unused(std::ostream& os) const
{
  for (ConstIterator i = params_.begin(); i != params_.end(); ++i) {
    if (!(entry(i).isUsed())) {
      os << "WARNING: Parameter \"" << name(i) << "\" " << entry(i)
         << " is unused" << std::endl;
    }
  }
}

std::string ParameterList::currentParametersString() const
{
  std::ostringstream oss;
  oss << "  {\n";
  ParameterList::ConstIterator itr;
  int i;
  for( itr = this->begin(), i = 0; itr != this->end(); ++itr, ++i ) {
    const std::string     &entryName = this->name(itr);
    const ParameterEntry  &theEntry = this->entry(itr);
    oss
      << "    \""<<entryName<<"\" : "<<theEntry.getAny().typeName()
      <<" = "<<filterValueToString(theEntry) << "\n";
  }
  oss << "  }\n";
  return oss.str();
}

bool ParameterList::isSublist(const std::string& name_in) const
{
  ConstIterator i = params_.find(name_in);
  if (i != params_.end())
    return (entry(i).isList());
  return false;
}

bool ParameterList::isParameter(const std::string& name_in) const
{
  return (params_.find(name_in) != params_.end());
}

bool ParameterList::remove(
  std::string const& name_in, bool throwIfNotExists
  )
{
  Iterator i = params_.find(name_in);
  TEST_FOR_EXCEPTION(
    throwIfNotExists && i == params_.end(), Exceptions::InvalidParameterName
    ,"Teuchos::ParameterList::remove(name,throwIfNotExists):"
    "\n\nError, the parameter \"" << name_in << "\" does not exist!"
    );
  if( i != params_.end() ) {
    params_.erase(i);
  }
  return false;
}

ParameterList& ParameterList::sublist(
  const std::string& name_in, bool mustAlreadyExist
  ,const std::string& docString
  )
{
  // Find name in list, if it exists.
  Iterator i = params_.find(name_in);

  // If it does exist and is a list, return the list value.
  // Otherwise, throw an error.
  if (i != params_.end()) {
#ifdef TEUCHOS_DEBUG
    const std::string actualName = this->name(i);
    TEST_FOR_EXCEPTION(
      name_in != actualName, std::logic_error,
      "Error, the sublist named \"" << name_in << "\" was said to be found\n"
      "but the actual parameter name is \"" << actualName << "\".\n"
      "This suggests some type of memory corruption in the list (try running a\n"
      "memory checking tool loke purify or valgrind)."
      );
#endif

    TEST_FOR_EXCEPTION_PURE_MSG(
      !entry(i).isList(), Exceptions::InvalidParameterType
      ,"Error, the parameter \"" << name_in << "\" is not a list, it is of type \""
      <<entry(i).getAny(false).typeName()<<"\"!" );
    return getValue<ParameterList>(entry(i));
  }

  // The list does not exist so create a new empty list and return its reference
  TEST_FOR_EXCEPTION_PURE_MSG(
    mustAlreadyExist, Exceptions::InvalidParameterName
    ,"The sublist "<<this->name()<<"->\""<<name_in<<"\" does not exist!"
    );
  const ParameterList newSubList(this->name()+std::string("->")+name_in);
  ParameterEntry &newParamEntry = params_.insert(
    Map::value_type(name_in,ParameterEntry(newSubList,false,true,docString))
    ).first->second;
  // Make sure we set the documentation std::string!
#ifdef TEUCHOS_DEBUG
  {
    ParameterEntry *newNewParamEntry = this->getEntryPtr(name_in);
    TEST_FOR_EXCEPTION(
      0 == newNewParamEntry, std::logic_error,
      "Error, the parameter was not set for sublist \"" << name_in << "\"!"
      );
    const std::string newDocString = newNewParamEntry->docString();
  TEST_FOR_EXCEPTION(
    newDocString != docString, std::logic_error,
    "Error, the set documentation std::string is not equal to the pass in std::string for\n"
    "the sublist \"" << name_in << "\"."
    );
  }
#endif
  return any_cast<ParameterList>(newParamEntry.getAny(false));
}

const ParameterList& ParameterList::sublist(const std::string& name_in) const
{
  // Find name in list, if it exists.
  ConstIterator i = params_.find(name_in);

  // If it does not exist, throw an error
  TEST_FOR_EXCEPTION_PURE_MSG(
    i == params_.end(), Exceptions::InvalidParameterName
    ,"Error, the sublist "<<this->name()<<"->\""<<name_in<<"\" does not exist!"
    );

  // If it does exist and is a list, return the list value.
  TEST_FOR_EXCEPTION_PURE_MSG(
    !entry(i).isList(), Exceptions::InvalidParameterType
    ,"Error, the parameter \""<<name_in<<"\" is not a list!  Instead it is of type"
    " \""<<entry(i).getAny(false).typeName()<<"\"!"
    );
  return getValue<ParameterList>(entry(i));
}
  
std::ostream& ParameterList::print(std::ostream& os, int indent, bool showTypes, bool showFlags) const
{
  return this->print(os,PrintOptions().indent(indent).showTypes(showTypes).showFlags(showFlags));
}
  
std::ostream& ParameterList::print(std::ostream& os, const PrintOptions &printOptions ) const
{
  const int   indent    = printOptions.indent();
  const bool  showTypes = printOptions.showTypes();
  const bool  showFlags = printOptions.showFlags();
  const bool  showDoc   = printOptions.showDoc();
  const std::string linePrefix(indent,' ');
  RCP<FancyOStream>
    out = getFancyOStream(rcp(&os,false));
  OSTab tab(out,indent);
  if (params_.begin() == params_.end()) {
    *out <<"[empty list]" << std::endl;
  }
  else { 
    // Print parameters first
    for (ConstIterator i = params_.begin(); i != params_.end(); ++i) 
    {
      const std::string &name_i = this->name(i);
      const ParameterEntry &entry_i = entry(i);
      RCP<const ParameterEntryValidator>
        validator = entry_i.validator();
      if(entry_i.isList())
        continue;
      *out << name_i;
      const std::string &docString = entry_i.docString();
      if(showTypes)
        *out << " : " << entry_i.getAny(false).typeName();
      *out << " = "; entry_i.leftshift(os,showFlags); *out << std::endl;
      if(showDoc) {
        if(validator.get()) {
          validator->printDoc(docString,OSTab(os).o());
        }
        else if( docString.length() ) {
          StrUtils::printLines(OSTab(out).o(),"# ",docString);
        }
      }
    }
    // Print sublists second
    for (ConstIterator i = params_.begin(); i != params_.end(); ++i) 
    {
      const ParameterEntry &entry_i = entry(i);
      if(!entry_i.isList())
        continue;
      const std::string &docString = entry_i.docString();
      const std::string &name_i = this->name(i);
      *out << name_i << " -> " << std::endl;
      if( docString.length() && showDoc ) {
        StrUtils::printLines(OSTab(out).o(),"# ",docString);
      }
      getValue<ParameterList>(entry_i).print(OSTab(out).o(), printOptions.copy().indent(0));
    }
  }
  return os;
}

ParameterList::ConstIterator ParameterList::begin() const
{
  return params_.begin();
}

ParameterList::ConstIterator ParameterList::end() const
{
  return params_.end();
}

#if defined(TFLOP)

const std::string& ParameterList::name(ConstIterator i) const
{
  return ((*i).first);
}

ParameterEntry& ParameterList::entry(Iterator i)
{
  return ((*i).second);
}

const ParameterEntry& ParameterList::entry(ConstIterator i) const
{
  return ((*i).second);
}

#else // defined(TFLOP)

const std::string& ParameterList::name(ConstIterator i) const
{
  return (i->first);
}

ParameterEntry& ParameterList::entry(Iterator i)
{
  return (i->second);
}

const ParameterEntry& ParameterList::entry(ConstIterator i) const
{
  return (i->second);
}

#endif // defined(TFLOP)


void ParameterList::validateParameters(
  ParameterList const& validParamList,
  int const depth,
  EValidateUsed const validateUsed,
  EValidateDefaults const validateDefaults
  ) const
{
  typedef std::deque<ListPlusValidList> sublist_list_t;
#ifdef TEUCHOS_PARAMETER_LIST_SHOW_TRACE
  RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();
  OSTab tab(out);
  *out << "\n*** Entering ParameterList::validateParameters(...) for "
    "this->name()=\""<<this->name()<<"\"...\n";
#endif
  //
  // First loop through and validate the parameters at this level.
  //
  // Here we generate a list of sublists that we will search next
  //
  sublist_list_t sublist_list;
  ConstIterator itr;
  for( itr = this->begin(); itr != this->end(); ++itr ) {
    const std::string    &entryName   = this->name(itr);
    const ParameterEntry &theEntry       = this->entry(itr);
#ifdef TEUCHOS_PARAMETER_LIST_SHOW_TRACE
    OSTab tab(out);
    *out << "\nentryName=\""<<entryName<<"\"\n";
#endif
    if(
      ( theEntry.isUsed() && validateUsed!=VALIDATE_USED_ENABLED )
      ||
      ( theEntry.isDefault() && validateDefaults!=VALIDATE_DEFAULTS_ENABLED )
      )
    {
      continue;
    }
    const ParameterEntry *validEntry = validParamList.getEntryPtr(entryName);
    TEST_FOR_EXCEPTION_PURE_MSG(
      !validEntry, Exceptions::InvalidParameterName
      ,"Error, the parameter {name=\""<<entryName<<"\","
      "type=\""<<theEntry.getAny(false).typeName()<<"\""
      ",value=\""<<filterValueToString(theEntry)<<"\"}"
      "\nin the parameter (sub)list \""<<this->name()<<"\""
      "\nwas not found in the list of valid parameters!"
      "\n\nThe valid parameters and types are:\n"
      <<validParamList.currentParametersString()
      );
    RCP<const ParameterEntryValidator> validator;
    if( (validator=validEntry->validator()).get() ) {
      validator->validate( theEntry, entryName, this->name() ); 
    }
    else {
      const bool validType =
        ( validEntry!=NULL
          ? theEntry.getAny(false).type() == validEntry->getAny(false).type()
          : false
          );
      TEST_FOR_EXCEPTION_PURE_MSG(
        !validType, Exceptions::InvalidParameterType
        ,"Error, the parameter {name=\""<<entryName<<"\","
        "type=\""<<theEntry.getAny(false).typeName()<<"\""
        ",value=\""<<filterValueToString(theEntry)<<"\"}"
        "\nin the parameter (sub)list \""<<this->name()<<"\""
        "\nexists in the list of valid parameters but has the wrong type."
        "\n\nThe correct type is \""
        << validEntry->getAny(false).typeName() << "\"."
        );
    }
    if( theEntry.isList() && depth > 0 ) {
      sublist_list.push_back(
        ListPlusValidList(
          &getValue<ParameterList>(theEntry),&getValue<ParameterList>(*validEntry)
          )
        );
    }
  }
  //
  // Now loop through the sublists and validate their parameters
  //
  for(
    sublist_list_t::const_iterator sl_itr = sublist_list.begin();
    sl_itr != sublist_list.end();
    ++sl_itr
    )
  {
    if (!sl_itr->validList->disableRecursiveValidation_) {
      sl_itr->list->validateParameters(
        *sl_itr->validList
        ,depth-1
        ,validateUsed
        ,validateDefaults
        );
    }
  }
#ifdef TEUCHOS_PARAMETER_LIST_SHOW_TRACE
  *out << "\n*** Existing ParameterList::validateParameters(...) for "
    "this->name()=\""<<this->name()<<"\"...\n";
#endif
}


void ParameterList::validateParametersAndSetDefaults(
  ParameterList const& validParamList,
  int const depth
  )
{
  typedef std::deque<ListPlusValidList> sublist_list_t;
#ifdef TEUCHOS_PARAMETER_LIST_SHOW_TRACE
  RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();
  OSTab tab(out);
  *out << "\n*** Entering ParameterList::validateParametersAndSetDefaults(...) "
    "for this->name()=\""<<this->name()<<"\"...\n";
#endif
  //
  // First loop through and validate the parameters at this level.
  //
  // Here we generate a list of sublists that we will search next
  //
  sublist_list_t sublist_list;
  {
    Iterator itr;
    for( itr = this->nonconstBegin(); itr != this->nonconstEnd(); ++itr ) {
      const std::string  &entryName = this->name(itr);
      ParameterEntry &theEntry = this->entry(itr);
#ifdef TEUCHOS_PARAMETER_LIST_SHOW_TRACE
      OSTab tab(out);
      *out << "\nentryName=\""<<entryName<<"\"\n";
#endif
      const ParameterEntry *validEntry = validParamList.getEntryPtr(entryName);
      TEST_FOR_EXCEPTION_PURE_MSG(
        !validEntry, Exceptions::InvalidParameterName
        ,"Error, the parameter {name=\""<<entryName<<"\","
        "type=\""<<theEntry.getAny(false).typeName()<<"\""
        ",value=\""<<filterValueToString(theEntry)<<"\"}"
        "\nin the parameter (sub)list \""<<this->name()<<"\""
        "\nwas not found in the list of valid parameters!"
        "\n\nThe valid parameters and types are:\n"
        <<validParamList.currentParametersString()
        );
      RCP<const ParameterEntryValidator> validator;
      if( (validator=validEntry->validator()).get() ) {
        validator->validateAndModify( entryName, this->name(), &theEntry );
        theEntry.setValidator(validator);
      }
      else {
        const bool validType =
          ( validEntry!=NULL
            ? theEntry.getAny(false).type() == validEntry->getAny(false).type()
            : false
            );
        TEST_FOR_EXCEPTION_PURE_MSG(
          !validType, Exceptions::InvalidParameterType
          ,"Error, the parameter {name=\""<<entryName<<"\","
          "type=\""<<theEntry.getAny(false).typeName()<<"\""
          ",value=\""<<filterValueToString(theEntry)<<"\"}"
          "\nin the parameter (sub)list \""<<this->name()<<"\""
          "\nexists in the list of valid parameters but has the wrong type."
          "\n\nThe correct type is \""
          << validEntry->getAny(false).typeName() << "\"."
          );
        // Note: If there is no validator for this item, then we can not
        // validate the value of the parameter, only its type!
      }
      if( theEntry.isList() && depth > 0 ) {
        sublist_list.push_back(
          ListPlusValidList(
            &getValue<ParameterList>(theEntry),
            &getValue<ParameterList>(*validEntry)
            )
          );
      }
    }
  }
  //
  // Second, loop through the valid parameters at this level that are not set
  // in *this, and set their defaults.
  //
  {
    ConstIterator itr;
    for( itr = validParamList.begin(); itr != validParamList.end(); ++itr ) {
      const std::string  &validEntryName = validParamList.name(itr);
      const ParameterEntry &validEntry = validParamList.entry(itr);
      const ParameterEntry *theEntry = this->getEntryPtr(validEntryName);
      if(!theEntry) {
        // This entry does not exist, so add it.  Here we will only set the
        // value of the entry and its validator and and leave off the
        // documentation.  The reason that the validator is set is so that it
        // can be used to extract and validate entries in the transformed list
        // *this without having to refer back to the valid parameter list.
        ParameterEntry newEntry;
        newEntry.setAnyValue(
          validEntry.getAny(),
          true // isDefault
          );
        newEntry.setValidator(validEntry.validator());
        this->setEntry(validEntryName,newEntry);
      }
    }
  }
  //
  // Now loop through the sublists and validate their parameters and set their
  // defaults!
  //
  for(
    sublist_list_t::iterator sl_itr = sublist_list.begin();
    sl_itr != sublist_list.end();
    ++sl_itr
    )
  {
    if (!sl_itr->validList->disableRecursiveValidation_) {
      sl_itr->list->validateParametersAndSetDefaults(*sl_itr->validList,depth-1);
    }
  }
#ifdef TEUCHOS_PARAMETER_LIST_SHOW_TRACE
  *out << "\n*** Existing ParameterList::validateParametersAndSetDefaults(...) "
    "for this->name()=\""<<this->name()<<"\"...\n";
#endif
}


// private


ParameterList::Iterator ParameterList::nonconstBegin()
{
  return params_.begin();
}


ParameterList::Iterator ParameterList::nonconstEnd()
{
  return params_.end();
}


void ParameterList::updateSubListNames(int depth)
{
  const std::string this_name = this->name();
  Map::iterator itr;
  for( itr = params_.begin(); itr != params_.end(); ++itr ) {
    const std::string &entryName = this->name(itr);
    const ParameterEntry &theEntry = this->entry(itr);
    if(theEntry.isList()) {
      ParameterList &sublistEntry = getValue<ParameterList>(theEntry);
      sublistEntry.setName(this_name+std::string("->")+entryName);
      if(depth > 0)
        sublistEntry.updateSubListNames(depth-1);
    }
  }
}


void ParameterList::validateEntryExists(
  const std::string & /*funcName*/, const std::string &name_in,
  const ParameterEntry *entry_in
  ) const
{
  TEST_FOR_EXCEPTION_PURE_MSG(
    entry_in==NULL, Exceptions::InvalidParameterName
    ,"Error!  The parameter \""<<name_in<<"\" does not exist"\
    "\nin the parameter (sub)list \""<<this->name()<<"\"."
    "\n\nThe current parameters set in (sub)list \""<<this->name()<<"\" are:\n\n"
    << this->currentParametersString()
    );
}


} // namespace Teuchos


bool Teuchos::operator==( const ParameterList& list1, const ParameterList& list2 )
{
  // Check that the top-level names of the two parameter lists are the same
  //const std::string &paramListName1 = list1.name();
  //const std::string &paramListName2 = list2.name();
  //if ( paramListName1 != paramListName2 ) {
  //  return false;
  //}
  ParameterList::ConstIterator itr1, itr2;
  for(
    itr1 = list1.begin(), itr2 = list2.begin();
    itr1 != list1.end() && itr2 != list2.end();
    ++itr1, ++itr2
    )
  {
    const std::string    &entryName1   = list1.name(itr1);
    const std::string    &entryName2   = list2.name(itr2);
    const ParameterEntry &entry1       = list1.entry(itr1);
    const ParameterEntry &entry2       = list2.entry(itr2);
    if( entryName1 != entryName2 ) {
      return false;
    }
    else if( entry1 != entry2 ) {
      return false;
    }
    // Note that the above statement automatically recursively compare the
    // sublists since ParameterList objects are stored in the 'any' variable
    // held by the ParameterEntry object and this same comparison operator will
    // be used.
  }
  // Check that the two parameter lists are the same length:
  if ((itr1 != list1.end()) || (itr2 != list2.end())) {
    return false;
  }
  return true;
}


bool Teuchos::haveSameValues( const ParameterList& list1, const ParameterList& list2 )
{
  // Check that the top-level names of the two parameter lists are the same
  //const std::string &paramListName1 = list1.name();
  //const std::string &paramListName2 = list2.name();
  //if ( paramListName1 != paramListName2 ) {
  //  return false;
  //}
  ParameterList::ConstIterator itr1, itr2;
  for(
    itr1 = list1.begin(), itr2 = list2.begin();
    itr1 != list1.end() && itr2 != list2.end();
    ++itr1, ++itr2
    )
  {
    const std::string    &entryName1   = list1.name(itr1);
    const std::string    &entryName2   = list2.name(itr2);
    const ParameterEntry &entry1       = list1.entry(itr1);
    const ParameterEntry &entry2       = list2.entry(itr2);
    if( entryName1 != entryName2 ) {
      return false;
    }
    if( entry1.isList() && entry2.isList() ) {
      if (
        !haveSameValues(
          getValue<ParameterList>(entry1),
          getValue<ParameterList>(entry2))
        )
      {
        // Note: Above we cast to a non-const ParameterList even through we
        // only need a const ParameterList.  We have to do this since a
        // non-const ParameterList is always added initially which determines
        // the value.
        return false;
      }
    }
    else {
      if( entry1.getAny() != entry2.getAny() ) {
        return false;
      }
    }
  }
  // Check that the two parameter lists are the same length:
  if ((itr1 != list1.end()) || (itr2 != list2.end())) {
    return false;
  }
  return true;
}
