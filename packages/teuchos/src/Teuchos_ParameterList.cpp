// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

//#define TEUCHOS_PARAMETER_LIST_SHOW_TRACE

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_VerboseObject.hpp"


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


// Constructors/Destructor/Info


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


ParameterList::~ParameterList() 
{}


Ordinal ParameterList::numParams() const
{
  return params_.numObjects();
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
    const std::string &name_i = this->name(i);
    const ParameterEntry &entry_i = this->entry(i);
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
    const std::string &name_i = this->name(i);
    const ParameterEntry &entry_i = this->entry(i);
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


void ParameterList::unused(std::ostream& os) const
{
  for (ConstIterator i = this->begin(); i != this->end(); ++i) {
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
    const std::string &entryName = this->name(itr);
    const ParameterEntry &theEntry = this->entry(itr);
    oss
      << "    \""<<entryName<<"\" : "<<theEntry.getAny().typeName()
      <<" = "<<filterValueToString(theEntry) << "\n";
  }
  oss << "  }\n";
  return oss.str();
}


bool ParameterList::isSublist(const std::string& name_in) const
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  const Ordinal param_idx = params_.getObjOrdinalIndex(name_in);
  if (param_idx != SIOVOCB::getInvalidOrdinal()) {
    return params_.getObjPtr(param_idx)->isList();
  }
  return false;
}


bool ParameterList::isParameter(const std::string& name_in) const
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  const Ordinal param_idx = params_.getObjOrdinalIndex(name_in);
  if (param_idx != SIOVOCB::getInvalidOrdinal()) {
    return true;
  }
  return false;
}


bool ParameterList::remove(
  std::string const& name_in, bool throwIfNotExists
  )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  const Ordinal param_idx = params_.getObjOrdinalIndex(name_in);
  if (param_idx != SIOVOCB::getInvalidOrdinal()) {
    // Parameter exists
    params_.removeObj(param_idx);
    return true;
  }
  // Parameter does not exist
  if (throwIfNotExists) {
    validateEntryExists("get", name_in, 0); // Will throw
  }
  return false; // Param does not exist but that is okay
}


ParameterList& ParameterList::sublist(
  const std::string& name_in, bool mustAlreadyExist,
  const std::string& docString
  )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;

  const Ordinal param_idx = params_.getObjOrdinalIndex(name_in);

  Ptr<ParameterEntry> sublist_entry_ptr;

  if (param_idx != SIOVOCB::getInvalidOrdinal()) {
    // Sublist parameter exists
    sublist_entry_ptr = params_.getNonconstObjPtr(param_idx);
    validateEntryIsList(name_in, *sublist_entry_ptr);
  }
  else {
    // Sublist does not exist so we need to create a new one
    validateMissingSublistMustExist(this->name(), name_in, mustAlreadyExist);
    const Ordinal new_param_idx =
      params_.setObj(
        name_in,
        ParameterEntry(
          ParameterList(this->name()+std::string("->")+name_in),
          false,
          true,
          docString
          )
        );
    sublist_entry_ptr = params_.getNonconstObjPtr(new_param_idx);
  }
  
  return any_cast<ParameterList>(sublist_entry_ptr->getAny(false));
}


const ParameterList& ParameterList::sublist(const std::string& name_in) const
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;

  const Ordinal param_idx = params_.getObjOrdinalIndex(name_in);
  if (param_idx == SIOVOCB::getInvalidOrdinal()) {
    validateMissingSublistMustExist(this->name(), name_in, true);
  }

  Ptr<const ParameterEntry>  sublist_entry_ptr = params_.getObjPtr(param_idx);
  validateEntryIsList(name_in, *sublist_entry_ptr);
  
  return any_cast<ParameterList>(sublist_entry_ptr->getAny(false));
}

  
void ParameterList::print() const
{
  this->print(*Teuchos::VerboseObjectBase::getDefaultOStream());
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
  if (this->begin() == this->end()) {
    *out <<"[empty list]" << std::endl;
  }
  else { 
    // Print parameters first
    for (ConstIterator i = this->begin(); i != this->end(); ++i) 
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
      if (showDoc) {
        if (nonnull(validator)) {
          validator->printDoc(docString,OSTab(os).o());
        }
        else if (docString.length()) {
          StrUtils::printLines(OSTab(out).o(),"# ",docString);
        }
      }
    }
    // Print sublists second
    for (ConstIterator i = this->begin(); i != this->end(); ++i) 
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

  
std::ostream& ParameterList::print(std::ostream& os, int indent, bool showTypes, bool showFlags) const
{
  return this->print(os,PrintOptions().indent(indent).showTypes(showTypes).showFlags(showFlags));
}


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
  for (itr = this->begin(); itr != this->end(); ++itr) {
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
    TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(
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
    if (nonnull(validator=validEntry->validator())) {
      validator->validate(theEntry, entryName, this->name()); 
    }
    else {
      const bool validType =
        ( validEntry!=NULL
          ? theEntry.getAny(false).type() == validEntry->getAny(false).type()
          : false
          );
      TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(
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
  // A) loop through and validate the parameters at this level.
  //
  // Here we generate a list of sublists that we will search next
  //
  sublist_list_t sublist_list;
  {
    Iterator itr;
    for (itr = this->nonconstBegin(); itr != this->nonconstEnd(); ++itr) {
      const std::string &entryName = this->name(itr);
      ParameterEntry &theEntry = this->nonconstEntry(itr);
#ifdef TEUCHOS_PARAMETER_LIST_SHOW_TRACE
      OSTab tab(out);
      *out << "\nentryName=\""<<entryName<<"\"\n";
#endif
      const ParameterEntry *validEntry = validParamList.getEntryPtr(entryName);
      TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(
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
      if (nonnull(validator=validEntry->validator())) {
        validator->validateAndModify(entryName, this->name(), &theEntry);
        theEntry.setValidator(validator);
      }
      else {
        const bool validType =
          ( validEntry!=NULL
            ? theEntry.getAny(false).type() == validEntry->getAny(false).type()
            : false
            );
        TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(
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
  // B) Loop through the valid parameters at this level that are not set in
  // *this, and set their defaults.
  //
  {
    ConstIterator itr;
    for (itr = validParamList.begin(); itr != validParamList.end(); ++itr) {
      const std::string &validEntryName = validParamList.name(itr);
      const ParameterEntry &validEntry = validParamList.entry(itr);
      const ParameterEntry *theEntry = this->getEntryPtr(validEntryName);
      if (!theEntry) {
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
  // C) Loop through the sublists and validate their parameters and set their
  // defaults!
  //
  for (
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


void ParameterList::updateSubListNames(int depth)
{
  const std::string this_name = this->name();
  Iterator itr;
  for( itr = this->nonconstBegin(); itr != this->nonconstEnd(); ++itr ) {
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
  TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(
    entry_in==NULL, Exceptions::InvalidParameterName
    ,"Error!  The parameter \""<<name_in<<"\" does not exist"\
    "\nin the parameter (sub)list \""<<this->name()<<"\"."
    "\n\nThe current parameters set in (sub)list \""<<this->name()<<"\" are:\n\n"
    << this->currentParametersString()
    );
}


void ParameterList::validateEntryIsList(
  const std::string &name_in, const ParameterEntry &entry_in
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(
    !entry_in.isList(), Exceptions::InvalidParameterType
    ,"Error, the parameter \"" << name_in << "\" is not a list, it is of type \""
    <<entry_in.getAny(false).typeName()<<"\"!" );
}


void ParameterList::validateMissingSublistMustExist(const std::string &baselist_name,
  const std::string &sublist_name, const bool mustAlreadyExist) const
{
  TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(
    mustAlreadyExist, Exceptions::InvalidParameterName
    ,"The sublist "<<baselist_name<<"->\""<<sublist_name<<"\" does not exist!"
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
