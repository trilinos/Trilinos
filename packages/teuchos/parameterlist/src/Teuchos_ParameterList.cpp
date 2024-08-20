// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//#define TEUCHOS_PARAMETER_LIST_SHOW_TRACE
#include <deque>
#include <functional>
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


ParameterList::ParameterList(const std::string &name_in,
    RCP<const ParameterListModifier> const& modifier_in)
  :name_(name_in), modifier_(modifier_in)
{}


ParameterList::ParameterList(const ParameterList& source)
{
  name_ = source.name_;
  params_ = source.params_;
  disableRecursiveValidation_ = source.disableRecursiveValidation_;
  disableRecursiveModification_= source.disableRecursiveModification_;
  disableRecursiveReconciliation_ = source.disableRecursiveReconciliation_;
  modifier_ = source.modifier_;
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
  disableRecursiveModification_= source.disableRecursiveModification_;
  disableRecursiveReconciliation_ = source.disableRecursiveReconciliation_;
  modifier_ = source.modifier_;
  return *this;
}


void ParameterList::setModifier(
    RCP<const ParameterListModifier> const& modifier_in)
{
  modifier_ = modifier_in;
}


ParameterList& ParameterList::setParameters(const ParameterList& source)
{
  for (ConstIterator i = source.begin(); i != source.end(); ++i) {
    const std::string &name_i = this->name(i);
    const ParameterEntry &entry_i = this->entry(i);
    if (entry_i.isList()) {
      ParameterList &pl = getValue<ParameterList>(entry_i);
      ParameterList &this_pl = this->sublist(name_i, false, entry_i.docString());
      this_pl.setParameters(pl);
      this_pl.setModifier(pl.getModifier());
    } else {
      this->setEntry(name_i, entry_i);
    }
  }
  this->updateSubListNames();
  return *this;
}


ParameterList& ParameterList::setParametersNotAlreadySet(
  const ParameterList& source
  )
{
  for (ConstIterator i = source.begin(); i != source.end(); ++i) {
    const std::string &name_i = this->name(i);
    const ParameterEntry &entry_i = this->entry(i);
    if (entry_i.isList()) {
      ParameterList pl = getValue<ParameterList>(entry_i);
      if (this->isSublist(name_i)){
        this->sublist(name_i, true).setParametersNotAlreadySet(pl);
      } else{
        this->sublist(name_i, pl.getModifier(), entry_i.docString())
            .setParametersNotAlreadySet(pl);
      }
    } else {
      const ParameterEntry *thisEntryPtr = this->getEntryPtr(name_i);
      // If the entry does not already exist, then set it.  Otherwise, leave the
      // existing entry alone.
      if (!thisEntryPtr)
        this->setEntry(name_i, entry_i);
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


ParameterList& ParameterList::disableRecursiveModification()
{
  disableRecursiveModification_ = true;
  return *this;
}


ParameterList& ParameterList::disableRecursiveReconciliation()
{
  disableRecursiveReconciliation_ = true;
  return *this;
}


ParameterList& ParameterList::disableRecursiveAll()
{
  this->disableRecursiveModification();
  this->disableRecursiveValidation();
  this->disableRecursiveReconciliation();
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


ParameterList& ParameterList::sublist(
  const std::string& name_in, RCP<const ParameterListModifier> const& modifier_in,
  const std::string& docString
  )
{
  bool alreadyExists = this->isParameter(name_in);
  TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(
    alreadyExists, Exceptions::InvalidParameterName
    ,"The parameter "<<this->name()<<"->\""<<name_in<<"\" already exists."
    );
  ParameterList &subpl = this->sublist(name_in, false, docString);
  subpl.setModifier(modifier_in);
  return subpl;
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
  const int   indent      = printOptions.indent();
  const bool  showTypes   = printOptions.showTypes();
  const bool  showFlags   = printOptions.showFlags();
  const bool  showDoc     = printOptions.showDoc();
  const bool  showDefault = printOptions.showDefault();
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
      if(!showDefault && entry_i.isDefault())
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


std::ostream& ParameterList::print(std::ostream& os, int indent, bool showTypes, bool showFlags, bool showDefault) const
{
  return this->print(os,PrintOptions().indent(indent).showTypes(showTypes).showFlags(showFlags).showDefault(showDefault));
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


void ParameterList::modifyParameterList(ParameterList & valid_pl,
    int const depth)
{
  RCP<const ParameterListModifier> modifier;
  if (nonnull(modifier = valid_pl.getModifier())) {
    modifier->modify(*this, valid_pl);
    this->setModifier(modifier);
  }
  ConstIterator itr;
  for (itr = valid_pl.begin(); itr != valid_pl.end(); ++itr){
    const std::string &entry_name = itr->first;
    const ParameterEntry &cur_entry = itr->second;
    if (cur_entry.isList() && depth > 0){
      ParameterList &valid_pl_sublist = valid_pl.sublist(entry_name, true);
      if(!valid_pl_sublist.disableRecursiveModification_){
        const ParameterEntry *validEntry = this->getEntryPtr(entry_name);
        TEUCHOS_TEST_FOR_EXCEPTION(
          !validEntry, Exceptions::InvalidParameterName
          ,"Error, the parameter {name=\""<<entry_name<<"\","
          "type=\""<<cur_entry.getAny(false).typeName()<<"\""
          ",value=\""<<filterValueToString(cur_entry)<<"\"}"
          "\nin the parameter (sub)list \""<<this->name()<<"\""
          "\nwas not found in the list of parameters during modification."
          "\n\nThe parameters and types are:\n"
          <<this->currentParametersString()
          );
        ParameterList &pl_sublist = this->sublist(entry_name, true);
        pl_sublist.modifyParameterList(valid_pl_sublist, depth-1);
      }
    }
  }
}


void ParameterList::reconcileParameterList(ParameterList & valid_pl,
    const bool left_to_right)
{
  // We do a breadth-first traversal of `valid_pl` and store references to all of the sublists
  // in `valid_pl` in a deque with a matching deque for `this`.
  std::deque<std::reference_wrapper<ParameterList>> refs, valid_refs, tmp, valid_tmp;
  tmp.push_back(*this);
  valid_tmp.push_back(valid_pl);
  while (!valid_tmp.empty()){
    ParameterList &cur_node = tmp.front();
    ParameterList &valid_cur_node = valid_tmp.front();
    tmp.pop_front();
    valid_tmp.pop_front();
    refs.push_back(cur_node);
    valid_refs.push_back(valid_cur_node);
    // Look for all sublists in valid_tmp
    for (auto itr = valid_cur_node.begin(); itr != valid_cur_node.end(); ++itr){
      const std::string &entry_name = itr->first;
      if (valid_cur_node.isSublist(entry_name)){
        const ParameterEntry &cur_entry = itr->second;
        ParameterList &valid_cur_node_sublist = valid_cur_node.sublist(entry_name);
        if (!valid_cur_node_sublist.disableRecursiveReconciliation_){
          TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(
            !cur_node.isSublist(entry_name), Exceptions::InvalidParameterName
            ,"Error, the parameter {name=\"" << entry_name <<"\","
            "type=\"" << cur_entry.getAny(false).typeName() << "\""
            ",value=\"" << filterValueToString(cur_entry) << "\"}"
            "\nin the parameter (sub)list \"" <<cur_node.name() << "\""
            "\nwas not found in the list of parameters during reconciliation."
            "\n\nThe parameters and types are:\n"
            <<cur_node.currentParametersString()
            );
          if (left_to_right){
            valid_tmp.push_back(valid_cur_node_sublist);
            tmp.push_back(cur_node.sublist(entry_name));
          } else{
            valid_tmp.push_front(valid_cur_node_sublist);
            tmp.push_front(cur_node.sublist(entry_name));
          }
        }
      }
    }
  }
  // We now apply the reconciliation from the bottom to the top of the parameter lists by
  // traversing the deques from the back to the front.
  RCP<const ParameterListModifier> modifier;
  std::deque<std::reference_wrapper<ParameterList>>::reverse_iterator ref, valid_ref;
  for(ref = refs.rbegin(), valid_ref = valid_refs.rbegin();
      ref != refs.rend() && valid_ref != valid_refs.rend();
      ++ref, ++valid_ref){
    if (nonnull(modifier = valid_ref->get().getModifier())) {
      modifier->reconcile(ref->get());
    }
  }
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
        RCP<const ParameterEntryValidator> validator;
        if (nonnull(validator=validEntry.validator())) {
#if defined(HAVE_TEUCHOS_MODIFY_DEFAULTS_DURING_VALIDATION)
          validEntry.validator()->validateAndModify(this->name(itr), validEntryName, &newEntry);
          // validateAndModify changes the default status so we reset it
          newEntry.setAnyValue(newEntry.getAny(), true);
#endif
          newEntry.setValidator(validator);
        }
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
  if (!Teuchos::haveSameModifiers(list1, list2)){
    return false;
  }
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
    // Note that the above statement automatically recursively compares the
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


bool Teuchos::haveSameModifiers(const ParameterList &list1, const ParameterList &list2) {
  // Check that the modifiers are the same
  const RCP<const ParameterListModifier> &modifier1 = list1.getModifier();
  const RCP<const ParameterListModifier> &modifier2 = list2.getModifier();
  // Compare the modifiers.
  const bool modifier1_is_null = is_null(modifier1);
  const bool modifier2_is_null = is_null(modifier2);
  if( modifier1_is_null || modifier2_is_null ){
    if ( modifier1_is_null != modifier2_is_null ){
      return false;
    }
  } else if ( *modifier1 != *modifier2 ){
    return false;
  }
  // Now look for more sublists
  ParameterList::ConstIterator itr1, itr2;
  for(
    itr1 = list1.begin(), itr2 = list2.begin();
    itr1 != list1.end() && itr2 != list2.end();
    ++itr1, ++itr2
    )
  {
    // Check the modifiers in each sublist.
    const ParameterEntry &entry1 = itr1->second;
    const ParameterEntry &entry2 = itr2->second;
    if (entry1.isList() && entry2.isList()){
      if ( !haveSameModifiers( Teuchos::getValue<ParameterList>(entry1),
                               Teuchos::getValue<ParameterList>(entry2) ) ){
        return false;
      }
    }
  }
  return true;
}


bool Teuchos::haveSameValues( const ParameterList& list1, const ParameterList& list2, bool verbose )
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
      if (verbose) std::cerr << "entryName1 \"" << entryName1 << "\" != entryName2 \"" << entryName2 << "\"\n";
      return false;
    }
    if( entry1.isList() && entry2.isList() ) {
      if (
        !haveSameValues(
          getValue<ParameterList>(entry1),
          getValue<ParameterList>(entry2),
          verbose)
        )
      {
        // Note: Above we cast to a non-const ParameterList even through we
        // only need a const ParameterList.  We have to do this since a
        // non-const ParameterList is always added initially which determines
        // the value.
        if (verbose) std::cerr << "sublists \"" << entryName1 << "\" differ\n";
        return false;
      }
    }
    else {
      if( entry1.getAny() != entry2.getAny() ) {
        if (verbose) std::cerr << "for key \"" << entryName1 << "\", value \"" << entry1.getAny() << "\" != \"" << entry2.getAny() << "\"\n";
        return false;
      }
    }
  }
  // Check that the two parameter lists are the same length:
  if ((itr1 != list1.end()) || (itr2 != list2.end())) {
    if (verbose) std::cerr << "lists are not the same size\n";
    return false;
  }
  return true;
}


bool Teuchos::haveSameValuesSorted( const ParameterList& list1, const ParameterList& list2, bool verbose )
{
  // Check that the top-level names of the two parameter lists are the same
  //const std::string &paramListName1 = list1.name();
  //const std::string &paramListName2 = list2.name();
  //if ( paramListName1 != paramListName2 ) {
  //  return false;
  //}
  ParameterList::ConstIterator itr1, itr2;
  Array<std::string> arr1, arr2;
  for(itr1 = list1.begin(); itr1 != list1.end(); ++itr1){
    arr1.push_back(list1.name(itr1));
  }
  for(itr2 = list2.begin(); itr2 != list2.end(); ++itr2){
    arr2.push_back(list2.name(itr2));
  }
  // Check that the two parameter lists are the same length:
  if (arr1.size() != arr2.size()) {
    if (verbose) std::cerr << "lists are not the same size\n";
    return false;
  }
  std::sort(arr1.begin(), arr1.end());
  std::sort(arr2.begin(), arr2.end());
  Array<std::string>::iterator iarr1, iarr2;
  for(
    iarr1 = arr1.begin(), iarr2 = arr2.begin();
    iarr1 != arr1.end() && iarr2 != arr2.end();
    ++iarr1, ++iarr2
    )
  {
    const std::string    &entryName1   = *iarr1;
    const std::string    &entryName2   = *iarr2;
    const ParameterEntry &entry1       = list1.getEntry(entryName1);
    const ParameterEntry &entry2       = list2.getEntry(entryName2);
    if( entryName1 != entryName2 ) {
      if (verbose) std::cerr << "entryName1 \"" << entryName1 << "\" != entryName2 \"" << entryName2 << "\"\n";
      return false;
    }
    if( entry1.isList() && entry2.isList() ) {
      if (
        !haveSameValuesSorted(
          getValue<ParameterList>(entry1),
          getValue<ParameterList>(entry2),
          verbose)
        )
      {
        // Note: Above we cast to a non-const ParameterList even through we
        // only need a const ParameterList.  We have to do this since a
        // non-const ParameterList is always added initially which determines
        // the value.
        if (verbose) std::cerr << "sublists \"" << entryName1 << "\" differ\n";
        return false;
      }
    }
    else {
      if( entry1.getAny() != entry2.getAny() ) {
        if (verbose) std::cerr << "for key \"" << entryName1 << "\", value \"" << entry1.getAny() << "\" != \"" << entry2.getAny() << "\"\n";
        return false;
      }
    }
  }
  return true;
}
