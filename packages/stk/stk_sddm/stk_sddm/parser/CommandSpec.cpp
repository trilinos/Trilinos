/*--------------------------------------------------------------------*/
/*    Copyright 2000 - 2010 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
/*        Copyright 2000 Sandia Corporation, Albuquerque, NM.         */
/*--------------------------------------------------------------------*/

#include <string>
#include <cstring>
#include <algorithm>

#include <parser/Prsr_CommandSpec.h>
#include <parser/Prsr_Parser.h>
#include <parser/Prsr_DiagWriter.h>
#include <parser/Prsr_SDDMLexer.h>
#include <stk_util/parallel/Exception.hpp>

using sierra::Diag::dendl;
using sierra::Diag::pop;
using sierra::Diag::push;

namespace sierra {
namespace Prsr {

std::ostream &
printCommandSpecForTaxonomy(
  std::ostream &        os,
  const std::string &   key_id)
{
  Identifier identifier = Identifier::createIdentifier(key_id);
  
  Parser::db_command_spec()[identifier]->printSpecification(os, true);

  return os;
}

namespace {

std::string
key(
  const String &        keywords) 
{
  std::string s;

  if (!*keywords.c_str()) {
    s = "EMPTY";
    return s;
  }

  for (const char *c = keywords.c_str(); *c; ++c) {
    if (c == keywords.c_str())
      s += std::toupper(*c);
    else if ((*c == ' ' || *c == '_' || *c == '-' || *c == '|') && *(c + 1)) 
      s += std::toupper(*++c);
    else if (*c == ':')
      ;
    else
      s += std::tolower(*c);
  }

  return s;
}


stk::sddm::AnyType *
type_of(
  const CommandSpec &                           command_spec, 
  const Lexer::LexemVector::const_iterator      it)
{
  String name = (*it).getString();

  stk::sddm::AnyType *any_type = 0;

  const CommandSpec::ParamSpecList &param_spec_list = command_spec.getParamSpecList();
  const CommandSpec::ParamSpec *param_spec = 0;
  CommandSpec::ParamSpecList::const_iterator param_it;
  for (param_it = param_spec_list.begin(); param_it != param_spec_list.end(); ++param_it)
    if ((*param_it).getId() != -1)
      break;
  
  for (; param_it != param_spec_list.end(); ++param_it) {
    if ((*param_it).getName() == name.c_str()) {
      param_spec = &(*param_it);
      break;
    }
  }

  if (param_spec) {
    switch (param_spec->getType()) {
    case CommandSpec::LINE_PARAMETER:
      any_type =  new stk::sddm::Type<std::string>();
      break;
      
    case CommandSpec::ENUMERATED:
      if (param_spec->getCount().second == 1)
        any_type = new stk::sddm::Type<stk::sddm::Enumerand>();
      else {
        any_type = new stk::sddm::Value<std::vector<stk::sddm::Enumerand> >();
      }
        
      break;
      
    case CommandSpec::INTEGER:
      if (param_spec->getCount().second == 1)
        any_type =  new stk::sddm::Type<int>();
      else
        any_type =  new stk::sddm::Type<std::vector<int> >();
      break;
      
    case CommandSpec::REAL:
      if (param_spec->getCount().second == 1)
        any_type = new stk::sddm::Type<double>();
      else
        any_type = new stk::sddm::Type<std::vector<double> >();
      break;
      
    case CommandSpec::STRING:
    case CommandSpec::USTRING:              
    case CommandSpec::QUOTED_STRING:
    case CommandSpec::EXPRESSION:
      if (param_spec->getCount().second == 1)
        any_type = new stk::sddm::Type<std::string>();
      else 
        any_type =  new stk::sddm::Type<std::vector<std::string> >();
      break;

    default: {    
        RuntimeWarning x;
    
        x << "Cannot process parameter " << command_spec.identifier() << "." << name << " for sddmid creation" << std::endl
          << "Parameters are: ";
    
        for (CommandSpec::ParamSpecList::const_iterator param_it = param_spec_list.begin(); param_it != param_spec_list.end(); ++param_it) 
          x << (*param_it).getName() << ", ";
        RuntimeWarning() << "Cannot process parameter " << command_spec.identifier() << "." << name << " of type " << param_spec->getType() << " for sddmid creation.  " << std::endl;
      }
      
      break;            
    }
  }
  else {
    RuntimeWarning x;
    
    x << "Cannot find parameter " << name << " for sddmid creation.  " << std::endl
      << "Candidates are: ";
    
    for (CommandSpec::ParamSpecList::const_iterator param_it = param_spec_list.begin(); param_it != param_spec_list.end(); ++param_it) 
      x << (*param_it).getName() << ", ";
  }

  prsrout << "Type of " << name << " is ";
  if (param_spec)
    prsrout << param_spec->getType() << " ";
  if (any_type)
    prsrout << any_type->demangle();
  prsrout << dendl;
  
  return any_type;
}

} // namespace <unnamed>


std::ostream &
CommandSpec::Count::printSpecification(
  std::ostream &		os) const
{
  /* %TRACE% */ Traceback trace__("sierra::Prsr::CommandSpec::Count::printSpecification( std::ostream & os) const"); /* %TRACE% */
  if (first == 0 && second == UInt_MAX)
    os << "*";
  else if (first == 1 && second == UInt_MAX)
    os << "+";
  else if (first > 1 || second > 1) {
    os << "[" << first;
    if (second != first) {
      os << ":";
      if (second != UInt_MAX)
	os << second;
    }

    os << "]";
  }
  return os;
}


CommandSpec::ParamSpec::ParamSpec(
  ParamType		type,
  int			id,
  const String &	name)
  : m_type(type),
    m_optional(false),
    m_id(id),
    m_name(key(name)),
    m_token(type == TOKEN ? name : ""),
    m_count(0, 0),
    m_enumeration(0),
    m_intRange(),
    m_realRange(),
    m_addTagList(),
    m_inTagList()
{ /* %TRACE% */ Traceback trace__("sierra::Prsr::CommandSpec::ParamSpec::ParamSpec( ParamType type, int id, const String & name)"); /* %TRACE% */}


CommandSpec::ParamSpec::ParamSpec(
  ParamType		type,
  int			id,
  const String &	name,
  const Count &	        count,
  const Range<Real> &	range)
  : m_type(type),
    m_optional(false),
    m_id(id),
    m_name(key(name)),
    m_token(),
    m_count(count),
    m_enumeration(0),
    m_intRange(),
    m_realRange(range),
    m_addTagList(),
    m_inTagList()
{
  /* %TRACE% */ Traceback trace__("sierra::Prsr::CommandSpec::ParamSpec::ParamSpec( ParamType type, int id, const String & name, const Count & count, const Range<Real> & range)"); /* %TRACE% */
}


CommandSpec::ParamSpec::ParamSpec(
  ParamType		type,
  int			id,
  const String &	name,
  const Count &	        count,
  const Range<Int> &	range)
  : m_type(type),
    m_optional(false),
    m_id(id),
    m_name(key(name)),
    m_token(),
    m_count(count),
    m_enumeration(0),
    m_intRange(range),
    m_realRange(),
    m_addTagList(),
    m_inTagList()
{
  /* %TRACE% */ Traceback trace__("sierra::Prsr::CommandSpec::ParamSpec::ParamSpec( ParamType type, int id, const String & name, const Count & count, const Range<Int> & range)"); /* %TRACE% */
}


CommandSpec::ParamSpec::ParamSpec(
  ParamType		type,
  int			id,
  const String &	name,
  const Count &	        count ,
  const Enumeration *	enumeration)
  : m_type(type),
    m_optional(false),
    m_id(id),
    m_name(key(name)),
    m_token(),
    m_count(count),
    m_enumeration(enumeration),
    m_intRange(),
    m_realRange(),
    m_addTagList(),
    m_inTagList()
{
  /* %TRACE% */ Traceback trace__("sierra::Prsr::CommandSpec::ParamSpec::ParamSpec( ParamType type, int id, const String & name, const Count & count , const Enumeration * enumerated)"); /* %TRACE% */
}


CommandSpec::ParamSpec::ParamSpec(
  ParamType		type,
  int			id,
  const String &	name,
  const Count &	        count)
  : m_type(type),
    m_optional(false),
    m_id(id),
    m_name(key(name)),
    m_token(),
    m_count(count),
    m_enumeration(0),
    m_intRange(),
    m_realRange(),
    m_addTagList(),
    m_inTagList()
{
  /* %TRACE% */ Traceback trace__("sierra::Prsr::CommandSpec::ParamSpec::ParamSpec( ParamType type, int id, const String & name, const Count & count)"); /* %TRACE% */
}


// bool
// CommandSpec::ParamSpec::matches(
//   const Lexeme &	lexeme) const
// {
//   switch (m_type) {
//     case ENDLINE:
//       return lexeme.getToken() == TOKEN_ENDLINE;
//     case TOKEN:
//       return lexeme.getToken() == TOKEN_IDENTIFIER && lexeme.getValue() == m_token;
//     case ENUMERATED:
//       return lexeme.getToken() == TOKEN_IDENTIFIER
//	&& m_enumeration->getEnumerand(lexeme.getValue()) != m_enumeration->end();
//     case INTEGER:
//       return lexeme.getToken() == TOKEN_INTEGER_LITERAL;
//     case REAL:
//       return lexeme.getToken() == TOKEN_REAL_LITERAL
//	|| lexeme.getToken() == TOKEN_INTEGER_LITERAL;
//     case STRING:
//     case USTRING:
//     case QUOTED_STRING:
//       return lexeme.getToken() == TOKEN_IDENTIFIER
//	|| lexeme.getToken() == TOKEN_STRING_LITERAL;
//     default:
//       return false;
//   }
// }


// bool
// CommandSpec::ParamSpec::terminates(
//   const Lexeme &	lexeme) const
// {
//   switch (m_type) {
//     case ENDLINE:
//       return false;
//     case TOKEN:
//       return lexeme.getToken() == TOKEN_IDENTIFIER && lexeme.getValue() == m_token;
//     case ENUMERATED:
//       return lexeme.getToken() == TOKEN_IDENTIFIER
//	&& m_enumeration->getEnumerand(lexeme.getValue()) != m_enumeration->end();
//     case INTEGER:
//       return false;
//     case REAL:
//       return false;
//     case STRING:
//     case USTRING:
//     case QUOTED_STRING:
//       return false;
//     default:
//       return false;
//   }
// }


template<class T>
std::ostream &
CommandSpec::Range<T>::printSpecification(
  std::ostream &		os) const
{
  /* %TRACE% */ Traceback trace__("sierra::Prsr::CommandSpec::Range<T>::printSpecification( std::ostream & os) const"); /* %TRACE% */
  if ((m_flags & (LOWER_BOUNDED | LOWER_INCLUSIVE | UPPER_BOUNDED | UPPER_INCLUSIVE))
      == (LOWER_BOUNDED | LOWER_INCLUSIVE | UPPER_BOUNDED | UPPER_INCLUSIVE)) {
    os << "(" << m_lowerBound << " - " << m_upperBound << ")";
  }
  else if (m_flags & (LOWER_BOUNDED | UPPER_BOUNDED)) {
    os << "(";
    if (m_flags & LOWER_BOUNDED) {
      os << (m_flags & LOWER_INCLUSIVE ? ">=" : ">") << m_lowerBound;
      if (m_flags & UPPER_BOUNDED)
	os << " ";
    }
    if (m_flags & UPPER_BOUNDED)
      os << (m_flags & UPPER_INCLUSIVE ? "<=" : "<") << m_upperBound;
    os << ")";
  }

  return os;
}


std::ostream &
CommandSpec::ParamSpec::printSpecification(
  std::ostream &		os,
  bool                          print_triplets) const
{
  /* %TRACE% */ Traceback trace__("sierra::Prsr::CommandSpec::ParamSpec::printSpecification( std::ostream & os) const"); /* %TRACE% */
  if (print_triplets)
    switch (m_type) {
    case ENUMERATED:
      if (!m_enumeration->isGeneral())
        os << specification(*m_enumeration) << " ";
      break;
    default:
      break;
    }
  else
    switch (m_type) {
    case TOKEN:
      os << title(m_token);
      break;
    case ENUMERATED:
      os << specification(*m_enumeration);
      break;
    case INTEGER:
      os << "<" << lower(m_name) << ": integer";
      m_intRange.printSpecification(os);
      os << specification(m_count) << ">";
      break;
    case REAL:
      os << "<" << lower(m_name) << ": real";
      m_realRange.printSpecification(os);
      os << specification(m_count) << ">";
      break;
    case STRING:
    case QUOTED_STRING:
    case USTRING:
      os << "<" << lower(m_name) << ": string" << specification(m_count) << ">";
      break;
    case EXPRESSION:
      os << "<" << lower(m_name) << ": expression" << specification(m_count) << ">";
      break;
    case LINE_PARAMETER:
      os << "<" << lower(m_name) << ">";
      break;
    case ENDLINE:
      break;
    default:
      break;
    }

  return os;
}


CommandSpec::CommandSpec(
  const Identifier &	id,
  bool			block_command)
  : m_identifier(id),
    m_blockCommand(block_command),
    m_requiredParameters(-1),
    m_optionalParameters(-1),
    m_summary("No summary provided"),
    m_description("No description provided"),
    m_keywords(),
    m_sddmkey(0),
    m_sddmid(),
    m_paramSpecList()
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::CommandSpec( const Identifier & id, bool block_command)"); /* %TRACE% */
}


CommandSpec::~CommandSpec()
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::~CommandSpec()"); /* %TRACE% */
}


const CommandSpec::ParamSpec &
CommandSpec::getParamSpec(
  int			id) const
{
  /* %TRACE% */ Traceback trace__("sierra::Prsr::CommandSpec::getParamSpec( int id) const"); /* %TRACE% */
  for (ParamSpecList::const_iterator it = m_paramSpecList.begin(); it != m_paramSpecList.end(); ++it) {
    if ((*it).getId() == id)
      return *it;
  }
  throw RuntimeError() << "No " << id << " parameter exists for command " << m_keywords;
}

std::vector<Identifier>
CommandSpec::s_activeIdentifiers;

std::string
active_identifiers() 
{
  std::ostringstream s;
  
  for (std::vector<Identifier>::const_iterator it = CommandSpec::s_activeIdentifiers.begin(); it != CommandSpec::s_activeIdentifiers.end(); ++it)
    s << (*it) << " ";

  return s.str();
}


std::ostream &
CommandSpec::printSpecification(
  std::ostream &		os,
  bool                          print_triplets) const
{
  /* %TRACE% */ Traceback trace__("sierra::Prsr::CommandSpec::printSpecification( std::ostream & os) const"); /* %TRACE% */
  if (print_triplets) {
    os << active_identifiers() << " ";
    
    bool optional = false;
    for (ParamSpecList::const_iterator it = m_paramSpecList.begin(); it != m_paramSpecList.end(); ++it) {
      if ((*it).getOptional() && !optional) {
        optional = true;
        os << "[ ";
      }
      os << specification(*it, true);
    }
    if (optional)
      os << "]";
    os << " , ";
  }
  
  if (m_blockCommand)
    os << "Begin ";

  bool optional = false;
  for (ParamSpecList::const_iterator it = m_paramSpecList.begin(); it != m_paramSpecList.end(); ++it) {
    if (it != m_paramSpecList.begin())
      os << " ";
    if ((*it).getOptional() && !optional) {
      optional = true;
      os << "[ ";
    }
    
    os << specification(*it, false);
  }
  if (optional)
    os << "]";
  
//  os << "  " << word_wrap(m_summary, 60, "  ", "  ");

  return os;
}


std::ostream &
CommandSpec::printTaxonomy(
  std::ostream &        os) const
{
  /* %TRACE% */ Traceback trace__("sierra::Prsr::CommandSpec::printTaxonomy( std::ostream & os) const"); /* %TRACE% */
//  os << m_identifier << " ";
  std::string s(m_sddmid);

  std::string::size_type i = s.rfind(';');
  if (i != std::string::npos)
    s.replace(0, i + 1, "");
  if (m_blockCommand) {  
    i = s.rfind(':');
    if (i != std::string::npos)
      s.replace(i, std::string::npos, "");
  }
  
  os << s;

  return os;
}


bool
CommandSpec::set_sddmid(
  const char * const    sddmid)
{
//   if (!sddmid || sddmid[0] =='\0') {  
//     if (keywords()[0]) {
//       if (number_required_parameters() > 1
//           && getParamSpec(0).getType() == ENUMERATED
//           && getParamSpec(0).getEnumeration().getName() == "Delimiter") 
//         m_sddmid = "$keywords:$1";
//       else
//         m_sddmid = "$keywords:$0";
//     }
//     else
//       m_sddmid = "$0";
//   }
//   else
//     m_sddmid = sddmid;
  
  if (sddmid && sddmid[0] != '\0')
    m_sddmid = sddmid;
  else if (m_sddmid.empty() && number_required_parameters() > 0)
    m_sddmid = "$0";
  
  return true;
}


bool
CommandSpec::set_keywords(
  const char * const    keywords)
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::set_keywords( const char * const keywords )"); /* %TRACE% */
  if (!keywords)
    return false;

  const char * k = keywords ;

  // Eliminate leading spaces
  for ( ; isspace(*k); ++k )
    ;

  m_keywords.clear();
  m_keywords.reserve(std::strlen(k));

  for ( m_keywords.push_back(*k) ; *k != 0 ; )
    if (!isspace(*++k) || !isspace(m_keywords[m_keywords.length() - 1]))
      m_keywords.push_back(*k);

  const char *t0=NULL;
  const char *t1=NULL;

  for (t0 = keywords; *t0; ) {
    while (*t0 && isspace(*t0))
      t0++;
    if (*t0) {
      t1 = t0 + 1;
      while (*t1 && !isspace(*t1))
	t1++;
      m_paramSpecList.addParamSpec(ParamSpec(TOKEN, -1, String(t0, t1)));
      t0 = t1;
    }
  }

  m_sddmid = key(m_keywords.c_str());

  return true;
}


bool
CommandSpec::set_sddm_key(
  const char * const	keywords)
{
  if (!keywords)
    return false;

  if (keywords[0] == '/') {
    m_sddmkey = Parser::db_taxonomy().createKey(identifier().to_string(), key(&keywords[1]));
    Parser::db_taxonomy().makeRoot(m_sddmkey);
  }
  else
    m_sddmkey = Parser::db_taxonomy().createKey(identifier().to_string(), key(keywords));

  return true;
}

bool
CommandSpec::set_summary(
  const char * const	text )
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::set_summary( const char * const text )"); /* %TRACE% */

  if (!text)
    return false;

  m_summary = text;
  return true;
}


bool
CommandSpec::set_description(
  const char * const	text )
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::set_description( const char * const text )"); /* %TRACE% */

  if (!text)
    return false;

  m_description = text;
  return true;
}

//----------------------------------------------------------------------

bool
CommandSpec::set_line_parameter(
  const char *		name)
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::set_line_parameter( const char * name)"); /* %TRACE% */
  if ( -1 != m_optionalParameters ) return false ;

  m_requiredParameters = 1 ;
  m_optionalParameters = 0 ;

  m_paramSpecList.addParamSpec(ParamSpec(LINE_PARAMETER, 0, name));

  if (m_blockCommand) {
    m_sddmid = key(m_keywords.c_str()) + "/" + key(name) + ":$" + key(name);
  }
  else {
    m_sddmid = key(m_keywords.c_str()) + ":$" + key(name);
  }
  
  return true;
}

bool
CommandSpec::sddm_cleanup()
{
  return true;
}


bool CommandSpec::set_identifier()
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::set_identifier()"); /* %TRACE% */
  if ( -1 != m_optionalParameters ) return false ;

  m_requiredParameters = 1 ;
  m_optionalParameters = 0 ;

  m_paramSpecList.addParamSpec(ParamSpec(IDENTIFIER, 0, "<identifier>"));

  return true ;
}

bool
CommandSpec::set_number_of_parameters(
  const UInt required ,
  const UInt optional )
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::set_number_of_parameters( const UInt required , const UInt optional )"); /* %TRACE% */
  if ( -1 != m_optionalParameters ) return false ;

  m_requiredParameters  = required ;
  m_optionalParameters  = optional ;

  return true ;
}


bool
CommandSpec::set_parameter_optional(
  const UInt		ith_parameter ,
  bool			optional)
{
  /* %TRACE% */ Traceback trace__("sierra::Prsr::CommandSpec::set_parameter_optional( const UInt ith_parameter , bool optional)"); /* %TRACE% */
  if (optional) {
    if (m_optionalParameters == -1) {
      m_requiredParameters = ith_parameter;
      m_optionalParameters = 1;
    }
    else
      ++m_optionalParameters;

    for (ParamSpecList::iterator it = m_paramSpecList.begin(); it != m_paramSpecList.end(); ++it)
      if ((UInt)(*it).getId() == ith_parameter)
	(*it).setOptional(optional);
  }

  return true;
}

bool CommandSpec::set_parameter_token(
  const UInt ith_parameter ,
  const char * const token )
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::set_parameter_token( const UInt ith_parameter , const char * const token )"); /* %TRACE% */

  const char *t0=NULL;
  const char *t1=NULL;

  for (t0 = token; *t0; ) {
    while (*t0 && isspace(*t0))
      t0++;
    if (*t0) {
      t1 = t0 + 1;
      while (*t1 && !isspace(*t1))
        t1++;
      m_paramSpecList.addParamSpec(ParamSpec(TOKEN, ith_parameter, String(t0, t1)));
      t0 = t1;
    }
  }

  return true ;
}

bool CommandSpec::set_parameter_real(
  const UInt		ith_parameter ,
  const char *		name,
  const Count &	        count,
  const Range<Real> &	range)
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::set_parameter_real( const UInt ith_parameter , const char * name, const Count & count, const Range<Real> & range)"); /* %TRACE% */
  m_paramSpecList.addParamSpec(ParamSpec(REAL, ith_parameter, name, count, range));

  if (ith_parameter == 0) {
    m_sddmid = key(m_keywords.c_str()) + ":$" + key(name);
  }

  else if (ith_parameter == 1 && getParamSpec(0).getType() == ENUMERATED && (getParamSpec(0).getName() == "Delimiter" || getParamSpec(0).getName() == "Assign")) {
    m_sddmid = key(m_keywords.c_str()) + ":$" + key(name);
  }
  
  return true ;
}

bool CommandSpec::set_parameter_int(
  const UInt		ith_parameter ,
  const char *		name,
  const Count &	        count ,
  const Range<Int> &	range)
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::set_parameter_int( const UInt ith_parameter , const char * name, const Count & count , const Range<Int> & range)"); /* %TRACE% */
  m_paramSpecList.addParamSpec(ParamSpec(INTEGER, ith_parameter, name, count, range));

  if (ith_parameter == 0) {
    m_sddmid = key(m_keywords.c_str()) + ":$" + key(name);
  }
  
  else if (ith_parameter == 1 && getParamSpec(0).getType() == ENUMERATED && (getParamSpec(0).getName() == "Delimiter" || getParamSpec(0).getName() == "Assign")) {
    m_sddmid = key(m_keywords.c_str()) + ":$" + key(name);
  }
  
  return true ;
}


bool CommandSpec::set_parameter_enum(
  const UInt		ith_parameter ,
  const char *		name,
  const Enumeration *	enumeration ,
  const Count &	        count )
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::set_parameter_enum( const UInt ith_parameter , const char * name, const Enumeration * enumeration , const Count & count )"); /* %TRACE% */
  m_paramSpecList.addParamSpec(ParamSpec(ENUMERATED, ith_parameter, name, count, enumeration));

  if (ith_parameter == 0 && String(name) != "Delimiter") {
    m_sddmid = key(m_keywords.c_str()) + ":$" + key(name);
  }

  else if (ith_parameter == 1 && getParamSpec(0).getType() == ENUMERATED && (getParamSpec(0).getName() == "Delimiter" || getParamSpec(0).getName() == "Assign")) {
    m_sddmid = key(m_keywords.c_str()) + ":$" + key(name);
  }
  
  return true ;
}

bool CommandSpec::set_parameter_string(
  const UInt		ith_parameter ,
  const char *		name,
  const Count &	        count )
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::set_parameter_string( const UInt ith_parameter , const char * name, const Count & count )"); /* %TRACE% */
  m_paramSpecList.addParamSpec(ParamSpec(STRING, ith_parameter, name, count));

  if (ith_parameter == 0)
    m_sddmid = key(m_keywords.c_str()) + ":$" + key(name);

  else if (ith_parameter == 1 && getParamSpec(0).getType() == ENUMERATED && (getParamSpec(0).getName() == "Delimiter" || getParamSpec(0).getName() == "Assign"))
    m_sddmid = key(m_keywords.c_str()) + ":$" + key(name);

  return true ;
}

bool CommandSpec::set_parameter_ustring(
  const UInt		ith_parameter ,
  const char *		name,
  const Count &	        count )
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::set_parameter_ustring( const UInt ith_parameter , const char * name, const Count & count )"); /* %TRACE% */
  m_paramSpecList.addParamSpec(ParamSpec(USTRING, ith_parameter, name, count));

  if (ith_parameter == 0) {
    m_sddmid = key(m_keywords.c_str()) + ":$" + key(name);
  }
  
  else if (ith_parameter == 1 && getParamSpec(0).getType() == ENUMERATED && (getParamSpec(0).getName() == "Delimiter" || getParamSpec(0).getName() == "Assign")) {
    m_sddmid = key(m_keywords.c_str()) + ":$" + key(name);
  }
  
  return true ;
}

bool CommandSpec::set_parameter_quoted_string(
  const UInt		ith_parameter ,
  const char *		name,
  const Count &	        count )
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::set_parameter_quoted_string( const UInt ith_parameter , const char * name, const Count & count )"); /* %TRACE% */
  m_paramSpecList.addParamSpec(ParamSpec(QUOTED_STRING, ith_parameter, name, count));

  if (ith_parameter == 0) {
    m_sddmid = key(m_keywords.c_str()) + ":$" + key(name);
  }
  
  else if (ith_parameter == 1 && getParamSpec(0).getType() == ENUMERATED && (getParamSpec(0).getName() == "Delimiter" || getParamSpec(0).getName() == "Assign")) {
    m_sddmid = key(m_keywords.c_str()) + ":$" + key(name);
  }
  
  return true ;
}

bool CommandSpec::set_parameter_expression(
  const UInt		ith_parameter ,
  const char *		name,
  const Count &	        count )
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::set_parameter_expression( const UInt ith_parameter , const char * name, const Count & count )"); /* %TRACE% */
  m_paramSpecList.addParamSpec(ParamSpec(EXPRESSION, ith_parameter, name, count));

  if (ith_parameter == 0) {
    m_sddmid = key(m_keywords.c_str()) + ":$" + key(name);
  }
  
  else if (ith_parameter == 1 && getParamSpec(0).getType() == ENUMERATED && (getParamSpec(0).getName() == "Delimiter" || getParamSpec(0).getName() == "Assign")) {
    m_sddmid = key(m_keywords.c_str()) + ":$" + key(name);
  }
  
  return true ;
}

bool CommandSpec::set_add_tag_list(const String &tag) {
  (*--(--m_paramSpecList.end())).setAddTagList(tag);
  return true;
}
  
bool CommandSpec::set_in_tag_list(const String &tag) {
  (*--(--m_paramSpecList.end())).setInTagList(tag);
  return true;
}
  
bool CommandSpec::add_nested_command( const Identifier & id )
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::add_nested_command( const Identifier & id )"); /* %TRACE% */
  ThrowAssert( m_blockCommand /* Only command blocks may have nested commands */ );

  const IdentifierVector::iterator i = std::lower_bound(m_nested.begin(), m_nested.end(), id);

  const bool result = (i == m_nested.end()) || (*i != id);

  if (result)
    m_nested.insert(i, id);

  return result ;
}


void
CommandSpec::fixup_parameters()
{
  std::string sddmid = m_sddmid;

  if (m_sddmkey)
    return;
  
  if (sddmid == "!")
    return;

  Lexer::LexemVector lex = Lexer::tokenize(sddmid);
    
  enum {START, SEPARATOR, NODE, VALUE};
  int state = START;
  stk::sddm::Key *start_key = 0;
  stk::sddm::Key *k = 0;
  bool make_root = false;
  
  for (Lexer::LexemVector::const_iterator it = lex.begin(); it != lex.end(); ) {
    switch (state) {
    case START:
      switch ((*it).getToken()) {
      case Lexer::TOKEN_DIVIDE:
        make_root = true;
        ++it;
        break;

      default:
        break;
      }
      state = NODE;
      break;

    case SEPARATOR:
      switch ((*it).getToken()) {
      case Lexer::TOKEN_DIVIDE:
        state = NODE;
        ++it;  
        break;

      case Lexer::TOKEN_COLON:
        state = VALUE;
        ++it;
        break;

      case Lexer::TOKEN_SEMI:
        k = start_key;
        state = SEPARATOR;
        ++it;
        break;

      case Lexer::TOKEN_IDENTIFIER:
        state = NODE;
        break;
        
      case Lexer::TOKEN_END:
        ++it;
        break;

      default:
        throw RuntimeError() << "Invalid Lexer token " << (*it).getToken() << " while processing state SEPARATOR" << std::endl
                             << "While processing sddmid " << sddmid << std::endl
                             << ErrorTrace;
        break;
      }
      break;

    case NODE:
      switch ((*it).getToken()) {
      case Lexer::TOKEN_IDENTIFIER: {
        stk::sddm::Key *new_key = 0;
        if (start_key == 0) {
          new_key = Parser::db_taxonomy().createKey(identifier().to_string(), (*it).getString().c_str());

          start_key = m_sddmkey = new_key;
          if (make_root)
            Parser::db_taxonomy().makeRoot(new_key);
        }
        else {
          if ((*(it + 1)).getToken() == Lexer::TOKEN_DIVIDE)
            new_key = Parser::db_taxonomy().createKey(identifier().to_string(), (*it).getString().c_str());
          else
            new_key = Parser::db_taxonomy().createTerminalKey(identifier().to_string(), (*it).getString().c_str());

          k->addChild(new_key);
        }
        
        k = new_key;
        ++it;
      }
        break;

      default:
        throw RuntimeError() << "Invalid Lexer token " << (*it).getToken() << " " << (*it).getString() << " while processing state NODE" << std::endl
                             << "While processing sddmid " << sddmid << std::endl
                             << ErrorTrace;
        break;
      }
      state = SEPARATOR;
      break;
        
    case VALUE:
      switch ((*it).getToken()) {
      case Lexer::TOKEN_IDENTIFIER:
        ++it;
        break;
          
      case Lexer::TOKEN_DOLLAR:
        ++it;
        k->setType(type_of(*this, it));
        ++it;
        break;

      default:
        throw RuntimeError() << "Invalid Lexer token " << (*it).getToken() << " while processing state VALUE" << std::endl
                             << "While processing sddmid " << sddmid << std::endl
                             << ErrorTrace;
        break;
      }
      state = SEPARATOR;
    }
  }
}


void
CommandSpec::fixup_nesting()
{
  CommandSpecMap &command_spec_map = Parser::db_command_spec();

//  prsrout << "For " << m_keywords << " adding" << dendl;

  for (IdentifierVector::const_iterator it2 = m_nested.begin(); it2 != m_nested.end(); ++it2) {
    const CommandSpec *child_command_spec = command_spec_map[*it2];

//    prsrout << "  " << child_command_spec->m_keywords << dendl;
    if (child_command_spec) {      
      stk::sddm::Key *child = child_command_spec->m_sddmkey;
      if (m_sddmkey && child) {
//        prsrout << "    " << m_sddmkey->getName() << " to " << child->getName() << dendl;
        m_sddmkey->addChild(child);
      }
    } 
  }
}


bool CommandSpec::is_valid() const
{
  /* %TRACE[ON]% */ Trace trace__("sierra::Prsr::CommandSpec::is_valid()"); /* %TRACE% */

  if (m_paramSpecList.size() == 0)
    return false;

  if (m_keywords.empty())
    return false;

  return true;
}


std::ostream &
operator<<(std::ostream &os, CommandSpec::ParamType type)
{
  /* %TRACE[SPEC]% */ Tracespec trace__("sierra::Prsr::operator<<(std::ostream &os, CommandSpec::ParamType type)"); /* %TRACE% */
  switch (type) {

  case CommandSpec::UNDEFINED:
    os << "UNDEFINED";
    break;
  case CommandSpec::LINE_PARAMETER:
    os << "LINE_PARAMETER ";
    break;
  case CommandSpec::TOKEN:
    os << "TOKEN";
    break;
  case CommandSpec::ENUMERATED:
    os << "ENUMERATED";
    break;
  case CommandSpec::INTEGER:
    os << "INTEGER";
    break;
  case CommandSpec::REAL:
    os << "REAL";
    break;
  case CommandSpec::STRING:
    os << "STRING";
    break;
  case CommandSpec::USTRING:
    os << "USTRING";
    break;
  case CommandSpec::QUOTED_STRING:
    os << "QUOTED_STRING";
    break;
  case CommandSpec::EXPRESSION:
    os << "EXPRESSION";
    break;
  default:
    os << "<undefined>";
  }
  return os;
}


Diag::Writer &
CommandSpec::verbose_print(
  Diag::Writer &		dout) const
{
  /* %TRACE[SPEC]% */ Tracespec trace__("sierra::Prsr::CommandSpec::verbose_print( Diag::Writer & dout) const"); /* %TRACE% */
  if (dout.shouldPrint()) {
    dout << "CommandSpec" << push << dendl;

    dout.m(LOG_MEMBERS) << "m_blockCommand, " << m_blockCommand << dendl;
    dout.m(LOG_MEMBERS) << "m_requiredParameters, " << m_requiredParameters << dendl;
    dout.m(LOG_MEMBERS) << "m_optionalParameters, " << m_optionalParameters << dendl;
    dout.m(LOG_MEMBERS) << "m_summary, " << m_summary << dendl;
    dout.m(LOG_MEMBERS) << "m_keywords, " << m_keywords << dendl;

    dout << pop;
  }

  return dout;
}

Diag::Writer &
CommandSpec::ParamSpec::verbose_print(
  Diag::Writer &		dout) const
{
  /* %TRACE[SPEC]% */ Tracespec trace__("sierra::Prsr::CommandSpec::ParamSpec::verbose_print( Diag::Writer & dout) const"); /* %TRACE% */
  if (dout.shouldPrint()) {
    dout << "CommandSpec::ParamSpec" << push << dendl;

    dout.m(LOG_MEMBERS) << "m_type, " << m_type << dendl;
    dout.m(LOG_MEMBERS) << "m_id, " << m_id << dendl;
    dout.m(LOG_MEMBERS) << "m_name, " << m_name << dendl;
    dout.m(LOG_MEMBERS) << "m_token, " << m_token << dendl;
    dout.m(LOG_MEMBERS) << "m_count, " << m_count << dendl;
    dout.m(LOG_MEMBERS) << "m_enumeration, " << Diag::c_ptr(m_enumeration) << dendl;

    dout << pop;
  }

  return dout;
}

} // namespace Prsr
} // namespace sierra
