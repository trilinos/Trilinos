/*--------------------------------------------------------------------*/
/*    Copyright 2000 - 2008, 2010 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
/*        Copyright 2000 Sandia Corporation, Albuquerque, NM.         */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_CommandSpec_h
#define SIERRA_CommandSpec_h

#include <list>
#include <vector>
#include <utility>
#include <limits>

#include <parser/Prsr_Identifier.h>
#include <parser/Prsr_Enumeration.h>

#include <stk_sddm/Tag.hpp>
#include <stk_sddm/DataTypes.hpp>

namespace stk {
namespace sddm {

class Lexeme;

class Occurs
{
public:
  Occurs(size_t count)
    : m_minOccurs(count),
      m_maxOccurs(count)
  {}

  Occurs(size_t min_occurs, size_t max_occurs)
    : m_minOccurs(min_occurs),
      m_maxOccurs(max_occurs)
  {}

  std::ostream &printSpecification(std::ostream &os) const;

private:
  size_t                m_minOccurs;
  
};


template <class BinOp>
class Validate
{
  Validate(const typename BinOp::argument_type &t)
    : m_op(std::binder2nd<BinOp>(BinOp(), t))
  {}

  typename BinOp::result_type test(const typename BinOp::argument_type &t) {
    return m_op(t);
  }

private:
  BinOp         m_op;
};


class CommandSpec
{
public:
  typedef std::vector<Identifier> IdentifierVector;

  ~CommandSpec();

  //: Construct a command specification with its unique identifier
  explicit CommandSpec( const Identifier & , bool is_block = false );

  //: Add nested command identifier to the list
  //  The method will fail (return false) is the command specification
  //  is not a command-block.
  bool add_nested_command( const Identifier & );

  void fixup_parameters();
  
  void fixup_nesting();
  
  /*--------------------------------------------------------------------*/
  /* Complete the command specification */

  bool set_sddmid(const char * const sddmid);
  
  //: Set the token sequence, all contiguous white spaces in the
  //: keywords will be replaced with a single space character.
  bool set_keywords( const char * const token_seq );

  //: SDDM taxonomy
  bool set_sddm_key(const char * const keywords);

  bool sddm_cleanup();
  
  //: Set the summary
  bool set_summary( const char * const summary );

  //: Set the description.
  bool set_description( const char * const description);

  //: Set the command specification to have a line value.
  bool set_line_parameter(const char * name);

  //: Set the identifier to have a line value.
  bool set_identifier();

  //: Set the number of parameters with specified value-types.
  bool set_number_of_parameters( const UInt required , const UInt optional );

  //: Set the ith parameter to be a required token.
  bool set_parameter_token( const UInt ith_parameter ,
			    const char * const token );

  bool set_parameter_optional(const UInt ith_parameter, bool optional);

  //: Set the ith parameter to be an list of real values.
  bool set_parameter_real( const UInt ith_parameter, const char *name,
			   const Count &count = Count(1, 1), const Range<Real> &range = Range<Real>(0, 0.0, 0.0));

  //: Set the ith parameter to be an list of integer values.
  bool set_parameter_int( const UInt ith_parameter, const char * name,
			  const Count &count = Count(1, 1), const Range<Int> &range = Range<Int>(0, 0, 0));

  //: Set the ith parameter to be an list of enumerated values.
  bool set_parameter_enum( const UInt ith_parameter, const char * name,
			   const Enumeration * enumerated, const Count & count = Count(1,1));

  //: Set the ith parameter to be an list of strings.
  bool set_parameter_string( const UInt ith_parameter ,
			     const char * name,
			     const Count & count = Count(1,1) );

  //: For upper-case strings
  bool set_parameter_ustring( const UInt ith_parameter ,
			      const char * name,
			      const Count & count = Count(1,1) );

  //: For quoted strings
  bool set_parameter_quoted_string( const UInt ith_parameter ,
				    const char * name,
				    const Count & count = Count(1,1) );

  bool set_parameter_expression(const UInt ith_parameter ,
				const char *		name,
				const Count &	count = Count(1,1) );

  bool set_add_tag_list(const String &tag);
  
  bool set_in_tag_list(const String &tag);
  
  /*--------------------------------------------------------------------*/
  /* QUERIES */

  //: Unique identifier of the command.
  const Identifier & identifier() const {
    return m_identifier;
  }

  //: Query if a block
  bool is_block() const {
    return m_blockCommand;
  }

  //: Query the nested command identifiers for a block
  const IdentifierVector &nested_identifiers() const  {
    return m_nested ;
  }

  //: sddmid
  const char *get_sddmid() const {
    return m_sddmid.c_str();
  }

  //: sddmkey
  stk::sddm::Key *get_sddmkey() const {
    return m_sddmkey;
  }

  //: Key word sequence of the command
  const char * keywords() const {
    return m_keywords.c_str();
  }

  //: Query summary
  const std::string &summary() const {
    return m_summary ;
  }

  //: Query description
  const std::string &description() const {
    return m_description ;
  }

  //: Query if the line following the token sequence is desired.
  bool uc_identifier() const {
    if (m_requiredParameters + m_optionalParameters == 0)
      return false;
    else
      return getParamSpec(0).getType() == IDENTIFIER;
  }


  //: Query total number of parameters
  size_t number_of_parameters() const {
    return m_optionalParameters + m_requiredParameters;
  }

  //: Query number of required parameters
  size_t number_required_parameters() const {
    return m_requiredParameters;
  }

  //: Query number of optional parameters
  size_t number_optional_parameters() const {
    return m_optionalParameters;
  }

  //: Types of parameters
  enum ParamType {
    UNDEFINED ,
    ENDLINE,		// Specification endline terminal
    LINE_PARAMETER ,	// Entire line as a single string
    IDENTIFIER ,	// Entire line as an identifier
    TOKEN ,		// One particular token is required
    ENUMERATED ,	// Any identifier within a list is required
    INTEGER ,		// Integer
    REAL ,		// Floating point
    STRING ,		// White space terminated character sequence
    USTRING ,		// White space terminated upper-case character sequence
    QUOTED_STRING,	// White space terminated quoted character sequence
    EXPRESSION		// Parenthesis bounded character sequence
  };

  typedef ParamType TYPE;

  //: Returns true if the command spec is defined correctly
  //: Returns false if there appear to be problems with the
  //: command spec.
  bool is_valid() const;

  /*--------------------------------------------------------------------*/

  Diag::Writer &verbose_print(Diag::Writer &dout) const;

  class ParamSpec
  {
  public:
    ParamSpec(ParamType type = UNDEFINED, int id = 0, const String &name = "") ;

    ParamSpec(ParamType type, int id, const String &name,
	      const Count & count, const Range<Int> &range);

    ParamSpec(ParamType type, int id, const String &name,
	      const Count & count, const Range<Real> &range);

    ParamSpec(ParamType type, int id, const String &name,
	      const Count & count, const Enumeration *enumerated);

    ParamSpec(ParamType type, int id, const String &name,
	      const Count & count);

    ParamSpec(const ParamSpec &param_spec)
      : m_type(param_spec.m_type),
        m_optional(param_spec.m_optional),
        m_id(param_spec.m_id),
        m_name(param_spec.m_name),
        m_token(param_spec.m_token),
        m_count(param_spec.m_count),
        m_enumeration(param_spec.m_enumeration),
        m_intRange(param_spec.m_intRange),
        m_realRange(param_spec.m_realRange),
        m_addTagList(param_spec.m_addTagList),
        m_inTagList(param_spec.m_inTagList)
    {}

    ParamType getType() const {
      return m_type;
    }

    ParamSpec &setType(ParamType type) {
      m_type = type;
      return *this;
    }

    bool getOptional() const {
      return m_optional;
    }

    const String &getAddTagList() const {
      return m_addTagList;
    }
    
    const String &getInTagList() const {
      return m_inTagList;
    }
    
    ParamSpec &setOptional(bool optional) {
      m_optional = optional;

      return *this;
    }

    int getId() const {
      return m_id;
    }

    ParamSpec &setId(int id) {
      m_id = id;
      return *this;
    }

    const String &getName() const {
      return m_name;
    }

    ParamSpec &setName(const char *name) {
      m_name = name;
      return *this;
    }

    const String &getToken() const {
      return m_token;
    }

    ParamSpec &setToken(const String &token) {
      m_token = token;
      return *this;
    }

    const Count &getCount() const {
      return m_count;
    }

    ParamSpec &setCount(const Count &count) {
      m_count = count;
      return *this;
    }

    const Enumeration &getEnumeration() const {
      return *m_enumeration;
    }

    ParamSpec &setEnumeration(const Enumeration *enumerated) {
      m_enumeration = enumerated;
      return *this;
    }

    const Range<Int> &getIntRange() const {
      return m_intRange;
    }

    ParamSpec &setIntRange(const Range<Int> &int_range) {
      m_intRange = int_range;
      return *this;
    }

    const Range<Real> &getRealRange() const {
      return m_realRange;
    }

    ParamSpec &setRealRange(const Range<Real> &real_range) {
      m_realRange = real_range;
      return *this;
    }

    void setAddTagList(const String &tag) {
      m_addTagList = tag;
    }
    
    void setInTagList(const String &tag) {
      m_inTagList = tag;
    }
    
//     bool matches(const Lexeme &lexeme) const;

//     bool terminates(const Lexeme &lexeme) const;

    bool checkCount(int count) const;

    bool inRange(Real x) const {
      return m_realRange.inRange(x);
    }

    bool inRange(Int i) const {
      return m_intRange.inRange(i);
    }

    std::ostream &printSpecification(std::ostream &os, bool print_triplets) const;

    Diag::Writer &verbose_print(Diag::Writer &dout) const;

  private:
    ParamType			m_type;
    bool			m_optional;
    int				m_id;
    String			m_name;
    String			m_token;
    Count			m_count;
    const Enumeration *		m_enumeration;
    Range<Int>			m_intRange;
    Range<Real>			m_realRange;
    String                      m_addTagList;
    String                      m_inTagList;
  };

  class ParamSpecList : public std::list<ParamSpec>
  {
  public:
    ParamSpecList()
      : std::list<ParamSpec>()
    {
      push_back(ParamSpec(ENDLINE, 100, "<endline>"));
    }

    void addParamSpec(const ParamSpec &param_spec) {
      insert(--end(), param_spec);
    }

//     void setLastOptional(bool optional) {
//       (*--(--end())).setOptional(optional);
//     }
  };


  const ParamSpecList &getParamSpecList() const {
    return m_paramSpecList;
  }

  const ParamSpec &getParamSpec(int id) const;

  std::ostream &printSpecification(std::ostream &os, bool print_triplets) const;
  std::ostream &printTaxonomy(std::ostream &os) const;

  inline friend Diag::Writer &operator<<(Diag::Writer &dout, const ParamSpec &param_spec);

public:
  static std::vector<Identifier> s_activeIdentifiers;
  
private:
  Identifier                    m_identifier;
  const bool                    m_blockCommand;         ///< Is / is not a block
  Int                           m_requiredParameters;   ///< Number of parameters required
  Int                           m_optionalParameters;   ///< Number of parameters optional

  std::string                   m_summary;              ///< Command summary
  std::string                   m_description;          ///< Command description
  std::string                   m_keywords;             ///< Keyword token sequence (will be obsolete)
  stk::sddm::Key *              m_sddmkey;              ///< SDDM key corresponding to keywords
  std::string                   m_sddmid;
  ParamSpecList                 m_paramSpecList;        ///< Parameter specifications
  IdentifierVector              m_nested;               ///< Children command block and line identifiers
};

typedef std::map<Identifier, CommandSpec *> CommandSpecMap;

/*----------------------------------------------------------------------*/

struct printCommandSpecification_
{
  printCommandSpecification_(const CommandSpec &command_spec, bool print_triplets)
    : m_commandSpec(command_spec),
      m_printTriplets(print_triplets)
  {}

  const CommandSpec &		m_commandSpec;
  bool                          m_printTriplets;
};

struct printParamSpecification_
{
  printParamSpecification_(const CommandSpec::ParamSpec &param_spec, bool print_triplets)
    : m_paramSpec(param_spec),
      m_printTriplets(print_triplets)
  {}

  const CommandSpec::ParamSpec &	m_paramSpec;
  bool                          m_printTriplets;
};

struct printCountSpecification_
{
  printCountSpecification_(const CommandSpec::Count &count)
    : m_count(count)
  {}

  const CommandSpec::Count &		m_count;
};

inline printCommandSpecification_ specification(const CommandSpec &command_spec, bool print_triplets) {
  return printCommandSpecification_(command_spec, print_triplets);
}

inline printParamSpecification_ specification(const CommandSpec::ParamSpec &param_spec, bool print_triplets) {
  return printParamSpecification_(param_spec, print_triplets);
}

inline printCountSpecification_ specification(const CommandSpec::Count &count) {
  return printCountSpecification_(count);
}

inline std::ostream &operator<<(std::ostream &os, const printCommandSpecification_ &command_spec) {
  return command_spec.m_commandSpec.printSpecification(os, command_spec.m_printTriplets);
}

inline std::ostream &operator<<(std::ostream &os, const printParamSpecification_ &param_spec) {
  return param_spec.m_paramSpec.printSpecification(os, param_spec.m_printTriplets);
}

inline std::ostream &operator<<(std::ostream &os, const printCountSpecification_ &count_spec) {
  return count_spec.m_count.printSpecification(os);
}

inline sierra::RuntimeDoomedSymmetric &operator<<(sierra::RuntimeDoomedSymmetric &es, const printCommandSpecification_ &command_spec) {
  std::ostringstream strout;
  command_spec.m_commandSpec.printSpecification(strout, false);
  return es << strout.str();
}

inline sierra::RuntimeDoomedSymmetric &operator<<(sierra::RuntimeDoomedSymmetric &es, const printParamSpecification_ &param_spec) {
  std::ostringstream strout;
  param_spec.m_paramSpec.printSpecification(strout, false);
  return es << strout.str();
}

inline sierra::RuntimeDoomedSymmetric &operator<<(sierra::RuntimeDoomedSymmetric &es, const printCountSpecification_ &count_spec) {
  std::ostringstream strout;
  count_spec.m_count.printSpecification(strout);
  return es << strout.str();
}

std::ostream &operator<<(std::ostream &os, CommandSpec::ParamType type);

inline Diag::Writer &operator<<(Diag::Writer &dout, const CommandSpec &command_spec) {
  return command_spec.verbose_print(dout);
}

inline Diag::Writer &operator<<(Diag::Writer &dout, const CommandSpec::ParamSpec &param_spec) {
  return param_spec.verbose_print(dout);
}

} // namespace Prsr
} // namespace sierra

typedef sierra::Prsr::CommandSpec Prsr_CommandSpec;

#endif
