#ifndef STK_ExprEvalFAD_hpp
#define STK_ExprEvalFAD_hpp

#include <string>
#include <limits>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <cctype>

#include <stk_util/util/string_case_compare.hpp>

#include <Sacado_Fad_DFad.hpp>

namespace stk_classic {
namespace expreval {
namespace fad {

typedef Sacado::Fad::DFad<double> FADDouble;

/**
 * @brief Class <b>CFunctionBase</b> is a base class for calling a function from the
 * expression evaluator.  Classes which inherit from this function must implement a
 * <b>operator()</b> function which accepts an argument count and an array of length
 * that count of pointer to the argument values.  This function should perform whatever
 * argument manipulation is required to call the actual implementation of the function.
 *
 */
class CFunctionBase
{
public:
  /**
   * Creates a new <b>CFunctionBase</b> instance.
   *
   * @param arg_count		an <b>int</b> value of the number of arguments
   *				expected when <b>operator()</b> is called.
   */
  explicit CFunctionBase(int arg_count)
    : m_argCount(arg_count)
  {}

  virtual ~CFunctionBase()
  {}

  /**
   * @brief Member function <b>operator()</b> is the pure virtual interface to the
   * function.  This function is overloaded to implement the calculation to be performed
   * when the function is called.
   *
   * @param argc		an <b>int</b> value of the number of arguments.
   *
   * @param argv		a <b>FADDouble</b> pointer to an array of FADDouble
   *				pointers to the actual function arguments.
   *
   * @return			a <b>FADDouble</b> value of the result of the execution
   *				of the function.
   */
  virtual FADDouble operator()(int argc, const FADDouble *argv) = 0;

  /**
   * @brief Member function <b>getArgCount</b> returns the argument count for this
   * function.
   *
   * @return			an <b>int</b> value of the argument count for this
   *				function.
   */
  int getArgCount() const {
    return m_argCount;
  }

private:
  int		m_argCount;			///< Argument count for this function
};


/**
 * @brief Class <b>CFunction</b> declares a base template for template instantiation
 * of C-like called functions.
 *
 */
template <class S>
class CFunction;

class Node;

/**
 * Class <b>Eval</b> parses and evaluates mathematical expressions.
 *
 */
class Eval
{
public:
  static const double s_e;			///< Constant "e"
  static const double s_pi;			///< Constant "pi"
  static const double s_false;			///< Constant "false"
  static const double s_true;			///< Constant "true"

public:
  /**
   * @brief Class <b>CFunctionMap</b> maps function names to c-style function pointers via
   * the <b>CFunctionBase</b> base class.  The mapping is case insensitive.
   *
   */
  class CFunctionMap : public std::map<std::string, CFunctionBase *, LessCase>
  {
  public:
    /**
     * Creates a new <b>CFunctionMap</b> instance.
     *
     * The <b>CFunctionMap</b> is populated with common C functions on construction.
     *
     */
    CFunctionMap();

    /**
     * Destroys a <b>CFunctionMap</b> instance.
     *
     */
    ~CFunctionMap();
  };

  /**
   * @brief Typedef <b>ConstantMap</b> maps a constant name to a double constant.
   * The mapping is case insensitive.
   *
   */
  typedef std::map<std::string, double, LessCase> ConstantMap;

  typedef std::set<std::string> UndefinedFunctionSet;

  /**
   * @brief Class <b>Variable</b> defines a variable for the expression evaluation
   * virtual machine.  Typing is limited to int and double.  While all operations are
   * performed as doubles, the value and results may be retrieved and stored as integer
   * values.
   *
   * Array indexing is also supported, but no bounds checkin is available and only double
   * type array may be addressed.
   *
   * Variable are implemented as a pointer, "typed" by the <b>Type</b> enumeration
   * to the actual value.  At construction, this may point to a "local" double value or
   * integer value, or the address may be specified.  In either case, a type and address
   * for the variable may be provided via variable "binding".
   *
   */
  class Variable
  {
  public:
    /**
     * @brief Enumeration <b>Type</b> lists the variable data types.  <b>double</b> and
     * <b>int</b> are currently supported.
     *
     */
    enum Type {FADDOUBLE};

    /**
     * Creates a new <b>Variable</b> instance.
     *
     */
    Variable()
      : m_type(FADDOUBLE),
	m_FADDoublePtr(&m_FADDoubleValue),
	m_FADDoubleValue(0.0)
    {}

    /**
     * Creates a new <b>Variable</b> instance.  The new variable will be local and
     * have the type specified by <b>type</b>.
     *
     * @param type		an <b>Type</b> value of the type of the new
     *				variable.
     *
     */
    explicit Variable(Type type)
      : m_type(type)
    {
      switch (type) {
	case FADDOUBLE:
	  m_FADDoublePtr = &m_FADDoubleValue;
	  m_FADDoubleValue = 0.0;
	  break;
      }
    }

    /**
     * Creates a new <b>Variable</b> instance.  The new variable will use the
     * address specified by <b>address</b> and be of type FADDouble.  This address must
     * remain in scope during the lifetime of the variable are until rebound to a new
     * address.
     *
     * @param address		an <b>FADDouble</b> reference to the value for this
     *				variable.
     *
     */
    explicit Variable(FADDouble &address)
      : m_type(FADDOUBLE),
	m_FADDoublePtr(&address),
	m_FADDoubleValue(0.0)
    {}

    /**
     * @brief Member function <b>operator=</b> assigns a new value to the variable.
     * If the variable type is integer, the value is converted to integer before
     * assignment.
     *
     * @param value		a <b>FADDouble</b> value that is to be assigned to the
     *				variable.
     *
     * @return			a <b>Variable</b> reference to the variable.
     */
    Variable &operator=(const FADDouble &value) {
      if (m_type == FADDOUBLE)
	*m_FADDoublePtr = value;
      return *this;
    }

    /**
     * @brief Member function <b>operator=</b> assigns a new value to the variable.
     * If the variable type is FADDouble, the value is converted to FADDouble before assignment.
     *
     * @param value		an <b>int</b> value that is to be assigned to the
     *				variable.
     *
     * @return			a <b>Variable</b> reference to the variable.
     */
    Variable &operator=(const int &value) {
      if (m_type == FADDOUBLE)
	*m_FADDoublePtr = (FADDouble) value;
      return *this;
    }

    /**
     * @brief Member function <b>operator[]</b> returns a value from an array of
     * FADDouble values.  No bounds checkin is performed.  Not even if the variable is and
     * array is checked.
     *
     * @param index		a <b>FADDouble</b> value of the zero based index into
     *				the array to retrieve the value.
     *
     * @return			a <b>FADDouble</b> reference to the value.
     */
    inline FADDouble &operator[](FADDouble index) {
      if (m_type != FADDOUBLE)
	throw std::runtime_error("Only FADDouble arrays allowed");

      if ((void *) m_FADDoublePtr == 0)
	throw std::runtime_error("Unbound variable");

      int i = 0; // (int) index;

      return m_FADDoublePtr[i];
    }

    /**
     * @brief Member function <b>operator[]</b> returns a value from an array of
     * FADDouble values.  No bounds checkin is performed.  Not even if the variable is and
     * array is checked.
     *
     * @param index		a <b>int</b> value of the zero based index into the
     *				array to retrieve the value.
     *
     * @return			a <b>FADDouble</b> reference to the value.
     */
    inline FADDouble &operator[](int index) {
      if (m_type != FADDOUBLE)
	throw std::runtime_error("Only FADDouble arrays allowed");

      if ((void *) m_FADDoublePtr == 0)
	throw std::runtime_error("Unbound variable");

      return m_FADDoublePtr[index];
    }

    /**
     * @brief Member function <b>bind</b> binds variable to the address of the
     * specified value.  The type is converted to FADDouble.  This address must remain in
     * scope during the lifetime of the variable are until rebound to a new address.
     *
     * @param value_ref	a <b>FADDouble</b> reference to be used for this variable.
     *
     * @return			a <b>Variable</b> reference to the variable.
     */
    inline Variable &bind(FADDouble &value_ref) {
      m_type = FADDOUBLE;
      m_FADDoublePtr = &value_ref;
      return *this;
    }

    /**
     * @brief Member function <b>unbind</b> binds the variable to the local value
     * and reinitializes that value to zero.  The type is left unchanged.
     *
     * @return			a <b>Variable</b> reference to the variable.
     */
    inline Variable &unbind() {
      switch (m_type) {
	case FADDOUBLE:
	  m_FADDoublePtr = &m_FADDoubleValue;
	  m_FADDoubleValue = 0.0;
	  break;
      }
      return *this;
    }

    /**
     * @brief Member function <b>getValue</b> returns the variable value as a FADDouble.
     *
     * @return			a <b>FADDouble</b> value of the value of the variable.
     */
    inline FADDouble getValue() const {
      switch (m_type) {
	case FADDOUBLE:
	  return *m_FADDoublePtr;
      }
      throw std::runtime_error("Invalid variable type");
    }

  private:
    Type	m_type;				///< Variable data type
    FADDouble *	m_FADDoublePtr;		        ///< Pointer to value as FADDouble
    FADDouble	m_FADDoubleValue;		///< Local variable value as FADDouble
  };

  /**
   * @brief Class <b>VariableMap</b> implements a mapping from name to a pointer to
   * a <b>Variable</b> object.  The mapping is case insensitive.
   *
   */
  class VariableMap : public std::map<std::string, Variable *, LessCase>
  {
  public:
    /**
     * @brief Typedef <b>value_type</b> is the value_type of the
     * <b>std::map</b> subclass.  The mapping is case insensitive.
     *
     */
    typedef std::map<std::string, Variable *, LessCase>::value_type value_type;

    /**
     * @brief Class <b>Resolver</b> is a base class for a variable name to value
     * resolver.  After parsing of an expression, a list of variable names is collected.
     * The <b>resolve</b> member function will be called for each of these during
     * the expression resolve stage.  This function may bind a new value to the variable,
     * ignore the variable by leaving the initialized local binding, or throw an exception
     * stating the variable cannot be resolved.
     *
     */
    class Resolver
    {
    public:
      /**
       * Creates a new <b>Resolver</b> instance.
       *
       */
      Resolver()
      {}

      /**
       * Destroys a <b>Resolver</b> instance.
       *
       */
      virtual ~Resolver()
      {}

      /**
       * @brief Member function <b>resolve</b> is the virtual member function for
       * the variable resolver.
       *
       * @param it		a <b>VariableMap::iterator</b> reference to the
       *			variable whose name is to be resolved.
       */
      virtual void resolve(VariableMap::iterator &it) = 0;
    };

  private:
    /**
     * @brief Member function <b>delete_variable</b> is a function which will delete
     * the variabel pointed to by <b>t</b>.
     *
     * @param t			a <b>value_type</b> reference to the variable to be
     *				deleted.
     */
    static void delete_variable(value_type &t) {
      delete t.second;
    }

    /**
     * @brief Class <b>DefaultResolver</b> implement a default resolver.
     *
     */
    class DefaultResolver : public Resolver
    {
    public:
      DefaultResolver()
      {}

      /**
       * Destroys a <b>DefaultResolver</b> instance.
       *
       */
      virtual ~DefaultResolver()
      {}

      /**
       * @brief Member function <b>resolve</b> implements the default resolvers
       * function, whihc does nothing.  I.E. lets the local variable values stand.
       *
       * @param it		a <b>VariableMap::iterator</b> variable ...
       */
      virtual void resolve(VariableMap::iterator &it)
      {}
    };

  public:
    /**
     * @brief Member function <b>getDefaultResolver</b> returns a reference to the
     * default resolver.
     *
     * @return a <b>Resolver</b> ...
     */
    static Resolver &getDefaultResolver();

    /**
     * Creates a new <b>VariableMap</b> instance with the specified variable name
     * resolver.
     *
     * @param resolver		a <b>Resolver</b> reference to the variable name
     *				resolver for this variable map.
     *
     */
    VariableMap(Resolver &resolver = getDefaultResolver())
      : std::map<std::string, Variable *, LessCase>(),
	m_resolver(resolver)
    {}

    /**
     * Destroys a <b>VariableMap</b> instance.  All variables are destroyed.
     *
     */
    ~VariableMap() {
      std::for_each(begin(), end(), &delete_variable);
    }

    /**
     * @brief Member function <b>operator[]</b> ...
     *
     * @param s			a <b>std::string</b> const reference to the variable's
     *				name.
     *
     * @return			a <b>Variable</b> pointer to the new variable.
     */
    Variable *operator[](const std::string &s) {
      std::pair<iterator,bool> i = insert(std::pair<const std::string, Variable *>(s, (Variable *) 0));
      if (i.second)
	(*i.first).second = new Variable();
      return (*i.first).second;
    }

    /**
     * @brief Member function <b>getResolver</b> returns a reference to the name
     * resolver.
     *
     * @return			a <b>Resolver</b> reference to the variable name
     *				resolver.
     */
    Resolver &getResolver() {
      return m_resolver;
    }

  private:
    Resolver &		m_resolver;		///< Reference to this variabl map name resolver
  };

  /**
   * @brief Member function <b>getCFunctionMap</b> returns a reference to the
   * C-style function map.
   *
   * @return			a <b>CFunctionMap</b> reference to the C-style
   *				function map.
   */
  static CFunctionMap &getCFunctionMap();

  /**
   * @brief Member function <b>getConstantMap</b> returns s reference to the defined
   * constants.
   *
   * @return			a <b>ConstantMap</b> reference to the defined
   *				constants.
   */
  static ConstantMap &getConstantMap();

  /**
   * @brief Member function <b>addFunction</b> adds a C-style function to the
   * C-style function map.
   *
   * @param name		a <b>std::string</b> const reference to the function
   *				name.
   *
   * @param function		a <b>CFunctionBase</b> pointer to the C-style
   *				function.
   */
  static void addFunction(const std::string &name, CFunctionBase *function) {
    getCFunctionMap()[name] = function;
  }

  /**
   * Creates a new <b>Eval</b> instance.
   *
   * @param resolver		a <b>Resolver</b> reference to the variable name
   *				resolver for this expression evaluator.
   *
   * @param expr		a <b>std::string</b> const reference to the
   *				expression to be parsed.
   */
  Eval(VariableMap::Resolver &resolver = VariableMap::getDefaultResolver(), const std::string &expr = "");

private:
  explicit Eval(const Eval &);
  Eval &operator=(const Eval &);

public:
  /**
   * Destroys a <b>Eval</b> instance.
   *
   */
  ~Eval();

  /**
   * @brief Member function <b>getExpression</b> returns the original text of the expression.
   *
   * @return			a <b>std::string</b> const reference to the
   *				expression.
   */
  inline const std::string &getExpression() const {
    return m_expression;
  }

  /**
   * @brief Member function <b>setExpression</b> gives the evaluator a new
   * expression to parse.
   *
   * @param expression		a <b>std::string</b> const refernce to the
   *				expression to be parsed.
   *
   * @return			an <b>Eval</b> reference to this expression
   *				evaluator.
   */
  inline Eval &setExpression(const std::string &expression) {
    m_expression = expression;
    m_parseStatus = false;
    return *this;
  }

  /**
   * @brief Member function <b>getVariableMap</b> returns a reference to the
   * variable map of this expression.
   *
   * @return			a <b>VariableMap</b> reference to the variable map.
   */
  inline VariableMap &getVariableMap() {
    return m_variableMap;
  }

  /**
   * @brief Member function <b>getVariableMap</b> returns a reference to the
   * variable map of this expression.
   *
   * @return			a <b>VariableMap</b> reference to the variable map.
   */
  inline UndefinedFunctionSet &getUndefinedFunctionSet() {
    return m_undefinedFunctionSet;
  }

  /**
   * @brief Member function <b>getSyntaxStatus</b> returns true if the expression has
   * been syntaxd successfully.
   *
   * @return			a <b>bool</b> value of true if the expression has
   *				been syntaxd successfully.
   */
  inline bool getSyntaxStatus() const {
    return m_syntaxStatus;
  }

  /**
   * @brief Member function <b>getParseStatus</b> returns true if the expression has
   * been parsed successfully.
   *
   * @return			a <b>bool</b> value of true if the expression has
   *				been parsed successfully.
   */
  inline bool getParseStatus() const {
    return m_parseStatus;
  }

  /**
   * @brief Member function <b>setValue</b> assigns a variables value in the
   * variable map.
   *
   * @param name		a <b>std::string</b> const reference of the name of the
   *				variable.
   *
   * @param value		a <b>FADDouble</b> value to be assigned to the
   *				variable.
   *
   * @return			an <b>Eval</b> reference to this expression
   *				evaluator.
   */
  inline Eval &setValue(const std::string &name, FADDouble value) {
    VariableMap::iterator it = m_variableMap.find(name);
    if (it != m_variableMap.end())
      *(*it).second = value;
    return *this;
  }

  /**
   * @brief Member function <b>getValue</b> returns the value of the variable
   * specified by <b>name</b>.
   *
   * @param name		a <b>std::string</b> const reference to the variable's
   *				name.
   *
   * @return			a <b>FADDouble</b> value of the variable.
   */
  inline FADDouble getValue(const std::string &name) {
    VariableMap::iterator it = m_variableMap.find(name);
    if (it == m_variableMap.end())
      throw std::runtime_error(std::string("Variable ") + name  + " not defined");
    return (*it).second->getValue();
  }

  /**
   * @brief Member function <b>newNode</b> allocates a new node.  The
   * new node is allocated on a node list so that it may be
   * deallocated properly on exception.
   *
   * @param opcode		a <b>int</b> value of the opcode for the node.
   *
   * @return                    a <b>Node</b> pointer to the newly allocated node.
   */
  Node *newNode(int op);
  
  /**
   * @brief Member function <b>bindVariable</b> binds the variable to the address of
   * the specified value.  This address must remain in scope during the lifetime of the
   * variable are until rebound to a new address.
   *
   * @param name		a <b>std::string</b> const reference to the variable's
   *				name.
   *
   * @param value_ref		a <b>FADDouble</b> reference to be used for this variable.
   *
   * @return			an <b>Eval</b> reference to this expression
   *				evaluator.
   */
  inline Eval &bindVariable(const std::string &name, FADDouble &value_ref) {
    VariableMap::iterator it = m_variableMap.find(name);
    if (it != m_variableMap.end())
      (*it).second->bind(value_ref);
    return *this;
  }

  /**
   * @brief Member function <b>getVariable</b> returns a reference to the variable
   * specified by <b>name</b>.
   *
   * @param name		a <b>std::string</b> const reference to the variable's
   *				name.
   *
   * @return			a <b>Variable</b> reference to the specified
   *				variable.
   */
  inline Variable &getVariable(const std::string &name) {
    VariableMap::iterator it = m_variableMap.find(name);
    if (it == m_variableMap.end())
      throw std::runtime_error(std::string("Variable ") + name  + " not defined");

    return *(*it).second;
  }

  /**
   * @brief Member function <b>parse</b> parses the expression.  If successful, the
   * parse status is set to true.
   *
   * @param expr		a <b>std::string</b> const reference to the
   *				expression to parse.
   *
   */
  inline void syntaxCheck(const std::string &expr) {
    setExpression(expr);
    syntax();
  }

  /**
   * @brief Member function <b>syntax</b> performs a syntax check on the current
   * expression.  If successful, the syntax status is set to true.
   *
   */
  void syntax();

  /**
   * @brief Member function <b>parse</b> parses the expression.  If successful, the
   * parse status is set to true.
   *
   * @param expr		a <b>std::string</b> const reference to the
   *				expression to parse.
   *
   */
  inline void parse(const std::string &expr) {
    setExpression(expr);
    parse();
  }

  /**
   * @brief Member function <b>parse</b> parses the current expression.  If
   * successful, the parse status is set to true.
   *
   */
  void parse();

  /**
   * @brief Member function <b>resolve</b> calls the variable name resolver for each
   * variable in the variable map.
   *
   */
  void resolve();

  /**
   * @brief Member function <b>evaluate</b> evaluates the expression.
   *
   * @return			a <b>FADDouble</b> value of the result on the
   *				expression evaluation.
   */
  FADDouble evaluate() const;

private:
  VariableMap		m_variableMap;		///< Variable map
  UndefinedFunctionSet	m_undefinedFunctionSet;	///< Vector of undefined functions

  std::string		m_expression;		///< Expression which was parsed.
  bool			m_syntaxStatus;		///< True if syntax is correct
  bool			m_parseStatus;		///< True if parsed successfully

  Node *		m_headNode;		///< Head of compiled expression
  std::vector<Node *>   m_nodes;                ///< Allocated nodes
};

} // namespace fad 
} // namespace expreval
} // namespace stk_classic

#endif // STK_ExprEvalFAD_hpp
