#ifndef STK_UTIL_DIAG_UserPlugin_h
#define STK_UTIL_DIAG_UserPlugin_h

#include <memory>
#include <map>
#include <vector>
#include <string>
#include <typeinfo>

#include <stk_util/diag/StringUtil.hpp>
#include <stk_util/util/Fortran.hpp>
#include <stk_util/diag/Writer_fwd.hpp>

/**
 * @file
 *
 * Historically, the term User Subroutine has referred to Fortran subroutines which are
 * generally called from within a procedure, allowing the user to create or select a
 * mathematical calculation to be applied to set of arguments.
 *
 * The purpose of this package is add User Plugins and User Functions as well as provide an
 * analogous implementation of the traditional User Subroutine.
 *
 * A User Plugin is a C++ class in which the application developer designs and implements
 * a base class that end users can create or select derivative classes to before the
 * desired operations.
 *
 * A User Subroutine is a function with a specific calling signature allowing any
 * registered subroutine with that signature to be called by its registered name.
 *
 * A User Function is a functional calculation which accepts one or more independent const
 * variables and returns a single dependent variable.  These variables may be of scalar,
 * vector or object quantities.  These functions can be created or selected by the end
 * user.
 *
 */

namespace sierra {
namespace Plugin {

std::string derived_id_name(int derived_id);

/**
 * Class <b>Registry</b> serves as a singleton for holding templatized
 * createInstance and UserSubroutine function pointers and pointer to class factory objects.
 * The registry is simply a mapping of name pairs to a void pointer.  The name pair
 * consists of the base class name and the derived class names.  And, since the only legal
 * way to get things into the registry is via the and UserSubroutine classes, there is no
 * signature checking performed at this level
 *
 * There should never be a need to instantiate this singleton as the UserSubroutine
 * registration process should perform that function when necessary.
 */

class Registry
{
public:

  /**
   * @brief Typedef <b>NamePair</b> is the derived class key.
   *
   */
  typedef std::pair<const std::type_info *, std::string> NamePair;

//   /**
//    * @brief Class <b>hash_nocase</b> implements a hash, case insensitive NamePair
//    * hash functor.
//    *
//    */
//   struct hash_nocase
//   {
//     size_t operator()(const NamePair &n) const {
//       return sierra::hash_string_nocase(n.second.c_str());
//     }
//   };

//   /**
//    * @brief Class <b>hash_nocase</b> implements a hash, case insensitive compare
//    * equal NamePair functor.
//    *
//    */
//   struct equal_nocase : public std::binary_function<NamePair, NamePair, bool>
//   {
//     bool operator()(const NamePair &lhs, const NamePair &rhs) const {
//       sierra::equal_nocase<NamePair::second_type> second_equal_nocase;

//       return *lhs.first == *rhs.first && second_equal_nocase(lhs.second, rhs.second);
//     }
//   };

  /**
   * @brief Class <b>less_nocase</b> implements a case insensitive NamePair compare
   * less functor.
   *
   */
  struct less_nocase : public std::binary_function<NamePair, NamePair, bool>
  {
    bool operator()(const NamePair &lhs, const NamePair &rhs) const {
      sierra::less_nocase<NamePair::second_type> second_less_nocase;

#ifdef SIERRA_TYPE_INFO_BEFORE_EQUALITY_BUG
      return (lhs.first->before(*rhs.first) && *lhs.first != *rhs.first)
	|| (*lhs.first == *rhs.first && second_less_nocase(lhs.second, rhs.second));
#else
      return lhs.first->before(*rhs.first)
	|| (*lhs.first == *rhs.first && second_less_nocase(lhs.second, rhs.second));
#endif
    }
  };

  /**
   * @brief Typedef <b>RegistryMap</b> is the registration map.
   *
   */
  typedef std::map<NamePair, void *, less_nocase> RegistryMap;

  static RegistryMap &getRegistryMap();

  /**
   * Creates a new <b>Registry</b> instance.
   *
   */
  Registry()
  {}

  /**
   * @brief Creates a new <b>Registry</b> instance and registers it, and more
   * importantly the derived class factory with the specified name pair.
   *
   * @param name_pair		a <b>NamePair</b> gives the name to the registry
   *				entry.
   */
  explicit Registry(const NamePair &name_pair) {
    registerIt(name_pair, this);
  }

  /**
   * Destructor <b>~Registry</b> is virtual to fake polymorphism so that the registry class can
   * utilize additional compiler/runtime checks that the registered class is indeed a class factory.
   *
   */
  virtual ~Registry()
  {}

  /**
   * @brief Member function <b>rootInstance</b> creates the singleton.
   *
   * @return			a <b>Registry</b> reference to the
   *				Registry singleton.
   */
  static Registry &rootInstance();

  /**
   * @brief Member function <b>registerDL</b> opens a dynamic library and optionally executes a "C"
   * registration function.
   *
   * If function is specified, and not zero length, the function must exist.  If the function name
   * specified with zero length, the function "dl_register" is executed if it exists.  If the
   * function name is not found, if platform specific fortran suffix is appended and the function is
   * searched again.
   *
   * If no function name is specified, no registration function is executed.
   *
   * NOTE: Loading C++ sharable objects results in static object construction upon load.
   *
   * @param so_path		a <b>char</b> const pointer to the path of the
   *				shareable object file.  If the path does not contains a
   *				'/' character, the file is searched for through the
   *				LD_LIBRARY_PATH envirnment variable.
   *
   * @param function_name	a <b>char</b> const pointer to the name of a registration function
   *				which will be called immediately after loading the sharable object.
   *
   */
  static void registerDL(const char *so_path, const char *function_name = 0);

  template <typename T>
  static T getsym(const char *sym);

  /**
   * @brief Member function <b>registerIt</b> registers a name pair with a void
   * pointer.
   *
   * If the name pair already exists within the registry, a std::invalid_argument exception
   * is thrown with the offending name pair called out.
   *
   * @param name_pair		a <b>NamePair</b> const reference to the name pair
   *				to be registered.
   *
   * @param func_ptr		a <b>void</b> pointer to the function to be
   *				registered.
   *
   * @throws			a <b>std::invalid_argument</b> exception is thrown
   *				if there is an instance creation function already
   *				registered for the derived class.
   *
   */
  void registerIt(const NamePair &name_pair, void *func_ptr);

  /**
   * @brief Member function <b>getPluginPtr</b> find the function with the name pair
   * specified.
   *
   * If the name pair does not exist within the registry, a std::invalid_argument
   * exception is thrown with the offending name pair called out.
   *
   * @param name_pair		a <b>NamePair</b> const reference to the name pair
   *				to be retrieved.
   *
   * @throws			a <b>std::invalid_argument</b> exception is thrown
   *				if the function is not found.
   *
   */
  void *getPluginPtr(const NamePair &name_pair) const;

  /**
   * @brief Member function <b>getFunctionPtr</b> find the function with the name pair
   * specified.
   *
   * If the name pair does not exist within the registry, a std::invalid_argument
   * exception is thrown with the offending name pair called out.
   *
   * @param name_pair		a <b>NamePair</b> const reference to the name pair
   *				to be retrieved.
   *
   * @throws			a <b>std::invalid_argument</b> exception is thrown
   *				if the function is not found.
   *
   */
  void *getFunctionPtr(const NamePair &name_pair) const;

  /**
   * @brief Member function <b>getFuncPtr</b> returns the function pointer with the
   * specfied <it>name_pair</i>.
   *
   * @param name_pair		a <b>NamePair</b> const reference to the registered
   *				name pair.
   *
   * @returns			a <b>void</b> function pointer with the specfied
   *				<it>name_pair</i>.
   */
  Registry *getFactoryPtr(const NamePair &name) const;

  /**
   * @brief Member function <b>getFuncPtr</b> returns the function pointer with the
   * specfied <it>name_pair</i>.
   *
   * @param name_pair		a <b>NamePair</b> const reference to the registered
   *				name pair.
   *
   * @returns			a <b>void</b> function pointer with the specfied
   *				<it>name_pair</i>.
   */
  void *getFuncPtr(const NamePair &name_pair) const;

  /**
   * @brief Member function <b>getDerivedNames</b> returns names assocaited with the
   * function pointers of the specified type.
   *
   * @param type		a <b>std::type_info</b> const reference to typeid to
   *				retrieve the derived names.
   *
   * @returns			a <b>std::vector<str::string></b> value of the derived
   *				names.
   */
  std::vector<std::string> getDerivedNames(const std::type_info &type) const;

  /**
   * Member template function <b>create</b> creates an instance of the desired
   * object by providing the factory responsible for generating that object type.  The
   * create factory is retrieved using the base class name specified by the factory base
   * class and the specified derived name.
   *
   * The derived factory is responsible for implementing the appropriate operator()
   * functions to constuct the derived object.
   *
   * @param derived_name	a <b>std::string</b> const reference to the derived
   *				object's name/
   *
   * @return			a <b>T</b> reference to the creation factory object.
   */
  template<class T>
  static T &create(const std::string &derived_name) {
    return static_cast<T &>(*Registry::rootInstance().getFactoryPtr(std::make_pair(&typeid(T), derived_name)));
  }

  /**
   * @brief Member function <b>dump</b> dumps the registry.
   *
   * @param os			a <b>std::ostream</b> reference to dump the registry
   *				to.
   *
   * @return			a <b>std::ostream</b> reference to <it>os</i>.
   */
  std::ostream &verbose_print(std::ostream &os) const;

  /**
   * @brief Member function <b>verbose_print</b> dumps the registry.
   *
   * @param os			a <b>std::ostream</b> reference to dump the registry
   *				to.
   *
   * @return			a <b>std::ostream</b> reference to <it>os</i>.
   */
  stk::diag::Writer &verbose_print(stk::diag::Writer &dout) const;
};

inline stk::diag::Writer &operator<<(stk::diag::Writer &dout, const Registry &registry) {
  return registry.verbose_print(dout);
}


/**
 * Template class <b>UserPlugin</b> is a template used for the association of base
 * and derived classed to be registered and created via the UserPlugin mechanism.  The
 * template traits enforces the signature matching of the base class constructor, derived class
 * constructor, derived class creator static function and the usage of the creator static
 * function.
 *
 * The registration requires a unique base class name for each base class type.  And,
 * since typeid is not reliable for that implementation, each base class is required to
 * implement a traits class with a Base typedef which specifies the base class, a
 * Signature typedef which specifies the signature of the create function and a static
 * function named getUserPluginCreatorName() which returns a const std::string reference to
 * the base classes name.  This name must be unique across users of the registry.  There
 * is no way to enforce this programatically, so it is recommended that a application
 * prefix be attached to the base class name.
 *
 */
template <class Creator, typename S = Creator *(*)()>
class UserPlugin
{
public:
  typedef S Signature;					///< Creator signature

private:
  UserPlugin();						///< Not implemented
  UserPlugin(const UserPlugin&);			///< Not implemented
  UserPlugin &operator=(const UserPlugin&);		///< Not implemented

public:
  /**
   * @brief Member function <b>instance</b> returns the instance of the registry,
   * cast as a UserPlugin registry..
   *
   * @return			a <b>UserPlugin</b> reference to the registry
   *				singleton.
   */
  static UserPlugin &instance() {
    return (UserPlugin<Creator, Signature> &) (Registry::rootInstance());
  }

  /**
   * @brief Member function <b>registerCreator</b> registers the base class name and
   * specified derived class name with the specified creator function.
   *
   * The base class name is determined by the BaseTraits template argument's
   * getUserPluginCreatorName() static member function.  The signature is defined by the
   * BaseTraits template argument's Signature typedef.
   *
   * @param derived_name	a <b>std::string</b> const reference to the derived
   *				class name.
   *
   * @param function		a <b>signature</b> function pointer to the creator
   *				function.
   *
   */
  static void registerCreator(const std::string &derived_name, Signature function) {
    Registry::rootInstance().registerIt(std::make_pair(&typeid(Signature), derived_name), (void *) function);
  }

  /**
   * @brief Member function <b>create</b> returns the createInstance() function
   * associated with the specified base class and derived_name.
   *
   * @param derived_name	a <b>std::string</b> const reference to the derived
   *				classes name.
   *
   * @throws			a <b>std::invalid_argument</b> exception is thrown
   *				if there is no instance creation function registered for
   *				the specified name of the base class.
   *
   * @return			a <b>Signature</b> function pointer to the instance
   *				create function.
   */
  static Signature create(const std::string &derived_name) {
    Signature creator_function = (Signature) Registry::rootInstance().getPluginPtr(std::make_pair(&typeid(Signature), derived_name));

    return (*creator_function);
  }

  /**
   * @brief Member function <b>create</b> returns the createInstance() function
   * associated with the specified base class and derived_name.
   *
   * @param derived_name	a <b>std::string</b> const reference to the derived
   *				classes name.
   *
   * @throws			a <b>std::invalid_argument</b> exception is thrown
   *				if there is no instance creation function registered for
   *				the specified name of the base class.
   *
   * @return			a <b>Signature</b> function pointer to the instance
   *				create function.
   */
  static Signature create(int derived_id) {
    return create(derived_id_name(derived_id));
  }

  /**
   * @brief Member function <b>exists</b> returns true if class of the type
   * specified by derived_name exists in BaseClass.
   *
   * @param derived_name	a <b>std::string</b> const reference to the derived
   *				classes name.
   *
   * @return			a <b>bool</b> of true if class of the type specified by derived_name exists in BaseClass.
   */
  static bool exists(const std::string &derived_name) {
    return Registry::rootInstance().getFuncPtr(std::make_pair(&typeid(Signature), derived_name)) != NULL;
  }

  static std::vector<std::string> getDerivedNames() {
    return Registry::rootInstance().getDerivedNames(typeid(Signature));
  }

  /**
   * @brief Class template <b>Register</b> registers the <i>createInstance()</i>
   * function with the <i>derived_name</i> on object creation.
   *
   * @param DerivedClass	a <b>class</b> which specifies the derived class
   *				which holds the <i>createInstance()</i> function.
   *
   */
  template <class DerivedClass>
  class Register
  {
  public:
    typedef DerivedClass XDerivedClass;

    /**
     * @brief Creates a new <b>Register</b> instance.  Upon creation, the
     * <i>DerivedClass::createInstance()</i> instance creation function is registered
     * with the <i>derived_name</i>.
     *
     * @param derived_name	a <b>std::string</b> const reference to the derived
     *				class' name.
     *
     */
    explicit Register(const std::string &derived_name)
      : m_function(DerivedClass::createInstance)
    {
      UserPlugin<Creator, Signature>::instance().registerCreator(derived_name, m_function);
    }

    /**
     * @brief Creates a new <b>Register</b> instance.  Upon creation, the
     * <i>DerivedClass::createInstance()</i> instance creation function is registered
     * with the <i>derived_name</i>.
     *
     * @param derived_name	a <b>std::string</b> const reference to the derived
     *				class' name.
     *
     */
    Register(const std::string &derived_name, Signature create_instance)
      : m_function(create_instance)
    {
      UserPlugin<Creator, Signature>::instance().registerCreator(derived_name, m_function);
    }

    /**
     * @brief Creates a new <b>Register</b> instance.  Upon creation, the
     * <i>DerivedClass::createInstance()</i> instance creation function is registered
     * with the <i>derived_id</i> and "enum id " derived_id.
     *
     * @param derived_id	a <b>int</b> to the derived class' id.
     *
     */
    explicit Register(int derived_id)
      : m_function(XDerivedClass::createInstance)
    {
      UserPlugin<Creator, Signature>::instance().registerCreator(derived_id_name(derived_id), m_function);
    }

    /**
     * @brief Creates a new <b>Register</b> instance.  Upon creation, the
     * <i>DerivedClass::createInstance()</i> instance creation function is registered
     * with the <i>derived_id</i> and "enum id " derived_id.
     *
     * @param derived_id	a <b>int</b> to the derived class' id.
     *
     */
    Register(int derived_id, Signature create_instance)
      : m_function(create_instance)
    {
      UserPlugin<Creator, Signature>::instance().registerCreator(derived_id_name(derived_id), m_function);
    }

  private:
    Signature		m_function;			///< Place to hold function pointer
  };
};


/**
 * Class template <b>UserSubroutine</b> is a template used for the association of
 * user function to be registered and acquired via the UserSubroutine mechanism.  The
 * template traits enforces the signature matching of the user function and the its
 * usage.
 *
 * The registration requires a unique function name and is required to implement a
 * traits class with a Signature typedef which specifies the signature of the user
 * function and a static function named getUserSubroutineName() which returns a const
 * std::string reference to the user function's name.  This name must be unique across users
 * of the registry.  There is no way to enforce this programatically, so it is recommended
 * that a application prefix be attached to the function name.
 *
 */
template <class S>
class UserSubroutine
{
private:
  UserSubroutine();					///< Not implemented
  UserSubroutine(const UserSubroutine&);		///< Not implemented
  UserSubroutine &operator=(const UserSubroutine&);	///< Not implemented

public:
  typedef S Signature;					///< Subroutine call signature

  /**
   * @brief Member function <b>instance</b> returns the instance of the registry,
   * cast as a UserSubroutine registry..
   *
   * @return			a <b>UserSubroutine</b> reference to the registry
   *				singleton.
   */
  inline static UserSubroutine &instance() {
    return (UserSubroutine<Signature> &) (Registry::rootInstance());
  }

  /**
   * @brief Member function <b>registerFunction</b> registers the user function's
   * name with the specified user function function pointer.
   *
   * The base class name is determined by the BaseClass template argument's
   * getCreatorName() static member function.  The signature is defined
   * UserSubroutineTraits template argument's Signature typedef.
   *
   * @param function_name	a <b>std::string</b> const reference to the user
   *				function's name.
   *
   * @param function		a <b>Signature</b> function pointer to the user
   *				function.
   *
   */
  inline static void registerFunction(const std::string &function_name, Signature *function) {
    Registry::rootInstance().registerIt(std::make_pair(&typeid(Signature), function_name), (void *) function);
  }

  /**
   * @brief Member function <b>execute</b> returns the user function function
   * associated with the specified signature and derived_name.
   *
   * @param function_name	a <b>std::string</b> const reference to the user
   *				function's name.
   *
   * @throws			a <b>std::invalid_argument</b> exception is thrown
   *				if there is no user function registered for the
   *				specified name.
   *
   * @return			a <b>Signature</b> user function.
   */
  static Signature *execute(const std::string &function_name) {
    Signature *user_function = (Signature *) Registry::rootInstance().getFunctionPtr(std::make_pair(&typeid(Signature), function_name));

    return (*user_function);
  }

  /**
   * @brief Member function <b>execute</b> returns the user function function
   * pointer associated with the specified signature and derived_name.
   *
   * @param function_name	a <b>std::string</b> const reference to the user
   *				function's name.
   *
   * @throws			a <b>std::invalid_argument</b> exception is thrown
   *				if there is no user function registered for the
   *				specified name.
   *
   * @return			a <b>Signature</b> user function pointer.
   */
  static Signature *getFunction(const std::string &function_name) {
    Signature *user_function = (Signature *) Registry::rootInstance().getFunctionPtr(std::make_pair(&typeid(Signature), function_name));

    return user_function;
  }

  /**
   * @brief Member function <b>exists</b> returns true if user function specified by
   * derived_name exists.
   *
   * @param function_name	a <b>std::string</b> const reference to the user
   *				function's name.
   *
   * @return			a <b>bool</b> of true if user function specified
   *				signature and <i>function_name</i> exists in BaseClass.
   */
  static bool exists(const std::string &derived_name) {
    return Registry::rootInstance().getFuncPtr(std::make_pair(&typeid(Signature), derived_name)) != NULL;
  }

  /**
   * @brief Class template <b>Register</b> registers the user function function
   * pointer with the <i>function_name</i> on object creation.
   *
   */
  class Register
  {
  public:
    /**
     * @brief Creates a new <b>Register</b> instance.  Upon creation, the
     * <i>func_ptr()</i> function is registered with the <i>function_name</i>.
     *
     * @param function_name	a <b>std::string</b> const reference to the user
     *				function's name.
     *
     */
    Register(const std::string &function_name, Signature *function)
      : m_function(function)
    {
      UserSubroutine<Signature>::instance().registerFunction(function_name, *m_function);
    }

  private:
    Signature *		m_function;			///< Holder of the function pointer
  };
};

template <>
void *Registry::getsym<void *>(const char *sym);

template <typename T>
inline T Registry::getsym(const char *sym) {
  return static_cast<T>(getsym<void *>(sym));
}

} // namespace Plugin
} // namespace sierra

typedef std::type_info *type_info_func();

/**
 * @brief FORTRAN compatible user subprogram registration routine.
 *
 * @par Description:
 * The first argument is for an application declared subprogram
 * that provides an example interface.  This example subprogram
 * is used to "type" the user-subprogram, e.g. the caller of
 * the registration routine guarantees that the interface of
 * the 'user_sub' exactly matches the interface of the 'type_sub'.
 */
extern "C" {
  void SIERRA_FORTRAN(register_user_subroutine)(
    type_info_func		type_id,
    void *			user_subroutine,
    const char *		name,
    int				name_length );
}


/**
 * Macro <b>FORTRAN_USER_SUBROUTINE</b> generates a FortranFunctionTraits template
 * specialization for the <i>RETURN</i> and <i>SIGNATURE</i> and creates a typedef
 * referencing the user function factory of <i>NAME</i>.
 *
 * Note that the user function has extern "C" linkage.
 *
 * @param NAME			name to be used by the Fortran EXTERNAL statement
 *
 * @param USER_SUB		user subroutine factory
 *
 */
#define FORTRAN_USER_SUBROUTINE(NAME, USER_SUB) extern "C" const std::type_info * SIERRA_FORTRAN(NAME)() {return &typeid(USER_SUB::Signature);}

#endif // STK_UTIL_DIAG_UserPlugin_h
