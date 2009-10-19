#ifndef STK_UTIL_UTIL_BOOTSTRAP_HPP
#define STK_UTIL_UTIL_BOOTSTRAP_HPP

namespace stk {

///
/// @addtogroup bootstrap_detail
/// @{
///

/**
 * @brief Class <b>Bootstrap</b> serves as a bootstrapping mechanism for products in the
 * sierra toolkit and elsewhere.
 *
 * Often, it is convenient to have a product perform registrations and other operations
 * when linked into an application.  One method of accomplishing this is to utilize a
 * static object whose constructor perform these operations.  However, static
 * constructions are executed prior to main and in non-deterministc order.  The later is
 * particularly problematic if the constructor results in the usage of another static
 * object which may not have been constructed yet.
 *
 * So, the Bootstrap class creates a stack of callback functions that are executed, by
 * main(), when Bootstrap::bootstrap() is called.  These functions are still executed in a
 * non-deterministic order, but all static constructions have occurred.
 *
 */
class Bootstrap
{
public:

  typedef void (*FunctionPtr)();

  /**
   * @brief Member function <b>bootstrap</b> runs through the stored bootstrap function
   * pointers and executes each function.
   *
   */
  static void bootstrap();
  
  /**
   * @brief Creates a new <b>Bootstrap</b> instance.
   *
   * The instance serves only to insert the specified function pointer to the bootstrap
   * function list.  If the bootstrapper has already been executed, the function is added
   * to the list and executed immediately.
   *
   */
  Bootstrap(void (*f)());

private:
  Bootstrap(const Bootstrap &);
  Bootstrap &operator=(const Bootstrap &);

public:
  ~Bootstrap()
  {}

private:
  static Bootstrap *    s_front;                        ///< Front of bootstrap registry
  static bool           s_bootstrapped;                 ///< Bootstart::bootstrap has been called
  
  Bootstrap *           m_next;                         ///< Next bootstrap object
  FunctionPtr           m_f;                            ///< Function to call during bootstrap()
};

///
/// @}
///

} // namespace stk

#endif // STK_UTIL_UTIL_BOOTSTRAP_HPP
