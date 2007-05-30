
#include "Teuchos_ParameterEntry.hpp"

// Boost initialization
#include <boost/python.hpp>
using namespace boost::python;


typedef Teuchos::ParameterEntry pentry;
void expose_parameter_entry()
{

    class_<pentry> ( "ParameterEntry" ,
        "This object is held as the \"value\" in the Teuchos::ParameterList map."
        " \n"
        "This structure holds a Teuchos::any value and information on the status of this \n"
        "parameter (isUsed, isDefault, etc.).  The type of parameter is chosen through the \n "
        "templated Set/Get methods. \n",
        init<>()
        )
        
        //   ParameterEntry(const ParameterEntry& source);
        .def( init<pentry>("Copy constructor") )
    
        // 
        //   //! Templated constructor
        //   template<typename T>
        //   explicit ParameterEntry(
        //     T value, bool isDefault = false, bool isList = false,
        //     const std::string &docString = "",
        //     RefCountPtr<const ParameterEntryValidator> const& validator = null
        //     );
        // 
        //   //! Destructor
        //   ~ParameterEntry();
        // 
        //   //@}
        // 
        //   //! @name Set Methods 
        //   //@{
        // 
        //   //! Replace the current parameter entry with \c source.
        //   ParameterEntry& operator=(const ParameterEntry& source);
        // 
        //   /*! \brief Templated set method that uses the input value type to determine the type of parameter.  
        //       
        //       \note <ul>
        // 	    <li> Invalidates any previous values stored by this object although it doesn't necessarily erase them.  
        //             <li> Resets 'isUsed' functionality.  
        // 	    </ul>
        //   */
        //   template<typename T>
        //   void setValue(
        //     T value, bool isDefault = false,
        //     const std::string &docString = "",
        //     RefCountPtr<const ParameterEntryValidator> const& validator = null
        //     );
        // 
        //   /*! \brief Set the value as an any object.
        //   *
        //   * This wipes all other data including documentation strings.
        //   *
        //   * Warning! Do not use function ths to set a sublist!
        //   */
        //   void setAnyValue(
        //     const any &value, bool isDefault = false
        //     );
        // 
        //   /*! \brief Set the validator. */
        //   void setValidator(
        //     RefCountPtr<const ParameterEntryValidator> const& validator
        //     );
        // 
        //   /*! \brief Set the documentation string. */
        //   void setDocString(const std::string &docString);
        // 
        //   //! Create a parameter entry that is an empty list.
        //   ParameterList& setList(
        //     bool isDefault = false,
        //     const std::string &docString = ""
        //     );
        // 
        //   //@}
        // 
        //   //! @name Get Methods 
        //   //@{
        //    
        //   /*! \brief Templated get method that uses the input pointer type to determine the type of parameter to return.  
        // 
        //       \note This method will cast the value to the type requested.  If that type is incorrect, 
        // 	    an exception will be thrown by the any_cast.
        //   */
        //   template<typename T>
        //   T& getValue(T *ptr) const;
        // 
        //   /*! \brief Direct access to the Teuchos::any data value underlying this
        //    *  object. The bool argument \c activeQry (default: true) indicates that the 
        //    *  call to getAny() will set the isUsed() value of the ParameterEntry to true.
        //    */
        //   any& getAny(bool activeQry = true);
        // 
        //   /*! \brief Constant direct access to the Teuchos::any data value underlying this
        //    *  object. The bool argument \c activeQry (default: true) indicates that the 
        //    *  call to getAny() will set the isUsed() value of the ParameterEntry to true.
        //    */
        //   const any& getAny(bool activeQry = true) const;
        // 
        //   //@}
        // 
        //   //! @name Attribute/Query Methods 
        //   //@{
        //   
        //   //! Return whether or not the value has been used; i.e., whether or not the value has been retrieved via a get function.
        //   bool isUsed() const;
        // 
        //   //! Return whether or not the value itself is a list.
        //   bool isList() const;
        // 
        //   //! Test the type of the data being contained.
        //   template <typename T>
        //   bool isType() const;
        // 
        //   //! Indicate whether this entry takes on the default value.
        //   bool isDefault() const;
        // 
        //   //! Return the (optional) documentation string
        //   std::string docString() const;

  ;
}