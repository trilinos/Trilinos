#ifndef _TSF_PARAMETER_H_
#define _TSF_PARAMETER_H_


#include "TSF_Object.h"

namespace TSF {

//! TSF::Parameter:  The Trilinos Solver Framework Parameter class.
/*! The TSF::Parameter class encapsulates information about solver parameters.
    Any common type of data can be used to construct an TSF::Parameter object.
    There is also a constructor that accepts a string that can be parsed into 
    parameter object.
*/
class Parameter {
    
  public:
  //@{ \name Constructor/Destructor Methods
  //! TSF::Parameter default constructor.
  /*! Allows use of an array of TSF::Parameter objects. Creates a uninitialized TSF::Parameter instance that
      must be initialized using one of the set() methods.
  */
  Parameter(void);

  //! TSF::Parameter Copy Constructor.
  /*! Makes an exact copy of an existing TSF::Parameter instance.
  */
  Parameter(const TSF::Parameter& parameter);

  //! TSF::Parameter character constructor.
  Parameter(char const * const & label, char parameter);

  //! TSF::Parameter character string constructor.
  Parameter(char const * const & label, char * parameter);

  //! TSF::Parameter int constructor.
  Parameter(char const * const & label, int parameter);

  //! TSF::Parameter int array constructor.
  Parameter(char const * const & label, int length, int * parameter);

  //! TSF::Parameter double constructor.
  Parameter(char const * const & label, double parameter);

  //! TSF::Parameter double array constructor.
  Parameter(char const * const & label, int length,  double * parameter);

  //! TSF::Parameter Destructor.
  /*! Completely deletes a TSF::Parameter object.  
  */
  virtual ~Parameter(void);
  //@}

  //@{ \name Methods to set parameter values

  //! TSF::Parameter character set function.
  int set(char const * const & label, char parameter);

  //! TSF::Parameter character string set function.
  int set(char const * const & label, char * parameter);

  //! TSF::Parameter int set function.
  int set(char const * const & label, int parameter);

  //! TSF::Parameter int array set function.
  int set(char const * const & label, int length, int * parameter);

  //! TSF::Parameter double set function.
  int set(char const * const & label, double parameter);

  //! TSF::Parameter double array set function.
  int set(char const * const & label, int length,  double * parameter);

  //@}

  //@{ \name Methods to test parameter values
  //! True if parameter is character
  bool isChar() const { return(isChar_);};
  //! True if parameter is character string
  bool isCharString() const { return(isCharString_);};
  //! True if parameter is integer
  bool isInt() const { return(isInt_);};
  //! True if parameter is integer array
  bool isIntArray() const { return(isIntArray_);};
  //! True if parameter is double
  bool isDouble() const { return(isDouble_);};
  //! True if parameter is double array
  bool isDoubleArray() const { return(isDoubleArray_);};
  //@}

  //@{ \name Methods to extract parameter values
  //! Return character parameter (assumes isChar() has been called and is true).
  char     getChar() const {return(char_);};
  //! Return character string parameter (assumes isCharString() has been called and is true).
  char *   getCharString() const {return(charString_);};
  //! Return integer parameter (assumes isInt() has been called and is true).
  int      getInt() const {return(int_);};
  //! Return integer array parameter (assumes isIntArray() has been called and is true).
  int *    getIntArray() const {return(intArray_);};
  //! Return double parameter (assumes isDouble() has been called and is true).
  double   getDouble() const {return(double_);};
  //! Return double array parameter (assumes isDouble() has been called and is true).
  double * getDoubleArray() const {return(doubleArray_);};

  //@}

  //@{ \name Utility methods.
  //! operator= performs a deep copy of a TSF::Parameter object
  Parameter & operator = (const Parameter & parameter) {assign(parameter);return(*this);};

  //! TSF::Parameter label accessor function.
  /*! Returns the label associated with this parameter.  A label is a description of the parameter, used
      both to describe the variable and identify it to any prospective Trilinos solver component.
   */
  char * & getLabel() const {return(label_);};

  //@}

 private:

  void initializeDefaults(char const * const & label);
  void initializeDefaults();
  void assign(const TSF::Parameter & parameter);
  char * label_;
  char char_;
  char * charString_;
  int int_;
  int * intArray_;
  double double_;
  double * doubleArray_;
  int length_;
  
  bool isChar_, isCharString_, isInt_, isIntArray_, isDouble_, isDoubleArray_;
};
} // TSF namespace

#endif /* _TSF_PARAMETER_H_ */
