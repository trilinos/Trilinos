#ifndef _TSF_PARAMETERLIST_H_
#define _TSF_PARAMETERLIST_H_


#include "TSF_Object.h"
#include "TSF_Parameter.h"

namespace TSF {
//! TSF::ParameterList:  The Trilinos Solver Framework ParameterList class.
/*! The TSF::ParameterList class encapsulates information about solver parameters.
    Any common type of data can be used to construct an TSF::ParameterList object.
    There is also a constructor that accepts a string that can be parsed into 
    parameter object.
*/
class ParameterList {
    
  public:

  //@{ \name Constructor/Destructor Methods
  //! TSF::ParameterList default constructor.
  /*! Allows use of an array of TSF::ParameterList objects. Creates a uninitialized TSF::ParameterList instance that
      must be initialized using one of the addParameter() methods.
  */
  ParameterList(void);

  //! TSF::Parameter Constructor; takes a list of TSF::Parameters.
  /*! Since ParameterList is a collection of Parameters, we allow one to submit multiple parameters simultaneously.
  */
  ParameterList(int numParameters, const TSF::Parameter *& parameters);

  //! TSF::ParameterList Copy Constructor.
  /*! Makes an exact copy of an existing TSF::ParameterList instance.
  */
  ParameterList(const TSF::ParameterList& parameterList);

  //! TSF::ParameterList Destructor.
  /*! Completely deletes a TSF::ParameterList object.  
  */
  virtual ~ParameterList(void);
  //@}

  //@{ \name Methods to set parameters to a ParameterList

  //! TSF::ParameterList set TSF::Parameter class instance.
  void setParameter(const TSF::Parameter parameter);

  //! TSF::ParameterList character set function.
  void setParameter(char const * const & key, bool parameter);

  //! TSF::ParameterList character set function.
  void setParameter(char const * const & key, char parameter);

  //! TSF::ParameterList character string set function.
  void setParameter(char const * const & key, char * parameter);

  //! TSF::ParameterList int set function.
  void setParameter(char const * const & key, int parameter);

  //! TSF::ParameterList int array set function.
  void setParameter(char const * const & key, int length, int * parameter);

  //! TSF::ParameterList double set function.
  void setParameter(char const * const & key, double parameter);

  //! TSF::ParameterList double array set function.
  void setParameter(char const * const & key, int length,  double * parameter);
  //@}

  //@{ \name Methods to query a ParameterList for a specific parameter

  //! TSF::ParameterList Parameter object extraction function.
  /*! Returns true if parameter is a char string, in which case, the argument parameter is returned.
   */
  bool hasParameter(char const * const & key) const;

  //! TSF::ParameterList character query function.
  /*! Returns true if the key matches a character parameter in the list.
   */
  bool hasChar(char const * const & key) const;

  //! TSF::ParameterList char string query function.
  /*! Returns true if key matches a char string in the list.
   */
  bool hasCharArray(char const * const & key) const;

  //! TSF::ParameterList integer query function.
  /*! Returns true if key matches an integer in the list.
   */
  bool hasInt(char const * const & key) const;

  //! TSF::ParameterList integer array query function.
  /*! Returns true if key matches an integer array in the list.
   */
  bool hasIntArray(char const * const & key) const;

  //! TSF::ParameterList double query function.
  /*! Returns true if key matches a double in the list.
   */
  bool hasDouble(char const * const & key) const;

  //! TSF::ParameterList double array query function.
  /*! Returns true if key matches a double array in the list.
   */
  bool hasDoubleArray(char const * const & key) const;
  //@}

  //@{ \name Methods to access parameters in a ParameterList

  //! TSF::ParameterList Parameter object access function.
  /*! Returns parameter corresponding to specified key.
   */
  Parameter getParameter(char const * const & key) const;

  //! TSF::ParameterList char string access function.
  /*! Returns char corresponding to specified key.
   */
  char getChar(char const * const & key) const;

  //! TSF::ParameterList char string access function.
  /*! Returns character array corresponding to specified key.
   */
  char * getCharArray(char const * const & key) const;

  //! TSF::ParameterList integer access function.
  /*! Returns integer corresponding to specified key.
   */
  int getInt(char const * const & key) const;

  //! TSF::ParameterList integer array access function.
  /*! Returns integer array corresponding to specified key,
      along with its length (on argument list).
   */
  int * getIntArray(char const * const & key, int & length) const;

  //! TSF::ParameterList double access function.
  /*! Returns double corresponding to specified key.
   */
  double getDouble(char const * const & key) const;

  //! TSF::ParameterList double array access function.
  /*! Returns double array corresponding to specified key,
      along with its length (on argument list).
   */
  double* getDoubleArray(char const * const & key, int & length) const;
  //@}

  //@{ \name Individual TSF::Parameter access methods
  //! Returns the number of parameters currently contained in this TSF::ParameterList object.
  /*! Use of this method and the overloaded operator[] method allows a user to access any particular
      TSF::Parameter in the TSF::ParameterList, or to process all TSF::Parameter objects, one at a time.
  */
  int getNumParameters(void) const {return(numParameters_);};
  //! Returns a specific TSF::Parameter in the TSF::ParameterList
  TSF::Parameter & operator[] (int i) const {if (i < numParameters_) return(parameters_[i]); else exit(-1);};
  //@}
 private:

  void initializeDefaults(void);
  void deleteArrays(void);
  void enlargeParameters(void);
  int getMaxNumParameters(void) const {return(maxNumParameters_);};
  int getKeyLoc(char const * const & key) const;
  int setKeyLoc(char const * const & key) const;
  int numParameters_;
  int maxNumParameters_;

  TSF::Parameter * parameters_;
  
};
} // TSF namespace
#endif /* _TSF_PARAMETERLIST_H_ */
