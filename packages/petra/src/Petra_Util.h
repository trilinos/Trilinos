#ifndef _PETRA_UTIL_H_
#define _PETRA_UTIL_H_

//! Petra_Util:  The Petra Util Wrapper Class.
/*! The Petra_Util class is a collection of useful functions that cut across a broad
    set of other classes.  Specifically, sorting is provided by this class.

    Petra_Util is a serial interface only.  This is appropriate since the standard 
    utilities are only specified for serial execution (or shared memory parallel).
*/

#include "Petra_Petra.h"

class Petra_Util {
    
  public:
  //! Petra_Util Constructor.
  /*! Builds an instance of a serial Util object.
   */
  Petra_Util(void);


  //! Petra_Util Copy Constructor.
  /*! Makes an exact copy of an existing Petra_Util instance.
  */
  Petra_Util(const Petra_Util& Util);

  //! Petra_Util Destructor.
  virtual ~Petra_Util(void);
  
  //! Petra_Util Sort Routine (Shell sort)
  /*! 

    This function sorts a list of integer values in ascending or descending order.  Additionally it sorts any
    number of companion lists of doubles or ints.  A shell sort is used, which is fast if indices are already sorted.
    
    \param In
           SortAscending - Sort keys in ascending order if true, otherwise sort in descending order..
    \param In
           NumKeys - Number of integer values to be sorted.
    \param In/Out
           Keys - List of integers to be sorted.
    \param In
           NumDoubleCompanions - Number of lists of double precision numbers to be sorted with the key.
    \param In
           DoubleCompanions - DoubleCompanions[i] is a pointer to the ith list of doubles to be sorted with key.
    \param In
           NumIntCompanions - Number of lists of integers to be sorted with the key.
    \param In
           IntCompanions - IntCompanions[i] is a pointer to the ith list of integers to be sorted with key.
	   
  */
  void Sort(bool SortAscending, int NumKeys, int * Keys, int NumDoubleCompanions,double ** DoubleCompanions, 
            int NumIntCompanions, int ** IntCompanions) const;
};

// Petra_Util constructor
inline Petra_Util::Petra_Util(void){}
// Petra_Util constructor
inline Petra_Util::Petra_Util(const Petra_Util& Util){}
// Petra_Util destructor
inline Petra_Util::~Petra_Util(){}

#endif /* _PETRA_UTIL_H_ */
