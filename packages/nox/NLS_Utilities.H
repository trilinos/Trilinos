
// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_UTILITIES_H
#define _NLS_UTILITIES_H
  
#include "NLS_ParameterList.H"

//! NLSPACK pure abstract class to provide parallel communication 
//! and printing utilities.

class NLS_Utilities {

public:

  //! Constructor
  NLS_Utilities();

  //! Destructor
  ~NLS_Utilities();

  //! Set the parameters in the utilities class.
  void setUtilities(NLS_ParameterList& p);

  //! Returns true if the current processor is designated as the processor 
  //! send output to screen and/or file.
  /*  This is used to prevent each processor from printing the same 
   *  information. */
  virtual bool isPrintProc();

  //! Returns true if printLevel is less than or equal to the
  //! "Output Level" flag.  
  /*! This is used to determine if an output line should be printed.
   *  The user sets an outputLevel in the ParameterList for the MethodManager 
   *  If the print level is less than the output level, this function 
   *  returns true.  The most information is printed at outputLevel 5.  
   *  The least at level 0. */ 
  virtual bool isOutput(int printLevel);

  //! Returns the processor ID number for parallel jobs. Returns 0 
  //! for serial jobs.
  /*! Returns the rank of calling process */
  virtual int getMyPID();
  
protected:
  //! Processor ID.  Defaults to 0 for serial.
  int myPID;
  
  //! Output Level. Defaults to level 2.
  int outputLevel;

  //! Processor responsible for output. Defaults to processor 0.
  int printProc;
};

#endif
