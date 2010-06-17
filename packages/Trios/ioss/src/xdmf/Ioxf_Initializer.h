/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef IOSS_Ioxf_Initializer_h
#define IOSS_Ioxf_Initializer_h

namespace Ioxf {
  class Initializer
    {
    public:
      Initializer();
      ~Initializer();
      // Copy constructor
      // Assignment operator
    private:
      static int useCount;
    };
} 
#endif
