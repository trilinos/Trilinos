/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_iosfwd_hpp_
#define _fei_iosfwd_hpp_

#include <fei_fwd.hpp>

#undef FEI_OSTREAM

#ifdef FEI_HAVE_IOSFWD
#include <iosfwd>
#define FEI_OSTREAM std::ostream
#elif defined FEI_HAVE_IOSTREAM
#include <iostream>
#define FEI_OSTREAM std::ostream
#else
#include <iostream.h>
#define FEI_OSTREAM ostream
#endif

#endif

