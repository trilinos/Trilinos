#ifndef _fei_ostream_ops_hpp_
#define _fei_ostream_ops_hpp_
/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_iosfwd.hpp"
namespace fei {
  class FillableMat;
  class CSRMat;
  class CSVec;
}//namespace fei

/*
 * This header is ONLY to be included from within other FEI headers or sources.
 * The macro FEI_OSTREAM must already be defined. (It is defined in
   fei_iostream.hpp and/or in fei_iosfwd.hpp.)
 */

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::Vector& vec);

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::Matrix& mat);

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::FillableMat& mat);

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::CSRMat& mat);

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::CSVec& vec);

#endif // _fei_ostream_ops_hpp_
