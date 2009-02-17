#ifndef _fei_fwd_hpp_
#define _fei_fwd_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_defs.h>

//
//Forward declarations for fei classes.
//

//
//First, the "old" classes that aren't in a namespace.
//

class BlockDescriptor;
class ConnectivityTable;
class Data;
class EqnBuffer;
class EqnCommMgr;
class FEI_Implementation;
class Filter;
class FiniteElementData;
class LibraryWrapper;
class LinearSystemCore;
class NodeCommMgr;
class NodeDatabase;
class NodeDescriptor;
class ProcEqns;
class SNL_FEI_Structure;
class SlaveVariable;

//
//Now the symbols in the fei namespace.
//

namespace fei {
  /** enumeration for various output levels */
  enum OutputLevel {
    //We want to make sure that BRIEF_LOGS < FULL_LOGS < ALL, so that
    //debug-output code can use statements like
    //  'if (output_level_ >= BRIEF_LOGS)'.
    //So be aware that re-arranging these values can change the behavior
    //of debug-output code... Don't do it!
    NONE = 0,
    STATS = 1,
    MATRIX_FILES = 2,
    BRIEF_LOGS = 3,
    FULL_LOGS = 4,
    ALL = 5
  };

  class Factory;
  class FillableMat;
  class Graph;
  class CSRMat;
  class CSVec;
  class LogFile;
  class Logger;
  class LogManager;
  class Reducer;
  class VectorSpace;
  class MatrixGraph;
  class Param;
  class ParameterSet;
  class SparseRowGraph;
  class Vector;
  class Matrix;
  class LinearSystem;
  template<typename T> class ctg_set;
  template<typename T> class Matrix_Impl;
  template<typename T> class Vector_Impl;
}//namespace fei

//
//Finally the symbols that are still in the soon-to-be-eliminated
//snl_fei namespace.
//
namespace snl_fei {
  template<class RecordType> class Constraint;
  class RecordCollection;

  class BlockDescriptor;
  class PointBlockMap;

  class Broker;
  class Factory;
}

#undef FEI_OSTREAM

#ifdef FEI_HAVE_IOSFWD
#include <iosfwd>
#define FEI_OSTREAM std::ostream
#else
#include <iostream.hpp>
#define FEI_OSTREAM ostream
#endif

#include <exception>
#include <vector>

#endif

