/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#ifndef _fei_fwd_hpp_
#define _fei_fwd_hpp_


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
  template<typename T> class SharedIDs;
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

