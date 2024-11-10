/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef _IFPACK_ILUK_GRAPH_H_
#define _IFPACK_ILUK_GRAPH_H_

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_ConfigDefs.h"
#include "Epetra_Object.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Import.h"

#include "Teuchos_RefCountPtr.hpp"

namespace Teuchos {
  class ParameterList;
}

//! Ifpack_IlukGraph: A class for constructing level filled graphs for use with ILU(k) class preconditioners.

/*! The Ifpack_IlukGraph class enable the construction matrix graphs using level-fill algorithms.  The only function
    required for construction is an ExtractRowView capability, i.e., the matrix that is passed in to the constructor
    must implement the Ifpack_CrsGraph interface defined in Ifpack_CrsMatrix.h


<b>Constructing Ifpack_IlukGraph objects</b>

Constructing Ifpack_IlukGraph objects is usually a two step process of passing in a Ifpack_CrsGraph object and
an integer indicating the desired level of fill and then calling the ConstructFilledGraph function to complete the
process.  This allows warning error codes to be returned to the calling routine.

It is worth noting that an Ifpack_IlukGraph object has two
Epetra_CrsGraph objects containing L and U, the graphs for the lower and upper triangular parts of the ILU(k) graph.
Thus, it is possible to manually insert and delete graph entries in L and U via the Epetra_CrsGraph
InsertIndices and RemoveIndices functions.  However, in this case FillComplete must be
called before the graph is used for subsequent operations.

*/
class Ifpack_IlukGraph {

  // Give ostream << function some access to private and protected data/functions.

  friend std::ostream& operator << (std::ostream& os, const Ifpack_IlukGraph& A);

 public:

  //! Ifpack_IlukGraph constuctor.
  /*! Creates a Ifpack_IlukGraph object using the input graph and specified level of fill.

    \param In
           Graph_in - An existing Ifpack_CrsGraph.  This object must implement the Ifpack_CrsGraph functions
           that provide graph dimension and pattern information.
    \param In
           LevelFill_in - The level of fill to compute via ILU(k) algorithm.
    \param In
           LevelOverlap_in - The level of between subdomains.

           \warning Actual construction occurs in ConstructFilledGraph.  This allows error codes to
                    be passed back to the user.
  */
  Ifpack_IlukGraph(const Epetra_CrsGraph & Graph_in, int LevelFill_in, int LevelOverlap_in);

  //! Copy constructor.
  Ifpack_IlukGraph(const Ifpack_IlukGraph & Graph_in);

  //! Ifpack_IlukGraph Destructor
  virtual ~Ifpack_IlukGraph();

  //!Set parameters using Teuchos::ParameterList object.
  /* This method is only available if the Teuchos package is enabled.
     This method recogizes two parameter names: Level_fill and Level_overlap.
     Both are case insensitive, and in both cases the ParameterEntry must
     have type int.
  */
  int SetParameters(const Teuchos::ParameterList& parameterlist,
                    bool cerr_warning_if_unused=false);

  //! Does the actual construction of the graph.
  /*
    \return Integer error code, set to 0 if successful.

  */
  virtual int ConstructFilledGraph();

  //! Does the actual construction of the overlap matrix graph.
  /*
    \return Integer error code, set to 0 if successful.

  */
  virtual int ConstructOverlapGraph();

  //! Returns the level of fill used to construct this graph.
  virtual int LevelFill() const {return(LevelFill_);};

  //! Returns the level of overlap used to construct this graph.
  virtual int LevelOverlap() const {return(LevelOverlap_);};

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  //! Returns the number of global matrix rows.
  int NumGlobalBlockRows() const {
    if(Graph_.RowMap().GlobalIndicesInt())
      return (int)(NumGlobalBlockRows_);
    else
      throw "Ifpack_IlukGraph::NumGlobalBlockRows: GlobalIndices not int.";
  }

  //! Returns the number of global matrix columns.
  int NumGlobalBlockCols() const {
    if(Graph_.RowMap().GlobalIndicesInt())
      return (int)(NumGlobalBlockCols_);
    else
      throw "Ifpack_IlukGraph::NumGlobalBlockCols: GlobalIndices not int.";
  }

  //! Returns the number of global matrix rows.
  int NumGlobalRows() const {
    if(Graph_.RowMap().GlobalIndicesInt())
      return (int)(NumGlobalRows_);
    else
      throw "Ifpack_IlukGraph::NumGlobalRows: GlobalIndices not int.";
  }

  //! Returns the number of global matrix columns.
  int NumGlobalCols() const {
    if(Graph_.RowMap().GlobalIndicesInt())
      return (int)(NumGlobalCols_);
    else
      throw "Ifpack_IlukGraph::NumGlobalCols: GlobalIndices not int.";
  }
  //! Returns the number of nonzero entries in the global graph.
  int NumGlobalNonzeros() const {
    if(Graph_.RowMap().GlobalIndicesInt())
      return (int)(NumGlobalNonzeros_);
    else
      throw "Ifpack_IlukGraph::NumGlobalNonzeros: GlobalIndices not int.";
  }

  //! Returns the number of diagonal entries found in the global input graph.
  virtual int NumGlobalBlockDiagonals() const {
    if(Graph_.RowMap().GlobalIndicesInt())
      return (int)(NumGlobalBlockDiagonals_);
    else
      throw "Ifpack_IlukGraph::NumGlobalBlockDiagonals: GlobalIndices not int.";
  }
#endif

  //! Returns the number of global matrix rows.
  long long NumGlobalBlockRows64() const {return(NumGlobalBlockRows_);};

  //! Returns the number of global matrix columns.
  long long NumGlobalBlockCols64() const {return(NumGlobalBlockCols_);};

  //! Returns the number of global matrix rows.
  long long NumGlobalRows64() const {return(NumGlobalRows_);};

  //! Returns the number of global matrix columns.
  long long NumGlobalCols64() const {return(NumGlobalCols_);};
  //! Returns the number of nonzero entries in the global graph.
  long long NumGlobalNonzeros64() const {return(NumGlobalNonzeros_);};

  //! Returns the number of diagonal entries found in the global input graph.
  virtual long long NumGlobalBlockDiagonals64() const {return(NumGlobalBlockDiagonals_);};

  //! Returns the number of local matrix rows.
  int NumMyBlockRows() const {return(NumMyBlockRows_);};

  //! Returns the number of local matrix columns.
  int NumMyBlockCols() const {return(NumMyBlockCols_);};


  //! Returns the number of local matrix rows.
  int NumMyRows() const {return(NumMyRows_);};

  //! Returns the number of local matrix columns.
  int NumMyCols() const {return(NumMyCols_);};

  //! Returns the number of nonzero entries in the local graph.
  int NumMyNonzeros() const {return(NumMyNonzeros_);};

  //! Returns the number of diagonal entries found in the local input graph.
  virtual int NumMyBlockDiagonals() const {return(NumMyBlockDiagonals_);};

  //! Returns the index base for row and column indices for this graph.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int IndexBase() const {
    if(Graph_.RowMap().GlobalIndicesInt())
      return (int) IndexBase64();
    throw "Ifpack_IlukGraph::IndexBase: GlobalIndices not int.";
  }
#endif
  long long IndexBase64() const {return(IndexBase_);};

  //! Returns the graph of lower triangle of the ILU(k) graph as a Epetra_CrsGraph.
  virtual Epetra_CrsGraph & L_Graph() {return(*L_Graph_);};

  //! Returns the graph of upper triangle of the ILU(k) graph as a Epetra_CrsGraph.
  virtual Epetra_CrsGraph & U_Graph() {return(*U_Graph_);};

  //! Returns the graph of lower triangle of the ILU(k) graph as a Epetra_CrsGraph.
  virtual Epetra_CrsGraph & L_Graph() const {return(*L_Graph_);};

  //! Returns the graph of upper triangle of the ILU(k) graph as a Epetra_CrsGraph.
  virtual Epetra_CrsGraph & U_Graph() const {return(*U_Graph_);};

  //! Returns the importer used to create the overlapped graph.
  virtual Epetra_Import * OverlapImporter() const  {return(&*OverlapImporter_);};

  //! Returns the the overlapped graph.
  virtual Epetra_CrsGraph * OverlapGraph() const  {return(&*OverlapGraph_);};

    //! Returns the Epetra_BlockMap object associated with the domain of this matrix operator.
    virtual const Epetra_BlockMap & DomainMap() const {return(DomainMap_);};

    //! Returns the Epetra_BlockMap object associated with the range of this matrix operator.
    virtual const Epetra_BlockMap & RangeMap() const{return(RangeMap_);};

    //! Returns the Epetra_BlockMap object associated with the range of this matrix operator.
    virtual const Epetra_Comm & Comm() const{return(Comm_);};

 private:


  const Epetra_CrsGraph & Graph_;
  const Epetra_BlockMap & DomainMap_;
  const Epetra_BlockMap & RangeMap_;
  const Epetra_Comm & Comm_;
  Teuchos::RefCountPtr<Epetra_CrsGraph> OverlapGraph_;
  Teuchos::RefCountPtr<Epetra_BlockMap> OverlapRowMap_;
  Teuchos::RefCountPtr<Epetra_Import> OverlapImporter_;
  int LevelFill_;
  int LevelOverlap_;
  Teuchos::RefCountPtr<Epetra_CrsGraph> L_Graph_;
  Teuchos::RefCountPtr<Epetra_CrsGraph> U_Graph_;
  long long IndexBase_;
  long long NumGlobalRows_;
  long long NumGlobalCols_;
  long long NumGlobalBlockRows_;
  long long NumGlobalBlockCols_;
  long long NumGlobalBlockDiagonals_;
  long long NumGlobalNonzeros_;
  long long NumGlobalEntries_;
  int NumMyBlockRows_;
  int NumMyBlockCols_;
  int NumMyRows_;
  int NumMyCols_;
  int NumMyBlockDiagonals_;
  int NumMyNonzeros_;
  int NumMyEntries_;


 };

//! << operator will work for Ifpack_IlukGraph.
std::ostream& operator << (std::ostream& os, const Ifpack_IlukGraph& A);

#endif /* _IFPACK_ILUK_GRAPH_H_ */
