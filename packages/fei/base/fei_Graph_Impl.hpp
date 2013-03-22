/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_Graph_Impl_hpp_
#define _fei_Graph_Impl_hpp_

#include <fei_iosfwd.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_Graph.hpp>
#include <fei_mpi.h>
#include <fei_EqnComm.hpp>

namespace fei {

  /** Basic algebraic matrix-graph representation. */
  class Graph_Impl : public fei::Graph {
  public:
    /** constructor */
    Graph_Impl(MPI_Comm comm, int firstLocalRow, int lastLocalRow);

    /** destructor */
    virtual ~Graph_Impl();

    /** Add indices to a specified row of the table */
    int addIndices(int row,
		   int len,
		   const int* indices);

    /** Add a symmetric block of indices. The array of indices will serve as
	both row-numbers, and as column-numbers in those rows.
    */
    int addSymmetricIndices(int numIndices,
			    int* indices,
			    bool diagonal=false);

    /** gather all remotely-owned table portions to owning processors */
    int gatherFromOverlap();

    /** Retrieve the local portion of the graph. i.e., The rows which
       correspond to locally-owned IDs. */
    table_type* getLocalGraph()
      {
	return( localGraphData_ );
      }

    /** Retrieve the remotely-owned portion of the graph. */
    std::vector<remote_table_type*>& getRemoteGraph()
      {
	return( remoteGraphData_);
      }

    /** Write locally-owned portion of the graph to a specified ostream.*/
    int writeLocalGraph(FEI_OSTREAM& os,
			bool debug=false,
			bool prefixLinesWithPoundSign=true);

    /** Write remotely-owned portion of the graph to a specified ostream.*/
    int writeRemoteGraph(FEI_OSTREAM& os);

    /** Get the number of locally-owned rows. */
    int getNumLocalRows();

    /** Get the number of nonzeros in the locally-owned portion of
       the graph. */
    int getNumLocalNonzeros() const;

    /** Get the length of a specified locally-owned row. */
    int getLocalRowLength(int row);

  private:
    void addDiagonals(int numIndices, int* indices);

    table_type* localGraphData_;
    std::vector<remote_table_type*> remoteGraphData_;
    fei::SharedPtr<fei::EqnComm> eqnComm_;

    int firstLocalRow_, lastLocalRow_;
    int localProc_, numProcs_;
    MPI_Comm comm_;
  };//class Graph_Impl

} //namespace fei

#endif // _fei_Graph_Impl_hpp_

