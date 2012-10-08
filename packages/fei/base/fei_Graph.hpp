/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_Graph_hpp_
#define _fei_Graph_hpp_

#include <fei_iosfwd.hpp>
#include <snl_fei_RaggedTable_specialize.hpp>

#include <map>

namespace fei {

  /** Basic algebraic matrix-graph representation. */
  class Graph {
  public:

    /** destructor */
    virtual ~Graph(){}

    /** alias for the 'table_type' data container */
    typedef snl_fei::RaggedTable<snl_fei::MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >
      table_type;

    /** alias for table_row_type, which is a row of the table */
    typedef fei::ctg_set<int> table_row_type;

    /** alias for the type of the remotely-owned portion of the table data */
    typedef snl_fei::RaggedTable<std::map<int,fei::ctg_set<int>*>,fei::ctg_set<int> > 
      remote_table_type;

    /** Add indices to a specified row of the table */
    virtual int addIndices(int row,
		   int len,
		   const int* indices) = 0;

    /** Add a symmetric block of indices. The array of indices will serve as
	both row-numbers, and as column-numbers in those rows.
    */
    virtual int addSymmetricIndices(int numIndices,
			    int* indices,
			    bool diagonal=false) = 0;

    /** gather all remotely-owned table portions to owning processors */
    virtual int gatherFromOverlap() = 0;

    /** Retrieve the local portion of the graph. i.e., The rows which
       correspond to locally-owned IDs. */
    virtual table_type* getLocalGraph() = 0;

    /** Retrieve the remotely-owned portion of the graph. */
    virtual std::vector<remote_table_type*>& getRemoteGraph() = 0;

    /** Write locally-owned portion of the graph to a specified ostream.*/
    virtual int writeLocalGraph(FEI_OSTREAM& os,
			bool debug=false,
			bool prefixLinesWithPoundSign=true) = 0;

    /** Write remotely-owned portion of the graph to a specified ostream.*/
    virtual int writeRemoteGraph(FEI_OSTREAM& os) = 0;

  };//class Graph

} //namespace fei

#endif // _fei_Graph_hpp_

