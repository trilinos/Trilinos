/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_GraphReducer_hpp_
#define _fei_GraphReducer_hpp_

#include <fei_iosfwd.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_Graph.hpp>
#include <fei_CommUtils.hpp>
#include <fei_Reducer.hpp>
#include <snl_fei_RaggedTable_specialize.hpp>

namespace fei {

  /** Basic algebraic matrix-graph representation. */
  class GraphReducer : public fei::Graph {
  public:
    /** constructor */
    GraphReducer(fei::SharedPtr<fei::Reducer> reducer,
                 fei::SharedPtr<fei::Graph> target);

    /** destructor */
    virtual ~GraphReducer();

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
	return( target_->getLocalGraph());
      }

    /** Retrieve the remotely-owned portion of the graph. */
    std::vector<remote_table_type*>& getRemoteGraph()
      {
	return( target_->getRemoteGraph());
      }

    /** Write locally-owned portion of the graph to a specified ostream.*/
    int writeLocalGraph(FEI_OSTREAM& os,
			bool debug=false,
			bool prefixLinesWithPoundSign=true);

    /** Write remotely-owned portion of the graph to a specified ostream.*/
    int writeRemoteGraph(FEI_OSTREAM& os);

  private:
    fei::SharedPtr<fei::Reducer> reducer_;
    fei::SharedPtr<fei::Graph> target_;
  };//class GraphReducer

} //namespace fei

#endif // _fei_GraphReducer_hpp_

