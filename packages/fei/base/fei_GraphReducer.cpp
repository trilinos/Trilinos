/*--------------------------------------------------------------------*/
/*    Copyright 2007 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_GraphReducer.hpp>
#include <fei_TemplateUtils.hpp>
#include <fei_VectorSpace.hpp>

#undef fei_file
#define fei_file "fei_GraphReducer.cpp"
#include <fei_ErrMacros.hpp>

//----------------------------------------------------------------------------
fei::GraphReducer::GraphReducer(fei::SharedPtr<fei::Reducer> reducer,
                              fei::SharedPtr<fei::Graph> target)
  : reducer_(reducer),
    target_(target)
{
}

//----------------------------------------------------------------------------
fei::GraphReducer::~GraphReducer()
{
}

//----------------------------------------------------------------------------
int fei::GraphReducer::addIndices(int row, int len, const int* indices)
{
  reducer_->addGraphIndices(1, &row, len, indices, *target_);
  return(0);
}

//----------------------------------------------------------------------------
int fei::GraphReducer::addSymmetricIndices(int numIndices, int* indices,
                                         bool diagonal)
{
  reducer_->addSymmetricGraphIndices(numIndices, indices, diagonal, *target_);
  return(0);
}

//----------------------------------------------------------------------------
int fei::GraphReducer::writeLocalGraph(FEI_OSTREAM& os, bool debug,
				    bool prefixLinesWithPoundSign)
{
  return(target_->writeLocalGraph(os, debug, prefixLinesWithPoundSign));
}

//----------------------------------------------------------------------------
int fei::GraphReducer::writeRemoteGraph(FEI_OSTREAM& os)
{
  return(target_->writeRemoteGraph(os));
}

//----------------------------------------------------------------------------
int fei::GraphReducer::gatherFromOverlap()
{
  reducer_->assembleReducedGraph(target_.get(), false);
  return(target_->gatherFromOverlap());
}

