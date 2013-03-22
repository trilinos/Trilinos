/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <test_utils/tester.hpp>

#undef fei_file
#define fei_file "tester.cpp"
#include <fei_ErrMacros.hpp>

tester::tester(MPI_Comm comm)
 : comm_(comm),
   numProcs_(1),
   localProc_(0),
   path_()
{
#ifndef FEI_SER
  MPI_Comm_rank(comm_, &localProc_);
  MPI_Comm_size(comm_, &numProcs_);
#endif
}

tester::~tester()
{
}

void tester::setPath(const std::string& path)
{
  path_ = path;
}

