/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#ifdef FEI_EXPLICIT_TEMPLATES

#include <feiArray.hpp>
#include <test_utils/DataReader.hpp>
#include <test_utils/tester.hpp>
#include <test_utils/feitester.hpp>
#include <test_utils/driverData.hpp>

template class feiArray<initElem*>;
template class feiArray<sumInElem*>;
template class feiArray<nodeBC*>;
template class feiArray<initCR*>;
template class feiArray<setIDLists*>;
template class feiArray<putBlockFieldNodeSolution*>;
template class feiArray<parameters*>;
template class feiArray<loadCR*>;
template class feiArray<sharedNodes*>;
template class feiArray<tester*>;
template class feitester<DataReader>;

#endif

