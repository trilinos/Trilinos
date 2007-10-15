/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_MatrixGraph.hpp>
#include <fei_LinearSystem.hpp>
#include <snl_fei_LinearSystem_General.hpp>
#include <snl_fei_Utils.hpp>

//----------------------------------------------------------------------------
fei::LinearSystem::LinearSystem()
 : matrix_(),
   soln_(),
   rhs_()
{
}

//----------------------------------------------------------------------------
fei::LinearSystem::~LinearSystem()
{
}

//----------------------------------------------------------------------------
fei::SharedPtr<fei::LinearSystem>
fei::LinearSystem::Factory::createLinearSystem(fei::SharedPtr<fei::MatrixGraph>& matrixGraph)
{
  fei::SharedPtr<fei::LinearSystem>
    linsys(new snl_fei::LinearSystem_General(matrixGraph));

  return(linsys);
}

//----------------------------------------------------------------------------
void fei::LinearSystem::setMatrix(fei::SharedPtr<fei::Matrix>& matrix)
{
  matrix_ = matrix;
}

//----------------------------------------------------------------------------
int fei::LinearSystem::putAttribute(const char* name,
                                    void* attribute)
{
  snl_fei::storeNamedAttribute(name, attribute,
                               attributeNames_, attributes_);
  return(0);
}

//----------------------------------------------------------------------------
int fei::LinearSystem::getAttribute(const char* name,
                                    void*& attribute)
{
  attribute = snl_fei::retrieveNamedAttribute(name, attributeNames_, attributes_);
  return(attribute==NULL ? -1 : 0);
}

