/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_ostream_ops.hpp>

#include <fei_Vector.hpp>
#include <fei_Matrix.hpp>
#include <fei_SSMat.hpp>
#include <fei_SSVec.hpp>


FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::Vector& vec)
{
  vec.writeToStream(os);
  return(os);
}

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::Matrix& mat)
{
  mat.writeToStream(os);
  return(os);
}

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, SSVec& vec)
{
  vec.writeToStream(os);
  return(os);
}

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, SSMat& mat)
{
  mat.writeToStream(os);
  return(os);
}
