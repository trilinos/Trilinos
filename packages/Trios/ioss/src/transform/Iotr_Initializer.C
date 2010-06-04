/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <transform/Iotr_Initializer.h>
#include <transform/Iotr_MinMax.h>
#include <transform/Iotr_Offset.h>
#include <transform/Iotr_Scale.h>
#include <transform/Iotr_Offset3D.h>
#include <transform/Iotr_Scale3D.h>
#include <transform/Iotr_Tensor.h>
#include <transform/Iotr_VectorMagnitude.h>

namespace Iotr {
  Initializer::Initializer()
  {
    MinMax_Factory::factory();
    Offset_Factory::factory();
    Scale_Factory::factory();
    Offset3D_Factory::factory();
    Scale3D_Factory::factory();
    Tensor_Factory::factory();
    VM_Factory::factory();
  }
}
