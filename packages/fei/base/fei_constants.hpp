#ifndef _fei_constants_hpp_
#define _fei_constants_hpp_

/*
   In this file we define some constants to use as parameters to
   some fei functions.
   These constants are primarily used as 'fieldType' arguments to
   fei::VectorSpace::defineFields and FEI::initFields.
   If defining a vector-field for displacement, use DISPLACEMENT. If
   defining separate scalar fields for the components of displacement,
   then use DISPLACEMENT_X, etc.
*/

namespace fei {

const int DISPLACEMENT     =  0;
const int DISPLACEMENT_X   =  0;
const int DISPLACEMENT_Y   =  1;
const int DISPLACEMENT_Z   =  2;
const int ROTATION         =  3;
const int ROTATION_X       =  3;
const int ROTATION_Y       =  4;
const int ROTATION_Z       =  5;
const int VELOCITY         =  6;
const int VELOCITY_X       =  6;
const int VELOCITY_Y       =  7;
const int VELOCITY_Z       =  8;
const int PRESSURE         =  9;
const int PRESSURE_X       =  9;
const int PRESSURE_Y       = 10;
const int PRESSURE_Z       = 11;
const int TEMPERATURE      = 12;

const int UNKNOWN          = 20;

}//namespace fei

#endif

