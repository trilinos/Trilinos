// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Common.H"
#include "MaterialProps.H"

MaterialPropFactory &
MaterialPropFactory::factory()
{
  static MaterialPropFactory * self = new MaterialPropFactory();
  return *self;
}

MaterialPropBase *
MaterialPropFactory::add_model( MaterialPropBase * model )
{
  property_models[model->get_property_type()] = Teuchos::rcp(model);
  return model;
}

const MaterialPropBase &
MaterialPropFactory::get_model(  PelletTransport::ACTIVE_REGIONS region )
{
  if( property_models.end() == property_models.find(region) )
    throw "No valid property model.";

  return *(property_models[region]);
}

MaterialProp_He::MaterialProp_He() :
  MaterialPropBase()
{
}

PelletTransport::ACTIVE_REGIONS
MaterialProp_He::get_property_type()
{
  return PelletTransport::HELIUM;
}

void
MaterialProp_He::computeProps( double T, double x, MaterialPropBase::PropData & props )
{
  props.density   = 0.0818 - 8.275e-5*(T - 600.0);
  props.k_thermal = 0.0468 + 3.81e-4*T - 6.79e-8*T*T;
  props.Cp        = 5190.0;
  props.Qstar     = 0.0;
  props.Qdot      = 0.0;
  props.thermoF   = 1.0;
  props.D_diff    = 0.0;

  return;
}

MaterialProp_UO2::MaterialProp_UO2() :
  MaterialPropBase()
{
}

PelletTransport::ACTIVE_REGIONS
MaterialProp_UO2::get_property_type()
{
  return PelletTransport::UO2;
}

void
MaterialProp_UO2::computeProps( double T, double x, MaterialPropBase::PropData & props )
{
  double a, b, c, d;
  if( 923.0 > T )
  {
    a =  0.997      ;
    b =  9.082e-6   ;
    c = -2.705e-10 ;
    d =  4.391e-13  ;
  }
  else
  {
    a =  0.997      ;
    b =  1.179e-5   ;
    c = -2.429e-9   ;
    d =  1.219e-12  ;
  }
  double term     = a + b*T + c*T*T + d*T*T*T;
  props.density   = 10960.0 / pow(term, 3);

  double lambda0  = 1.0 / (3.24e-2 + 2.51e-4*T);
  double theta    = 3.67*exp(-4.73e-4*T)*sqrt(2.0*x*lambda0);
  props.k_thermal = lambda0*atan(theta)/theta + 5.95e-11*T*T*T;
  //props.k_thermal = lambda0* 1.0 + 5.95e-11*T*T*T;

  props.Cp        = 264256.0 + 47.0*T;
  props.Qstar     = -1380.8 - 134435.5*exp(-x/0.0261);
  props.Qdot      = qdot;
  props.thermoF   = (2.0 + x) / ( 2.0*(1.0-3.0*x)*(1.0-2.0*x) );
  props.D_diff    = exp( -9.386 - 4.26e3/T + 1.2e-3*T*x + 7.5e-4*T*log((2.0+x)/x) );
  //props.D_diff    = exp( log(10) * (-9.386 - 4.26e3/T + 1.2e-3*T*x + 7.5e-4*T*log((2.0+x)/x)/log(10) ) );

  return;
}

MaterialProp_Clad::MaterialProp_Clad() :
  MaterialPropBase()
{
}

PelletTransport::ACTIVE_REGIONS
MaterialProp_Clad::get_property_type()
{
  return PelletTransport::CLAD;
}

void
MaterialProp_Clad::computeProps( double T, double x, MaterialPropBase::PropData & props )
{
  props.density   = 7817.0;
  props.k_thermal = 10.98 + 0.014*T - 7.44e-6*T*T;
  props.Cp        = 420.0;
  props.Qstar     = 0.0;
  props.Qdot      = 0.0;
  props.thermoF   = 1.0;

  return;
}

