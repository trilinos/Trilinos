// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// 1D Linear basis function for finite element method

#ifndef _NOX_TEST_EPETRA_PELLETTRANSPORT_MATERIALPROPS_H
#define _NOX_TEST_EPETRA_PELLETTRANSPORT_MATERIALPROPS_H

#include "Teuchos_RCP.hpp"
#include <string>

#include "PelletTransport.H"

class MaterialPropBase {

 public:

 class PropData
 {
  public:

   double density   ;
   double k_thermal ;
   double Cp        ;
   double Qstar     ;
   double Qdot      ;
   double thermoF   ;
   double D_diff    ;
 };

  // Constructor
  MaterialPropBase() { }

  // Destructor
  virtual ~MaterialPropBase() { }

  // Calculates the values of u and x at the specified gauss point
  virtual void computeProps( double T, double x, PropData & props ) = 0;

  virtual PelletTransport::ACTIVE_REGIONS get_property_type()
  { return PelletTransport::NONE; }

};

class MaterialPropFactory {

 public:

  static MaterialPropFactory & factory();

  static void computeProps( PelletTransport::ACTIVE_REGIONS propType, double T, double x, MaterialPropBase::PropData & props )
  { return factory().property_models[propType]->computeProps( T, x, props ); }

  MaterialPropBase *  add_model( MaterialPropBase * );
  const MaterialPropBase &  get_model( PelletTransport::ACTIVE_REGIONS );

 private:

   MaterialPropFactory() { };
   ~MaterialPropFactory();
   MaterialPropFactory( const MaterialPropFactory & );

   map< PelletTransport::ACTIVE_REGIONS, Teuchos::RCP<MaterialPropBase> > property_models;

};

class MaterialProp_He : public MaterialPropBase
{

  public:

  // Constructor
  MaterialProp_He();

  // Destructor
  virtual ~MaterialProp_He() { }

  // Calculates the values of u and x at the specified gauss point
  void computeProps( double T, double x, MaterialPropBase::PropData & props );

  PelletTransport::ACTIVE_REGIONS get_property_type();

};

class MaterialProp_UO2 : public MaterialPropBase
{

  public:

  // Constructor
  MaterialProp_UO2();

  // Destructor
  virtual ~MaterialProp_UO2() { }

  // Calculates the values of u and x at the specified gauss point
  void computeProps( double T, double x, MaterialPropBase::PropData & props );

  PelletTransport::ACTIVE_REGIONS get_property_type();

  void set_qdot( double qdot_ )
    { qdot = qdot_; }

  void set_xB( double xB_ )
    { xB = xB_; }

  const double & get_xB( ) const
    { return xB; }

 private:

  double qdot;
  double xB;

};

class MaterialProp_Clad : public MaterialPropBase
{

  public:

  // Constructor
  MaterialProp_Clad();

  // Destructor
  virtual ~MaterialProp_Clad() { }

  // Calculates the values of u and x at the specified gauss point
  void computeProps( double T, double x, MaterialPropBase::PropData & props );

  PelletTransport::ACTIVE_REGIONS get_property_type();

};

#endif
