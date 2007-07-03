//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
