//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
