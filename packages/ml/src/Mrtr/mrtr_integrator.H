/*
#@HEADER
# ************************************************************************
#
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifdef TRILINOS_PACKAGE

#ifndef MOERTEL_INTEGRATOR_H
#define MOERTEL_INTEGRATOR_H

#include <ctime>
#include <iostream>
#include <iomanip>
#include <vector>

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_RefCountPtr.hpp"

#include "mrtr_overlap.H"

using namespace std;
using namespace Teuchos;

// ----------   User Defined Includes   ----------

namespace MOERTEL
{

// forward declarations
class MOERTEL::Interface;
class MOERTEL::Segment;
class MOERTEL::Node;

class Integrator 
{
public:
  
  // ctor
  explicit Integrator(int ngp, bool oneD, int outlevel);
  
  // dtor
  virtual ~Integrator();
  
  // return the level of output
  int OutLevel() { return outputlevel_; }
  
  // get the order of integration (== # gaussian points)
  inline int Ngp() { return ngp_; }
  
  // return the coordinate of a given gaussian point gp
  inline double Coordinate(int gp) { return coords_[gp]; }
  inline double* Coordinate(int* gp) { return &coords_[*gp*2]; }
  
  // return the weight of a given gaussian point gp
  inline double Weight(int gp) { return weights_[gp]; }


  // integrate a 1D overlap between a slave and a master segment
  Epetra_SerialDenseMatrix* Integrate(MOERTEL::Segment& sseg, double sxia, double sxib,
                                      MOERTEL::Segment& mseg, double mxia, double mxib);
  // assemble the result -Mdense from the integration above into M
  bool Assemble(MOERTEL::Interface& inter, MOERTEL::Segment& sseg, MOERTEL::Segment& mseg, 
                Epetra_CrsMatrix& M, Epetra_SerialDenseMatrix& Mdense);


  // integrate a 1D slave segment
  Epetra_SerialDenseMatrix* Integrate(MOERTEL::Segment& sseg, double sxia, double sxib);
  // assemble the result Ddense from the integration above into D
  bool Assemble(MOERTEL::Interface& inter, MOERTEL::Segment& sseg,  
                Epetra_CrsMatrix& D, Epetra_SerialDenseMatrix& Ddense);


  // integrate modification of M in 1D
  Epetra_SerialDenseMatrix* Integrate_2D_Mmod(MOERTEL::Segment& sseg, double sxia, double sxib,
                                              MOERTEL::Segment& mseg, double mxia, double mxib);
  // assemble the result -Mdense from the integration above into M
  bool Assemble_2D_Mod(MOERTEL::Interface& inter, MOERTEL::Segment& sseg, MOERTEL::Segment& mseg, 
                Epetra_CrsMatrix& M, Epetra_SerialDenseMatrix& Mmod);


  // integrate a 2D triangle overlap segment, master part M and slave part D
  bool Integrate(RefCountPtr<MOERTEL::Segment> actseg,
                 MOERTEL::Segment& sseg, 
                 MOERTEL::Segment& mseg,
                 Epetra_SerialDenseMatrix** Ddense, 
                 Epetra_SerialDenseMatrix** Mdense, 
                 MOERTEL::Overlap& overlap, double eps);

  // Assemble D into nodes
  bool Assemble(MOERTEL::Interface& inter,MOERTEL::Segment& sseg,Epetra_SerialDenseMatrix& Ddense);

  // Assemble M into nodes
  bool Assemble(MOERTEL::Interface& inter,MOERTEL::Segment& sseg,MOERTEL::Segment& mseg,
                Epetra_SerialDenseMatrix& Mdense);

private:  
  // don't want = operator
  Integrator operator = (const Integrator& old);
  // don't want copy-ctor
  Integrator(MOERTEL::Integrator& old);

private:

  bool           oneD_;    // dimension of the integration, true if 1-dimensional
  int            ngp_;     // order of the integration
  int            outputlevel_; // output level (0-10)
  vector<double> coords_;  // vector holding gaussian point coordinates
  vector<double> weights_; // vector holding gaussian point weights


};

} // namespace MOERTEL

#endif // MOERTEL_INTEGRATOR_H
#endif // TRILINOS_PACKAGE
