// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include <MFLOCA.H>

extern "C" {
#include <stdio.h>
#include <MFNRegion.h>
#include <stdlib.h>

int LOCATest(MFNVector, void *);

MFNRegion MFNRegionCreateLOCA(LOCAData* data)
 {
  static char RoutineName[]={"MFNRegionCreateLOCA"};
  MFNRegion loca;

  loca=MFNRegionCreateBaseClass("LOCA");
  MFNRegionSetTest(loca,LOCATest);
  MFNRegionSetData(loca,(void *)data);

  return(loca);
 }

int LOCATest(MFNVector u, void *d)
{
   
   LMCEV* v = (LMCEV *)MFNVectorGetData(u);
   LOCAData* data = (LOCAData*) d;

   list<ParamData>::iterator it = data->paramData.begin();
   for (unsigned int i=0; i<data->paramData.size(); i++) {

     if (v->getScalar(i) < it->minValue)
       return 0;

     if (v->getScalar(i) > it->maxValue)
       return 0;

     ++it;

   }

   if (v->getXVec().norm(NOX::Abstract::Vector::MaxNorm) > data->solutionMax)
     return 0;
   
   return 1;
}


}
