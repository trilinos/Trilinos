/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


//========================================================================================================================
//========================================================================================================================
//========================================================================================================================

#include <iostream>
#include <Intrepid_FieldContainer.hpp>
#include "Intrepid_DefaultCubatureFactory.hpp"

int main()
{
  using namespace Intrepid;
  using namespace shards;
  
  Intrepid::FieldContainer<double> aa(1,2);
  aa(0,1)=1.0;
  std::cout << "aa= \n" << aa << std::endl;

  CellTopology triangle_3( shards::getCellTopologyData<shards::Triangle<3> >() );

  DefaultCubatureFactory<double> cubFactory;                                              // create cubature factory
  unsigned cubDegree = 2;                                                                      // set cubature degree, e.g. 2
  Teuchos::RCP<Cubature<double> > myCub = cubFactory.create(triangle_3, cubDegree);         // create default cubature



}


