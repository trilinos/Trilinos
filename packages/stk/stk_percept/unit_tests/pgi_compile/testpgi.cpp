
/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <Intrepid_FieldContainer.hpp>
//#include <Intrepid_CellTools.hpp>
//#include "Intrepid_CellToolsDef.hpp"
#include "Shards_CellTopology.hpp"


#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"

// #ifdef __PGI
// #include "Intrepid_HGRAD_HEX_C1_FEMDef.hpp"
// #endif



int mymainpgi(int argc, char **argv);
int main(int argc, char **argv) { 
#if 1
  if( new Intrepid::Basis_HGRAD_HEX_C1_FEM<double, Intrepid::FieldContainer<double> >())
    {
    }
#endif
#if 0
  {
    using namespace Intrepid;
    using namespace shards;
    CellTopology cell_topo(shards::getCellTopologyData<shards::Hexahedron<8> >() );

    Teuchos::RCP< Basis< double, FieldContainer<double> > > HGRAD_Basis = getHGRAD_Basis<double >(cell_topo);
  }
#endif
  mymainpgi(argc, argv);
}
