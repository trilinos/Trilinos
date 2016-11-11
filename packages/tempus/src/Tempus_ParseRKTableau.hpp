#ifndef __Drekar_ParseRKTableau_hpp__
#define __Drekar_ParseRKTableau_hpp__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Tempus_RKButcherTableau.hpp"

namespace Tempus {

/** Parse an RK Butcher Tableau. The expected format
  * for the input ParameterList is:

    \verbatim
      <Parameter name="A" type="string" value="# # # ;
                                               # # # ;
                                               # # #"/>
      <Parameter name="b" type="string" value="# # #"/>
      <Parameter name="c" type="string" value="# # #"/>
    \endverbatim

  * Where "c" is the time step nodes and "A" and "b" are
  * the other components of the Butcher tableau.
  */
Teuchos::RCP<RKButcherTableau<double> >
parseRKTableau(const Teuchos::ParameterList & pl);

}

#endif
