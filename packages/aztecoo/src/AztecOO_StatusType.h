
/* Copyright (2001) Sandia Corportation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 *
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef AZTECOO_STATUSTYPE_H
#define AZTECOO_STATUSTYPE_H
/*! \file AztecOO_StatusType.h 
    \brief AztecOO StatusType: Used to return convergence status information for AztecOO_StatusTest objects.
 */

/*! \enum AztecOO_StatusType 
    When the CheckStatus and GetStatus methods of AztecOO_StatusTest objects are called a variable
    of type AztecOO_StatusType is returned.
*/



enum AztecOO_StatusType { Unchecked = 2,   /*!< Initial state of status */
			  Unconverged = 1, /*!< Convergence is not reached. */
			  Converged = 0,   /*!< Convergence is reached. */
			  Failed = -1,      /*!< Some failure occured.  Should stop */
			  NaN = -2         /*!< Result from test contains a NaN value.  Should stop */
			  
};
#endif /* AZTECOO_STATUSTYPE_H */
