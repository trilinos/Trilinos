
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

#ifndef _EPETRA_COMBINEMODE_H_
#define _EPETRA_COMBINEMODE_H_
/*! \file Epetra_CombineMode.h 
    \brief Epetra_Combine Mode enumerable type
 */

/*! \enum Epetra_CombineMode 
    If set to Add, components on the receiving processor will be added
    together.    If set to Zero, off-processor components will be ignored.
    If set to Insert, off-processor components will replace existing
    components on the receiving processor.
    If set to Average, off-processor components will be averaged with
    existing components on the receiving processor. (Recursive Binary Average)
    If set to AbsMax, magnitudes of off-processor components will be maxed
    with magnitudes of existing components of the receiving processor.
    { V = Supported by Epetra_Vector and Epetra_MultiVector,
      M = Supported by Epetra_CrsMatrix and Epetra_VbrMatrix }
*/

enum Epetra_CombineMode {Add,    /*!< Components on the receiving processor
                                     will be added together. (V,M) */
                        Zero,   /*!< Off-processor components will be
                                     ignored. (V,M) */
                        Insert, /*!< Off-processor components will
                                     be inserted into locations on
                                     receiving processor replacing existing values. (V,M) */
                        Average,/*!< Off-processor components will be
                                     averaged with existing components 
                                     on the receiving processor. (V) */
                        AbsMax  /*!< Magnitudes of Off-processor components will be
                                     maxed with magnitudes of existing components 
                                     on the receiving processor. (V) */
                        };

#endif // _EPETRA_COMBINEMODE_H_
