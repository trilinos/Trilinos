#ifndef _IFPACK_SCALINGTYPE_H_
#define _IFPACK_SCALINGTYPE_H_
/*! \file Ifpack_ScalingType.h
    \brief Ifpack_ScalingType enumerable type
 */

//! Ifpack scaling type selector.
/*! Selects the type of scaling used (if any) for Ifpack preconditioners.
*/
enum Ifpack_ScalingType {None, LeftDiagonal, RightDiagonal, 
                     SymmetricDiagonal, RowSum, ColSum, 
		           RowAndColSum};

#endif /* _IFPACK_SCALINGTYPE_H_ */
