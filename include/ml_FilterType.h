#ifndef ML_FILTERTYPE_H
#define ML_FILTERTYPE_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/* \file ml_FilterType.h
 *
 * \brief Enum for filtering.
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

namespace ML_Epetra {

/*! \enum FilterType
 *
 * \brief Defined the type of filter to be applied after each
 *  ExtractMyRowCopy().
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \date Last updated on 15-Mar-05.
 */

enum FilterType {
  ML_NO_FILTER,           /*< no filter is applied */
  ML_EQN_FILTER,          /*< decouples the equations */
  ML_TWO_BLOCKS_FILTER,   /*< decoupled the system in two blocks */
  ML_THREE_BLOCKS_FILTER, /*< decoupled the system in three blocks */
  ML_MASK_FILTER          /*< general approach */
};

} // namespace ML_Epetra
#endif
