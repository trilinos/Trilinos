// 27-May-2002 General cleanup. Checked for newNamingConvention (nothing changed).

#ifndef _TPETRA_COMBINEMODE_H_
#define _TPETRA_COMBINEMODE_H_
/*! \file Tpetra_CombineMode.h 
    \brief Tpetra_CombineMode enumerable type
*/

/*! \enum Tpetra_CombineMode 
    If set to Add, components on the receiving processor will be added
    together.    If set to Zero, off-processor components will be ignored.
    If set to Insert, off-processor components will replace existing
    components on the receiving processor.
    If set to Average, off-processor components will be averaged with
    existing components on the receiving processor.
*/

enum Tpetra_CombineMode {Tpetra_Add,    /*!< Components on the receiving processor
                                     will be added together. */
                        Tpetra_Zero,   /*!< Off-processor components will be
                                     ignored. */
                        Tpetra_Insert, /*!< Off-processor components will
                                     be inserted into locations on
                                     receiving processor. */
                        Tpetra_Replace, /*!< Off-processor components will
                                     replace existing components on the 
                                     receiving processor. */
                        Tpetra_Average /*!< Off-processor components will be
                                     averaged with existing components 
                                     on the receiving processor. */
                        };

#endif /* _TPETRA_COMBINEMODE_H_ */
