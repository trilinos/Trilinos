/*Paul
27-May-2002 General cleanup. Checked for newNamingConvention (nothing changed).
06-August-2002 Changed to images.
*/

#ifndef _TPETRA_COMBINEMODE_H_
#define _TPETRA_COMBINEMODE_H_
/*! \file Tpetra_CombineMode.h 
    \brief Tpetra_CombineMode enumerable type
*/

/*! \enum Tpetra_CombineMode 
    If set to Add, components on the receiving image will be added
    together.    If set to Zero, off-image components will be ignored.
    If set to Insert, off-image components will replace existing
    components on the receiving image.
    If set to Average, off-image components will be averaged with
    existing components on the receiving image.
*/

enum Tpetra_CombineMode {Tpetra_Add,    /*!< Components on the receiving image
                                     will be added together. */
                        Tpetra_Zero,   /*!< Off-image components will be
                                     ignored. */
                        Tpetra_Insert, /*!< Off-image components will
                                     be inserted into locations on
                                     receiving image. */
                        Tpetra_Replace, /*!< Off-image components will
                                     replace existing components on the 
                                     receiving image. */
                        Tpetra_Average /*!< Off-image components will be
                                     averaged with existing components 
                                     on the receiving image. */
                        };

#endif /* _TPETRA_COMBINEMODE_H_ */
