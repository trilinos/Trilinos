/*Paul
27-May-2002 General cleanup. Checked for newNamingConvention (nothing changed).
06-August-2002 Changed to images.
04-Feb-2003 Moved into Tpetra namespace & renamed
*/

#ifndef _TPETRA_COMBINEMODE_HPP_
#define _TPETRA_COMBINEMODE_HPP_
/*! \file Tpetra_CombineMode.hpp 
    \brief Tpetra::CombineMode enumerable type
*/

namespace Tpetra {

/*! \enum Tpetra::CombineMode 
    If set to Add, components on the receiving image will be added
    together.    If set to Zero, off-image components will be ignored.
    If set to Replace, off-image components will replace existing
    components on the receiving image.
    If set to Average, off-image components will be averaged with
    existing components on the receiving image.
*/

	enum CombineMode {Add,     /*!< Components on the receiving image
															 will be added together. */
										Zero,    /*!< Off-image components will be
															 ignored. */
										Insert,  /*!< Off-image components will
															 be inserted into locations on
															 receiving image. */
										Replace, /*!< Off-image components will
															 replace existing components on the 
															 receiving image. */
										Average  /*!< Off-image components will be
															 averaged with existing components 
															 on the receiving image. */
	};
	
} // namespace Tpetra
#endif /* _TPETRA_COMBINEMODE_HPP_ */
