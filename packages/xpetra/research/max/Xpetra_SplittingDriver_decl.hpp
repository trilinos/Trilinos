#ifndef XPETRA_SPLITTINGDRIVER_DECL_HPP
#define XPETRA_SPLITTINGDRIVER_DECL_HPP

#include "Xpetra_Map.hpp"
#include <string>
#include <fstream>
#include <vector>

namespace Xpetra{

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	class Splitting_MapsInfo{
		std::vector<std::tuple<GlobalOrdinal,GlobalOrdinal> > node_;
		std::vector<GlobalOrdinal> global_map_; //used as RowMap for global matrices
		std::vector<std::vector<GlobalOrdinal> > local_maps_; //used as RowMap for local matrices
	};

	 template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class SplittingDriver{

	public:
    //! @name Constructor/Destructor Methods
    //@{

		//! Constructor specifying the file name containing regional information.
		SplittingDriver (const std::string, Teuchos::RCP< const Teuchos::Comm<int> >);

		void printView();

	private:
		std::vector<std::tuple<GlobalOrdinal, GlobalOrdinal> > node_;
		Xpetra::Splitting_MapsInfo<Scalar, LocalOrdinal, GlobalOrdinal, Node> maps_;

  }; //class SplittingDriver

} //namespace Xpetra

#endif
