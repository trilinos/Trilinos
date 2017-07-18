#ifndef XPETRA_SPLITTINGDRIVER_DECL_HPP
#define XPETRA_SPLITTINGDRIVER_DECL_HPP

#include "Xpetra_Map.hpp"
#include <string>
#include <fstream>
#include <vector>

namespace Xpetra{

template<class GlobalOrdinal>
bool compareRegions(const std::tuple<GlobalOrdinal, GlobalOrdinal> &, const std::tuple<GlobalOrdinal, GlobalOrdinal> &);

template<class GlobalOrdinal>
bool compareNodes(const std::tuple<GlobalOrdinal, GlobalOrdinal> &, const std::tuple<GlobalOrdinal, GlobalOrdinal> &);

//Definition of the predicate for the node_ structude
template<class GlobalOrdinal>
class checkerNode { 
 
	public:  

		//Constructor
		checkerNode( GlobalOrdinal region_index){region_index_ = region_index;};

		//Unary Operator
  		bool operator()(const std::tuple<GlobalOrdinal, GlobalOrdinal> &node)  
  		{ return (std::get<1>(node) == region_index_); }  

	private:

		GlobalOrdinal region_index_;

};

 
//Definition of the predicate for the nodesToRegion_ structude
template<class GlobalOrdinal>
class checkerNodesToRegion { 
 
	public:  

		//Constructor
		checkerNodesToRegion( GlobalOrdinal node_index){node_index_ = node_index;};

		//Unary Operator
  		bool operator()(const std::tuple<GlobalOrdinal, Teuchos::Array<GlobalOrdinal> > &node)  
  		{ return (std::get<0>(node) == node_index_); }  

	private:

		GlobalOrdinal node_index_;

}; 

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	class Splitting_MapsInfo{
	public:
		Teuchos::Array<Teuchos::Array< std::tuple<GlobalOrdinal,GlobalOrdinal> > > regionToAll_;//used as a map for a RegionToAll node index
		Teuchos::Array<GlobalOrdinal> global_map_; //used as RowMap for global matrices
		Teuchos::Array<Teuchos::Array<GlobalOrdinal> > local_maps_; //used as RowMap for local matrices
	};

	template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  	class SplittingDriver{

	public:

	//! @name Constructor/Destructor Methods
	//@{

		//! Constructor specifying the file name containing regional information.
		SplittingDriver (){num_regional_nodes_.clear();};
		SplittingDriver (const std::string &, Teuchos::RCP< const Teuchos::Comm<int> >);

	//}
	//! @Interface methods
	//@{
		GlobalOrdinal GetNumGlobalElements(){return num_total_nodes_;};
		GlobalOrdinal GetNumTotalRegions(){return num_total_regions_;};
		Teuchos::Array<GlobalOrdinal> GetGlobalRowMap(){return maps_.global_map_;};
		Teuchos::Array<GlobalOrdinal> GetLocalRowMap(GlobalOrdinal region_index);
		Teuchos::Array<Teuchos::Array<GlobalOrdinal> > GetLocalRowMaps(){return maps_.local_maps_;};
	//}
	//! @Printout methods
		void printView();
		void printNodesToRegion();
		void printInactive();
	//}

		void CreateRowMaps();

		Teuchos::Array<GlobalOrdinal> num_regional_nodes_;
		
	private:
		Teuchos::RCP< const Teuchos::Comm<int> > comm_;
		bool nodes_sorted_by_regions_ = false;

		//Global information
		GlobalOrdinal num_total_nodes_ = 0;
		GlobalOrdinal num_total_regions_ = 0;	
		Teuchos::Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > nodes_;//basic structure that imports the information from the input file
		Teuchos::Array<GlobalOrdinal> regions_per_proc_;//if num_proc > num_regions, then it says how many regions are owned by a single process, empty otherwise
		Teuchos::Array<std::tuple<int, Teuchos::Array<GlobalOrdinal> > > procs_per_region_; //lists of processes instantiated for each region
		Teuchos::Array<std::tuple<int, Teuchos::Array<GlobalOrdinal> > > nodesToRegion_; //for each node it lists the regions it belongs to

		//Maps used for global and local operators
		Splitting_MapsInfo<Scalar, LocalOrdinal, GlobalOrdinal, Node> maps_;

		//Interface routines
		void ReadFileInfo(const std::string &);
		void ComputeProcRegions();
		void NodesToRegion();

  	}; //class SplittingDriver


	template<class GlobalOrdinal>
	bool compareRegions(const std::tuple<GlobalOrdinal, GlobalOrdinal> &lhs, const std::tuple<GlobalOrdinal, GlobalOrdinal> &rhs)
	{
		//First we prioritize the sorting according to the region label
		//If the region is the same, then the sorting looks at the global node index
		if( std::get<1>(lhs) < std::get<1>(rhs) )
			return true;
		else if( std::get<1>(lhs) == std::get<1>(rhs) )
			return std::get<0>(lhs) < std::get<0>(rhs); 
		else
			return false;
	}

	template<class GlobalOrdinal>
	bool compareNodes(const std::tuple<GlobalOrdinal, GlobalOrdinal> &lhs, const std::tuple<GlobalOrdinal, GlobalOrdinal> &rhs)
	{
		return std::get<0>(lhs) < std::get<0>(rhs); 
	}

} //namespace Xpetra

#endif
