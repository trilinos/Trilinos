#ifndef XPETRA_SPLITTINGDRIVER_DEF_HPP
#define XPETRA_SPLITTINGDRIVER_DEF_HPP

#include "Xpetra_SplittingDriver_decl.hpp"
#include <algorithm>

namespace Xpetra{

	template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	SplittingDriver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SplittingDriver(const std::string &file_name, Teuchos::RCP< const Teuchos::Comm<int> > comm): comm_(comm)
	{
		ReadFileInfo(file_name);

		//Nodes are shuffled so that region are sorted in ascending labeling order
		std::sort(nodes_.begin(), nodes_.end(), compareNodesRegions<GlobalOrdinal>);
		nodes_sorted_by_regions_ = true;

		NodesToRegion();
		ComputeProcRegions();
		CreateRowMaps();
	}


	template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	void SplittingDriver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReadFileInfo(const std::string &file_name)
	{
		std::ifstream input_file_(file_name, std::ifstream::in);
		std::string  line;
		TEUCHOS_TEST_FOR_EXCEPTION( !input_file_.good(), Exceptions::RuntimeError, "Can not read \"" << file_name << "\"");

		GlobalOrdinal line_index = 0;
	
		while ( std::getline (input_file_,line) )
		{
			std::istringstream is( line );
			GlobalOrdinal number;
 			Teuchos::Array<GlobalOrdinal> node;
			std::tuple<GlobalOrdinal, GlobalOrdinal> node_region;
			Teuchos::Array<GlobalOrdinal> global_info;

			node.clear();
			global_info.clear();

			if( 1==line_index )
			{	
				while(is>>number)
					global_info.push_back(number);

				TEUCHOS_TEST_FOR_EXCEPTION( global_info.size()!=2, Exceptions::RuntimeError, "The global information must be a couple of integers: nTotal of nodes + nTotal of regions \n");
				num_total_nodes_ = global_info[0];
				num_total_regions_ = global_info[1];
			}
			else if( line_index>2 )
			{
				while(is>>number)
				{
					node.push_back(number);
				}
				TEUCHOS_TEST_FOR_EXCEPTION( node.size()!=2, Exceptions::RuntimeError, "The node information must be a couple of integers: Node index + Region idnex \n");
				node_region = std::make_tuple(node[0], node[1]);
				nodes_.push_back(node_region);
				node.clear();
			}
			line_index++;
		}
		input_file_.close();
	}


	template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	void SplittingDriver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ComputeProcRegions()
	{
		int tot_num_proc = comm_->getSize();
		int myPID = comm_->getRank();

		regions_per_proc_.clear();
		procs_per_region_.clear();	

		if( tot_num_proc < num_total_regions_ )
		{
			int min_nregions_proc = std::floor( static_cast<double>(num_total_regions_)/static_cast<double>(tot_num_proc) );
			int num_leftover_regions = num_total_regions_ % tot_num_proc;

			for( int i=1; i<=min_nregions_proc; ++i )
				regions_per_proc_.push_back( myPID*min_nregions_proc+i );		
			
			if( num_leftover_regions<=myPID+1 )
				regions_per_proc_.push_back( min_nregions_proc*tot_num_proc + (myPID+1) );

			for( int procID = 0; procID<comm_->getSize(); ++procID )
			{
				Teuchos::Array<GlobalOrdinal> proc;
				proc.clear();
				proc.push_back(procID);
				for( int i = 1; i<=min_nregions_proc; ++i )	
				{
					GlobalOrdinal region_index = (procID)*min_nregions_proc + i;
					std::tuple<GlobalOrdinal, Teuchos::Array<GlobalOrdinal> > tuple_aux = std::make_tuple(region_index, proc);
					procs_per_region_.push_back( tuple_aux );
				}

				if( num_leftover_regions<=procID+1 )			
				{
					GlobalOrdinal region_index = min_nregions_proc*tot_num_proc + (procID+1);
					std::tuple<GlobalOrdinal, Teuchos::Array<GlobalOrdinal> > tuple_aux = std::make_tuple(region_index, proc);
					procs_per_region_.push_back( tuple_aux );	
				}
			}
			TEUCHOS_TEST_FOR_EXCEPTION( !( procs_per_region_.size()==num_total_regions_ ), Exceptions::RuntimeError, "Number of regions detected does not match with the initially declared one \n");
		}
		else if( tot_num_proc == num_total_regions_ )
		{
			regions_per_proc_.push_back( myPID+1 );

			for( int i = 0; i<num_total_regions_; ++i )
			{
				GlobalOrdinal region_index = i + 1;
				Teuchos::Array<GlobalOrdinal> proc;
				proc.clear();
				proc.push_back(i);
				std::tuple<GlobalOrdinal, Teuchos::Array<GlobalOrdinal> > tuple_aux = std::make_tuple(region_index, proc);
				procs_per_region_.push_back( tuple_aux );	
			}
		}
		else if( tot_num_proc > num_total_regions_ )
		{
			int num_procs_region = std::ceil( static_cast<double>(tot_num_proc)/static_cast<double>(num_total_regions_) );	
			int num_regions_extra_proc = tot_num_proc % num_total_regions_;
			int proc_count = 0;
			std::tuple<int, Teuchos::Array<GlobalOrdinal> > region_tuple;

			for( int i = 1; i<=num_total_regions_; ++i )			
			{
				Teuchos::Array<GlobalOrdinal> procs;
				procs.clear();
				if( i<=num_regions_extra_proc || num_regions_extra_proc==0 )				
					for( int j=1; j<=num_procs_region; ++j )
					{
						procs.push_back(proc_count);
						proc_count++;
					}
				else
					for( int j=1; j<=num_procs_region-1; ++j )
					{
						procs.push_back(proc_count);
						proc_count++;
					}
				std::sort(procs.begin(), procs.end());
				region_tuple = std::make_tuple(i, procs);
				procs_per_region_.push_back(region_tuple);
			}
			regions_per_proc_.clear();
		}
	}


	template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	void SplittingDriver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NodesToRegion()
	{
		nodesToRegion_.clear();
		Teuchos::Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > nodes_reordered;
		nodes_reordered = nodes_;
		std::sort(nodes_reordered.begin(), nodes_reordered.end(), compareNodes<GlobalOrdinal>);

		typename Teuchos::Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator node_iterator;
		node_iterator = nodes_reordered.begin();
	
		while( node_iterator != nodes_reordered.end() )
		{
			GlobalOrdinal current_node = std::get<0>( *(node_iterator) );
			Teuchos::Array<GlobalOrdinal> regions;
			regions.clear();
			regions.push_back( std::get<1>(*node_iterator) );
			
			typename Teuchos::Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator next_node_iterator = node_iterator + 1;
			
			while( next_node_iterator != nodes_reordered.end() )
			{
				GlobalOrdinal next_node = std::get<0>( *(next_node_iterator) );
				if( current_node == next_node )
				{
					regions.push_back( std::get<1>( *(next_node_iterator) ) );
					next_node_iterator++;
				}
				else
				{
					node_iterator = next_node_iterator;
					break;
				}
			}
			std::tuple<GlobalOrdinal, Teuchos::Array<GlobalOrdinal> > new_tuple;
			std::sort(regions.begin(), regions.end());
			new_tuple = std::make_tuple(current_node, regions);
			nodesToRegion_.push_back(new_tuple);
			
			if( next_node_iterator == nodes_reordered.end() )
				break;
		}
		TEUCHOS_TEST_FOR_EXCEPTION( !( nodesToRegion_.size()==num_total_nodes_ ), Exceptions::RuntimeError, "Number of nodes detected does not match with the initially declared one \n");
	}

	template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	void SplittingDriver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateRowMaps()
	{
		TEUCHOS_TEST_FOR_EXCEPTION( !( procs_per_region_.empty() || regions_per_proc_.empty() ), Exceptions::RuntimeError, "Missing information about region partitioning across processors instantiated \n");
		TEUCHOS_TEST_FOR_EXCEPTION( ( procs_per_region_.empty() && regions_per_proc_.empty() ), Exceptions::RuntimeError, "Information about region partitioning across processors is not consistent: incorrect values for number of processors or number of regions \n");
		Teuchos::Array<GlobalOrdinal> elements;
		Teuchos::Array<GlobalOrdinal> regional_elements;
		Teuchos::Array<Teuchos::Array<GlobalOrdinal> > elements_per_region;
		int myPID = comm_->getRank(); 

		elements.clear();
		regional_elements.clear();
		elements_per_region.resize(num_total_regions_);

		TEUCHOS_TEST_FOR_EXCEPTION( !nodes_sorted_by_regions_, Exceptions::RuntimeError, "Nodes are not sorted by regions in ascending order \n");
		TEUCHOS_TEST_FOR_EXCEPTION( num_total_nodes_>nodes_.size(), Exceptions::RuntimeError, "Number of nodes declared in input file does not match with the effective number of nodes provided\n");
		TEUCHOS_TEST_FOR_EXCEPTION( num_total_regions_!=std::get<1>( *(nodes_.end()-1) ), Exceptions::RuntimeError, "Number of regions declared in input file does not match with the effective number of regions provided\n");

		if( !( regions_per_proc_.empty() ) )
		{
			typename Teuchos::Array<GlobalOrdinal>::iterator iter_array; 
			for( iter_array=regions_per_proc_.begin(); iter_array!=regions_per_proc_.end(); ++iter_array )
			{
				regional_elements.clear();
				checkerNode<GlobalOrdinal> unaryPredicate(*iter_array);
				typename Teuchos::Array< std::tuple<GlobalOrdinal,GlobalOrdinal> >::iterator nodes_iterator1;
				typename Teuchos::Array< std::tuple<GlobalOrdinal,GlobalOrdinal> >::iterator nodes_iterator2;
				nodes_iterator1 = std::find_if<typename Teuchos::Array< std::tuple<GlobalOrdinal,GlobalOrdinal> >::iterator, checkerNode<GlobalOrdinal> >(nodes_.begin(), nodes_.end(), unaryPredicate);
				nodes_iterator2 = std::find_if_not<typename Teuchos::Array< std::tuple<GlobalOrdinal,GlobalOrdinal> >::iterator, checkerNode<GlobalOrdinal> >(nodes_iterator1, nodes_.end(), unaryPredicate);
						
				typename Teuchos::Array< std::tuple<GlobalOrdinal,GlobalOrdinal> >::iterator nodes_iterator_aux;
				GlobalOrdinal local_node_label = 0;
				for( nodes_iterator_aux=nodes_iterator1; nodes_iterator_aux!=nodes_iterator2; ++nodes_iterator_aux )
				{
					GlobalOrdinal node = std::get<0>(*nodes_iterator_aux);
					checkerNodesToRegion<GlobalOrdinal> unaryPredicateNode(node);
					typename Teuchos::Array< std::tuple<GlobalOrdinal, Teuchos::Array<GlobalOrdinal> > >::iterator nodes_to_region_iterator;
					nodes_to_region_iterator = std::find_if<typename Teuchos::Array< std::tuple< GlobalOrdinal,Teuchos::Array<GlobalOrdinal> > >::iterator, checkerNodesToRegion<GlobalOrdinal> >(nodesToRegion_.begin(), nodesToRegion_.end(), unaryPredicateNode);	
					Teuchos::Array<GlobalOrdinal> nodal_regions =  std::get<1>(*nodes_to_region_iterator);

					if( *iter_array==nodal_regions[0] )
						elements.push_back( node );

					regional_elements.push_back(local_node_label);
					local_node_label++;
				}
				elements_per_region[*iter_array-1] = regional_elements;
			}
		}
		else
		{
			bool region_found = false;			
			GlobalOrdinal myRegion = -1;
			TEUCHOS_TEST_FOR_EXCEPTION( !( procs_per_region_.size() == num_total_regions_ ), Exceptions::RuntimeError, "Number of total regions does not match with driver structures \n");

			Teuchos::Array<GlobalOrdinal> regional_procs;
			while( !region_found )
			{
				typename Teuchos::Array<GlobalOrdinal>::iterator iter_proc;
				for( GlobalOrdinal region_index=1; region_index<=procs_per_region_.size(); ++region_index )
				{
					regional_procs = std::get<1>(procs_per_region_[region_index-1]);
					iter_proc = std::find( regional_procs.begin(), regional_procs.end(), myPID );
					if( iter_proc!=regional_procs.end() )
					{
						myRegion = region_index;
						region_found = true;
					}
				}
			}

			TEUCHOS_TEST_FOR_EXCEPTION( ( myRegion == -1 || !region_found ), Exceptions::RuntimeError, ( "Region containing PROC ID: "+  std::to_string(myPID) + " NOT FOUND \n" ) );
			regional_procs = std::get<1>(procs_per_region_[myRegion-1]);

			checkerNode<GlobalOrdinal> unaryPredicate(myRegion);
			typename Teuchos::Array< std::tuple<GlobalOrdinal,GlobalOrdinal> >::iterator nodes_iterator1;
			typename Teuchos::Array< std::tuple<GlobalOrdinal,GlobalOrdinal> >::iterator nodes_iterator2;
			nodes_iterator1 = std::find_if<typename Teuchos::Array< std::tuple<GlobalOrdinal,GlobalOrdinal> >::iterator, checkerNode<GlobalOrdinal> >(nodes_.begin(), nodes_.end(), unaryPredicate);	
			nodes_iterator2 = std::find_if_not<typename Teuchos::Array< std::tuple<GlobalOrdinal,GlobalOrdinal> >::iterator, checkerNode<GlobalOrdinal> >(nodes_iterator1, nodes_.end(), unaryPredicate);

			int num_regional_nodes = nodes_iterator2 - nodes_iterator1;
			int num_regional_procs = regional_procs.size();

			if( num_regional_nodes < num_regional_procs )
			{
				Teuchos::Array<GlobalOrdinal> regional_procs_reduced;
				regional_procs_reduced.clear();
				for(int i = 0; i<num_regional_nodes; ++i)
					regional_procs_reduced.push_back( regional_procs[i] );

				typename Teuchos::Array<GlobalOrdinal>::iterator proc_iterator;
				proc_iterator = std::find<typename Teuchos::Array<GlobalOrdinal>::iterator, GlobalOrdinal>(regional_procs_reduced.begin(), regional_procs_reduced.end(), myPID);	

				if( proc_iterator!=regional_procs_reduced.end() )//This reasoning works because the PROC ID for each region has been previously sorted in ascending order
				{
					GlobalOrdinal node = std::get<0>( *( nodes_iterator1+(proc_iterator-regional_procs_reduced.begin()+1) ) );
					GlobalOrdinal local_node_label = proc_iterator-regional_procs_reduced.begin();
					checkerNodesToRegion<GlobalOrdinal> unaryPredicateNode(node);
					typename Teuchos::Array< std::tuple<GlobalOrdinal, Teuchos::Array<GlobalOrdinal> > >::iterator nodes_to_region_iterator;
					nodes_to_region_iterator = std::find_if<typename Teuchos::Array< std::tuple< GlobalOrdinal,Teuchos::Array<GlobalOrdinal> > >::iterator, checkerNodesToRegion<GlobalOrdinal> >(nodesToRegion_.begin(), nodesToRegion_.end(), unaryPredicateNode);	
					Teuchos::Array<GlobalOrdinal> nodal_regions =  std::get<1>(*nodes_to_region_iterator);
					if( myRegion == nodal_regions[0] )
						elements.push_back( node );

					regional_elements.push_back(local_node_label);
				}
			}
			else if( num_regional_nodes == num_regional_procs )		
			{
				typename Teuchos::Array<GlobalOrdinal>::iterator proc_iterator;
				proc_iterator = std::find<typename Teuchos::Array<GlobalOrdinal>::iterator, GlobalOrdinal>(regional_procs.begin(), regional_procs.end(), myPID);	

				if( proc_iterator!=regional_procs.end() )//This reasoning works because the PROC ID for each region has been previously sorted in ascending order
				{
					GlobalOrdinal node = std::get<0>( *( nodes_iterator1+(proc_iterator-regional_procs.begin()+1) ) );
					GlobalOrdinal local_node_label = proc_iterator-regional_procs.begin();
					checkerNodesToRegion<GlobalOrdinal> unaryPredicateNode(node);
					typename Teuchos::Array< std::tuple<GlobalOrdinal, Teuchos::Array<GlobalOrdinal> > >::iterator nodes_to_region_iterator;
					nodes_to_region_iterator = std::find_if<typename Teuchos::Array< std::tuple< GlobalOrdinal,Teuchos::Array<GlobalOrdinal> > >::iterator, checkerNodesToRegion<GlobalOrdinal> >(nodesToRegion_.begin(), nodesToRegion_.end(), unaryPredicateNode);	
					Teuchos::Array<GlobalOrdinal> nodal_regions =  std::get<1>(*nodes_to_region_iterator);
					if( myRegion == nodal_regions[0] )
						elements.push_back( node );

					regional_elements.push_back(local_node_label);
				}
			}
			else
			{
				typename Teuchos::Array<GlobalOrdinal>::iterator proc_iterator;
				proc_iterator = std::find<typename Teuchos::Array<GlobalOrdinal>::iterator, GlobalOrdinal>(regional_procs.begin(), regional_procs.end(), myPID);	

				int num_nodes_proc = std::ceil( static_cast<double>(num_regional_nodes)/static_cast<double>(num_regional_procs) );	
				int num_procs_extra_node = num_regional_nodes % num_regional_procs;

				if( proc_iterator-regional_procs.begin()+1 <= num_procs_extra_node || num_procs_extra_node == 0 )
				{
					int init_node = num_nodes_proc * ( proc_iterator-regional_procs.begin() );
					for( int i=0; i<num_nodes_proc; ++i )
					{
						GlobalOrdinal node = std::get<0>( *( nodes_iterator1 + init_node + i ) );
						GlobalOrdinal local_node_label = init_node + i;
						checkerNodesToRegion<GlobalOrdinal> unaryPredicateNode(node);
						typename Teuchos::Array< std::tuple<GlobalOrdinal, Teuchos::Array<GlobalOrdinal> > >::iterator nodes_to_region_iterator;
						nodes_to_region_iterator = std::find_if<typename Teuchos::Array< std::tuple< GlobalOrdinal,Teuchos::Array<GlobalOrdinal> > >::iterator, checkerNodesToRegion<GlobalOrdinal> >(nodesToRegion_.begin(), nodesToRegion_.end(), unaryPredicateNode);	
						Teuchos::Array<GlobalOrdinal> nodal_regions =  std::get<1>(*nodes_to_region_iterator);
						if( myRegion == nodal_regions[0] )
							elements.push_back( node );

						regional_elements.push_back(local_node_label);
					}
				}
				else
				{
					int init_node = num_nodes_proc * num_procs_extra_node + (proc_iterator - regional_procs.begin() - num_procs_extra_node) * (num_nodes_proc-1); 
					for( int i=0; i<num_nodes_proc-1; ++i )
					{
						GlobalOrdinal node = std::get<0>( *( nodes_iterator1 + init_node + i ) );
						GlobalOrdinal local_node_label = init_node + i;
						checkerNodesToRegion<GlobalOrdinal> unaryPredicateNode(node);
						typename Teuchos::Array< std::tuple<GlobalOrdinal, Teuchos::Array<GlobalOrdinal> > >::iterator nodes_to_region_iterator;
						nodes_to_region_iterator = std::find_if<typename Teuchos::Array< std::tuple< GlobalOrdinal,Teuchos::Array<GlobalOrdinal> > >::iterator, checkerNodesToRegion<GlobalOrdinal> >(nodesToRegion_.begin(), nodesToRegion_.end(), unaryPredicateNode);	
						Teuchos::Array<GlobalOrdinal> nodal_regions =  std::get<1>(*nodes_to_region_iterator);
						if( myRegion == nodal_regions[0] )
							elements.push_back( node );

						regional_elements.push_back(local_node_label);
					}
				}		
			}	
			elements_per_region[myRegion-1] = regional_elements;
		}

		for( typename Teuchos::Array<GlobalOrdinal>::iterator iter = elements.begin(); iter!=elements.end(); ++iter )
			*iter = *iter - 1;

		maps_.global_map_ = elements;
		maps_.local_maps_ = elements_per_region;
	}

	template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	Teuchos::Array<GlobalOrdinal> SplittingDriver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetLocalRowMap(GlobalOrdinal region_index)
	{
		TEUCHOS_TEST_FOR_EXCEPTION( region_index>num_total_regions_, Exceptions::RuntimeError, "Value of region index exceeds total number of regions stored \n");
		return maps_.local_maps_[region_index-1];
	}

	template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	void SplittingDriver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::printView()
	{
		if( 0==comm_->getRank() )
		{
			std::cout<<"Total number of mesh nodes: "<<num_total_nodes_<<std::endl;
			std::cout<<"Total number of mesh regions: "<<num_total_regions_<<std::endl;
			std::cout<<"Number of rows in nodes_ structure: "<<nodes_.size()<<std::endl;
			for( int i = 0; i < nodes_.size(); ++i )
			{
				std::cout<< std::get<0>(nodes_[i]) <<"\t"<< std::get<1>(nodes_[i]) <<std::endl;
			}
		}
	}


	template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	void SplittingDriver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::printNodesToRegion()
	{
		if( 0==comm_->getRank() )
		{
			std::cout<<"Total number of mesh nodes: "<<num_total_nodes_<<std::endl;
			std::cout<<"Total number of mesh regions: "<<num_total_regions_<<std::endl;
			std::cout<<"Number of rows in nodes_ structure: "<<nodes_.size()<<std::endl;
			for( int i = 0; i < nodesToRegion_.size(); ++i )
			{
				std::cout<<"Node "<< std::get<0>(nodesToRegion_[i]) <<"\t belongs to regions: "<< std::get<1>(nodesToRegion_[i]) <<std::endl;
			}
		}
	}

	template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	void SplittingDriver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::printInactive()
	{
		if( maps_.global_map_.empty() )
			std::cout<<"INACTIVE PROC ID: "<<comm_->getRank()<<std::endl;
	}

} //namespace Xpetra

#endif
