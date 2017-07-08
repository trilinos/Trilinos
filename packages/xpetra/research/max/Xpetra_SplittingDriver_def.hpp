#ifndef XPETRA_SPLITTINGDRIVER_DEF_HPP
#define XPETRA_SPLITTINGDRIVER_DEF_HPP

#include "Xpetra_SplittingDriver_decl.hpp"

namespace Xpetra{

	template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	SplittingDriver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SplittingDriver(const std::string file_name, Teuchos::RCP< const Teuchos::Comm<int> > comm)
	{
		std::ifstream input_file_(file_name, std::ifstream::in);
		std::string  line;
		TEUCHOS_TEST_FOR_EXCEPTION( !input_file_.good(), Exceptions::RuntimeError, "Can not read \"" << file_name << "\"");

		GlobalOrdinal line_index = 0;
	
		if(comm->getRank()==0)
		{
		while ( std::getline (input_file_,line) )
		{
			if( 0!=line_index )
			{
				std::istringstream is( line );
				GlobalOrdinal number;
 				std::vector<GlobalOrdinal> node;
				std::tuple<GlobalOrdinal, GlobalOrdinal> node_region;
				while(is>>number)
				{
					node.emplace_back(number);
				}
				TEUCHOS_TEST_FOR_EXCEPTION( node.size()!=2, Exceptions::RuntimeError, "The node information must be a couple of integers: Node index + Region idnex \n");
				node_region = std::make_tuple(node[0], node[1]);
				node_.emplace_back(node_region);
				node.empty();
			}
			line_index++;
		}
		}
		input_file_.close();

	}

	template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
	void SplittingDriver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::printView()
	{
		for( int i = 0; i < node_.size(); ++i )
		{
			std::cout<< std::get<0>(node_[i]) <<"\t"<< std::get<1>(node_[i]) <<std::endl;
		}
	}

} //namespace Xpetra

#endif
