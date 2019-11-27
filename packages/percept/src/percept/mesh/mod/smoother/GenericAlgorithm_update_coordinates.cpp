#include <percept/mesh/mod/smoother/GenericAlgorithm_update_coordinates.hpp>

#include <percept/PerceptUtils.hpp>
#include <percept/PerceptMesh.hpp>

#include <percept/structured/BlockStructuredGrid.hpp>

#include <percept/mesh/mod/smoother/ReferenceMeshSmootherConjugateGradient.hpp>

namespace percept {

template<>
GenericAlgorithm_update_coordinates<STKMesh>::
GenericAlgorithm_update_coordinates(RefMeshSmoother *rms, PerceptMesh *eMesh, Double alpha_in)
  : m_rms(rms), m_eMesh(eMesh), alpha(alpha_in)
{
	 stk::diag::Timer     cumulative_timer(eMesh->getProperty("in_filename"), rootTimerStructured());
	 stk::diag::Timer constrGATM("GATM <STKMesh> constructor",cumulative_timer);
	 stk::diag::TimeBlock my_time_block(constrGATM);

  dmax=m_rms->m_dmax;
  dmax_relative=m_rms->m_dmax_relative;
  m_calc_deltas=false;

  coord_field = m_eMesh->get_coordinates_field();
  coord_field_current   = coord_field;
  coord_field_lagged  = m_eMesh->get_field(stk::topology::NODE_RANK, "coordinates_lagged");

  cg_s_field = m_eMesh->get_field(stk::topology::NODE_RANK, "cg_s");
  cg_edge_length_field = eMesh->get_field(stk::topology::NODE_RANK, "cg_edge_length");

  on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
  on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );
  spatialDim = m_eMesh->get_spatial_dim();

  {
    std::vector< const stk::mesh::FieldBase *> fields;
    fields.push_back(cg_s_field);
    fields.push_back(coord_field);
    // stk::mesh::copy_owned_to_shared(*m_eMesh->get_bulk_data(), fields);
    // stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->aura_ghosting(), fields);
    MTcommFields<STKMesh>(fields, m_eMesh);
  }

  // cache coordinates
  m_eMesh->copy_field("coordinates_lagged", "coordinates");

  nodes.resize(0);

  // node loop - update positions
  {
    const stk::mesh::BucketVector & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->node_rank() );
    for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        // update local and globally shared
        if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_nodes_in_bucket = bucket.size();

            for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
              {
                stk::mesh::Entity node = bucket[i_node];
                nodes.push_back(node);
              }
          }
      }
  }

  //filter out fixed nodes
  std::vector<typename STKMesh::MTNode> unFixedNodes;
  unFixedNodes.resize(nodes.size());
  std::pair<bool,int> fixed;
  int64_t numUnFixed = 0;

  for (int64_t iNode = 0; iNode < (int64_t)nodes.size(); ++iNode)
	{
		fixed = m_rms->get_fixed_flag(nodes[iNode]);
		if(!fixed.first)
		{
			unFixedNodes[numUnFixed]=nodes[iNode];
			numUnFixed++;
		}
	}
 nodes.resize(numUnFixed);

 for (int64_t iNode = 0; iNode < numUnFixed; ++iNode)
	{ nodes[iNode]=unFixedNodes[iNode]; } //only operate on unFixedNodes

}//GATM constrcr <STKMesh>


template<>
KOKKOS_INLINE_FUNCTION
void GenericAlgorithm_update_coordinates<STKMesh>::
operator()(const int64_t& index) const
{
	  STKMesh::MTNode node = nodes[index];
	  double coord_current[3];
	  double cg_s[3];
	  get_field<STKMesh>(coord_current, spatialDim, m_eMesh, coord_field_current, node);
	  get_field<STKMesh>(cg_s, spatialDim, m_eMesh, cg_s_field, node);
	  double cg_edge_length[1];

	  //double coord_project[3] = {0,0,0};
	  for (int i=0; i < spatialDim; i++)
	    {
	      Double dt = alpha * cg_s[i];
	      coord_current[i] += dt;
	      if(m_calc_deltas)
	      {
	          get_field<STKMesh>(cg_edge_length, 1, m_eMesh, cg_edge_length_field, node);
	          m_rms->m_dmax = std::max(std::fabs(dt), m_rms->m_dmax);
	          m_rms->m_dmax_relative = std::max(std::fabs(dt)/cg_edge_length[0], m_rms->m_dmax_relative);
	      }
	    }
	  set_field<STKMesh>(coord_current, spatialDim, m_eMesh, coord_field_current, node);
}


template<>
void GenericAlgorithm_update_coordinates<STKMesh>::
run(bool calc_deltas)
{
      m_calc_deltas=calc_deltas;
      if(m_calc_deltas)
      {//ensure deltas are updated before calculating them
          m_rms->m_dmax=0;
          m_rms->m_dmax_relative=0;
          dmax=m_rms->m_dmax;
          dmax_relative=m_rms->m_dmax_relative;
      }

	  {
	   	  std::vector< const typename STKMesh::MTField *> fields;
	   	  fields.push_back(cg_s_field);
	   	  fields.push_back(coord_field);
	      MTcommFields<STKMesh>(fields, m_eMesh);
	  }

	  // cache coordinates
	  m_eMesh->copy_field("coordinates_lagged", "coordinates");

	  for (int64_t index = 0; index < (int64_t)nodes.size(); ++index)
	   {
	     (*this)(index);
	   }


}

template<>
GenericAlgorithm_update_coordinates<StructuredGrid>::
GenericAlgorithm_update_coordinates(
		RefMeshSmoother *rms, PerceptMesh *eMesh, Double alpha_in) :
		m_rms(rms), m_eMesh(eMesh), alpha(alpha_in) {
	std::vector<StructuredGrid::MTNode> nodesVec;
	//madbrew: for some reason using the member nodes to populate my A4MD's causes backward scaling in PF
	//using this however does not
	stk::diag::Timer cumulative_timer(eMesh->getProperty("in_filename"),
			rootTimerStructured());
	stk::diag::Timer constrGATM("GATM <StrctGrd> constructor",
			cumulative_timer);
	stk::diag::TimeBlock my_time_block(constrGATM);


	dmax=m_rms->m_dmax;
	dmax_relative=m_rms->m_dmax_relative;
	m_calc_deltas=false;

	std::shared_ptr<BlockStructuredGrid> bsg =
			m_eMesh->get_block_structured_grid();
	coord_field = bsg->m_fields["coordinates"].get();	//original coordinates
	coord_field_current = coord_field; //most up to date coords
	coord_field_lagged = bsg->m_fields["coordinates_lagged"].get(); //prev iteration's coords
	cg_s_field = bsg->m_fields["cg_s"].get();
	cg_edge_length_field= bsg->m_fields["cg_edge_length"].get();

	dmax_candidates.resize(bsg->m_sblocks.size());
	dmax_relative_candidates.resize(bsg->m_sblocks.size());

	//on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
	//on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );
	spatialDim = 3;    //madbrew: should this always be 3?

	{
		std::vector<const typename StructuredGrid::MTField *> fields;
		fields.push_back(cg_s_field);
		fields.push_back(coord_field);
		// stk::mesh::copy_owned_to_shared(*m_eMesh->get_bulk_data(), fields);
		// stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->aura_ghosting(), fields);
		MTcommFields<StructuredGrid>(fields, m_eMesh); //ensure consistency between stk::mesh parts
	}

	// cache coordinates
	m_eMesh->copy_field("coordinates_lagged", "coordinates");

	for (unsigned iter =0;iter<bsg->m_sblocks.size();iter++) { //for each unique blockID
			A4DMD<StructuredGrid> blockMD;
			Kokkos::View<unsigned**, DataLayout,MemSpace>::HostMirror interimNodeIDs;
			blockMD.blkIndex = (/**iter*/iter);

			//get unfixed nodes
			nodesVec.resize(0);
			bsg->get_nodes_of_sb(nodesVec,iter);
			std::vector<StructuredGrid::MTNode> unFixedNodes;
			unFixedNodes.resize(nodesVec.size());
			unsigned numUnFixed = 0;
			for(unsigned iNode=0;iNode<nodesVec.size();iNode++)
			{
			    std::pair<bool,int> fixed = get_fixed_flag_sgrid(nodesVec[iNode],m_rms->get_sgrid_select());
		        if (!fixed.first) {
		            unFixedNodes[numUnFixed] = nodesVec[iNode];
		            numUnFixed++;
		        }
			}
			nodesVec.resize(numUnFixed);
			for (int64_t iNode = 0; iNode < numUnFixed; ++iNode)
			        nodesVec[iNode] = unFixedNodes[iNode];
			//get unfixed nodes

			blockMD.numNodes = numUnFixed;
			Kokkos::resize(blockMD.blockNodes,blockMD.numNodes,3/*dimension*/); //resize view to fit number of nodes inside of it
			Kokkos::resize(dmax_candidates[/**iter*/iter],blockMD.numNodes);//resize delta maxes to account for all node movements inside a block
			Kokkos::resize(dmax_relative_candidates[/**iter*/iter],blockMD.numNodes);

			interimNodeIDs = Kokkos::create_mirror_view(blockMD.blockNodes);

			for (int64_t iNode = 0; iNode < numUnFixed; ++iNode) {
					interimNodeIDs(iNode,0) = nodesVec[iNode][0];
					interimNodeIDs(iNode,1) = nodesVec[iNode][1];
					interimNodeIDs(iNode,2) = nodesVec[iNode][2];
			}

			Kokkos::deep_copy( blockMD.blockNodes, interimNodeIDs);

			blockMetaDatas.push_front(blockMD); //finally, add to our collection of block meta datas. This will be used to do a blockwise coordinate update in the run function
	} //for each unique blockID
} //GATM constructor SGrid

template<>
KOKKOS_INLINE_FUNCTION
void GenericAlgorithm_update_coordinates<StructuredGrid>::
operator()(const int64_t& index) const
{
//	printf("calling operator\n");
	int i = nodesThisBlock(index,0);//)[0];
	int j = nodesThisBlock(index,1);//)[1];
	int k = nodesThisBlock(index,2);//)[2];

	Double maxDelta = 0.0;

	for(int iDim=0;iDim<spatialDim;iDim++)
	{
		Double dt = alpha * cgsf(i,j,k,iDim);
		cfc(i,j,k,iDim) += dt;
		//this is predicated on the assumption that cgsf, cfl, and cfc have identical i,j,k,block offsets into their respective fields
		if(m_calc_deltas && maxDelta<device_safe_abs_flt(dt))
		        maxDelta = device_safe_abs_flt(dt);
	}

    if(m_calc_deltas)
    {
        Double cg_edge_length = cgelf(i,j,k,0);
        dmax_cands_this_block(index) = device_safe_max<Double>(maxDelta, dmax); //gather all the max delta candidates from each node movement
        dmax_rel_cands_this_block(index) = device_safe_max<Double>(maxDelta/cg_edge_length, dmax_relative);
    }
} //GATM1 sgrid operator

template<>
void GenericAlgorithm_update_coordinates<StructuredGrid>::
run(bool calc_deltas)
{
	stk::diag::Timer cumulative_timer(m_eMesh->getProperty("in_filename"),
			rootTimerStructured());
	stk::diag::Timer runGATM("GATM1 run()", cumulative_timer);
	stk::diag::TimeBlock my_time_block(runGATM);

	m_calc_deltas=calc_deltas;
	if(m_calc_deltas)
	{//ensure deltas are up to date before calculating them
	    m_rms->m_dmax=0;
	    m_rms->m_dmax_relative=0;
        dmax=m_rms->m_dmax;
        dmax_relative=m_rms->m_dmax_relative;
	}
	// cache coordinates
	{
		stk::diag::Timer cpyfldT("copy field within GATM1 run",cumulative_timer);
		stk::diag::TimeBlock cpyfldTB(cpyfldT);
		m_eMesh->copy_field("coordinates_lagged", "coordinates");
	}

	for (auto iBlock = blockMetaDatas.begin(); iBlock != blockMetaDatas.end();
			iBlock++) {
		cfc = *(coord_field_current->m_block_fields[iBlock->blkIndex].get()); //orient the A4D's to the block you want to operate on
		cfl = *(coord_field_lagged->m_block_fields[iBlock->blkIndex].get());
		cgsf = *(cg_s_field->m_block_fields[iBlock->blkIndex].get());
		cgelf = *(cg_edge_length_field->m_block_fields[iBlock->blkIndex].get());
		nodesThisBlock = iBlock->blockNodes;

		dmax_rel_cands_this_block = dmax_relative_candidates[iBlock->blkIndex];
		dmax_cands_this_block = dmax_candidates[iBlock->blkIndex];

		{
			stk::diag::Timer pf("GATM1 parallel_for",cumulative_timer);
			stk::diag::TimeBlock pftb(pf);
			Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,(int64_t)iBlock->numNodes),*this);
		} //kkpf timer

	    if(m_calc_deltas)
	    {   //find max deltas from candidates
	        max_scanner<Double> mscan(dmax_cands_this_block);
	        dmax = std::max( mscan.find_max(),dmax );
	        mscan.candidates = dmax_rel_cands_this_block;
	        dmax_relative = std::max(mscan.find_max(),dmax_relative);
	    }
	}//loop over blocks

    //propogate deltas back to smoother
    m_rms->m_dmax=dmax;
    m_rms->m_dmax_relative=dmax_relative;
} //run GATM1 Sgrid

simplified_gatm_1_BSG::simplified_gatm_1_BSG(RefMeshSmoother *rms, PerceptMesh *eMesh, Double alpha_in)
	  	: m_rms(rms), m_eMesh(eMesh), alpha(alpha_in)
{
	  std::shared_ptr<BlockStructuredGrid> bsg = m_eMesh->get_block_structured_grid();
	  coord_field                              = bsg->m_fields["coordinates"].get();	//original coordinates
	  coord_field_current                      = coord_field; //most up to date coords
	  coord_field_lagged                       = bsg->m_fields["coordinates_lagged"].get(); //prev iteration's coords
	  cg_s_field                               = bsg->m_fields["cg_s"].get();
	  //somehow adding the above stuff made it run like the current code!

	  std::vector<StructuredGrid::MTNode> nodesVec;
	  bsg->get_nodes(nodesVec);
	  //bsg->get_nodes(nodes);
	  numNodes = nodesVec.size();


	  Kokkos::View< StructuredCellIndex*, DataLayout , MemSpace > nds("nodes",nodesVec.size());
	  numBlocks=0;
	  for(unsigned iNode=0;iNode<numNodes;iNode++){
		  nds(iNode)=nodesVec[iNode];
		  if(nodesVec[iNode][3]>numBlocks)
			  numBlocks=nodesVec[iNode][3];//find max block number, that is the number of blocks-1
		  	  	  	  	  	  	  	  	  //assuming block numbers always start at 0
	  }


	  numBlocks++;
	  nodesV = nds;
	  spatialDim=3;
}

void simplified_gatm_1_BSG::run(){

    stk::diag::Timer     cumulative_timer(m_eMesh->getProperty("in_filename"), rootTimerStructured());
    stk::diag::Timer runGATMs("GATM_simple run()",cumulative_timer);
    stk::diag::TimeBlock my_time_block(runGATMs);
    m_eMesh->copy_field("coordinates_lagged", "coordinates");
    {
    	stk::diag::Timer runGATMs_F("GATM_simple total time to loop over all blocks",cumulative_timer);
    	stk::diag::TimeBlock F_time_block(runGATMs_F);
    	for (unsigned iBlockF = 0; iBlockF<numBlocks;iBlockF++){
    		cfcB = *(coord_field_current->m_block_fields[iBlockF].get());
    		cflB = *(coord_field_lagged->m_block_fields[iBlockF].get());
    		cgsfB = *(cg_s_field->m_block_fields[iBlockF].get());
    		{
    			stk::diag::Timer runGATMs_PF("GATM_simple parallel_for",cumulative_timer);
    			stk::diag::TimeBlock PF_time_block(runGATMs_PF);
//        	    			std::cout<<"Size of view" <<nodes.size() <<std::endl;	//THESE YIELD THE SAME RESULTS
//        	    			std::cout<<"number nodes collected" <<nodes.size() <<std::endl;
    			Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,(int64_t)nodesV.size()),*this);
    		}//parallel_for timer													//how to give extra thread private variables?
    	}
    }//all blocks timer
    // cache coordinates

}//run

template class GenericAlgorithm_update_coordinates<STKMesh>;
template class GenericAlgorithm_update_coordinates<StructuredGrid>;

} // namespace percept

