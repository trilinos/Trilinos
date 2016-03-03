
Case 1
MappingProblem(
  InputAdapter
  partitioningSolution
  MachineRepresentation=NULL
)
{
  // Create MachineRepresentation if not provided
  // User would have called partitioning problem and provides a solution
  // Mapping vertices are the parts from the partitioning solution
  // Create MappingSolution that can return getRankForPart(part)
}


Case 2
MappingProblem(
  InputAdapter
  MachineRepresentation=NULL
)
{
  // Create MachineRepresentation if not provided
  // Compute mapping vertices based on InputAdapter's partition
  // Assuming input adapter's partition should be used.
  // KDD would use with Exodus/Nemesis input files or PamGen meshes

}


Case 3
MappingProblem(
  InputAdapter
  MachineRepresentation=NULL
)
{
  // Create MachineRepresentation if not provided
  // Call a partitioning algorithm to get mapping vertices that are the parts
  // Mapping vertices are computed from this internal partitioning solution
  // Maybe relevant models can be shared.
  // Similar to what's in PartitioningProblem now and to what LibTopoMap does

}


Case 4
MappingProblem(
  InputAdapter
  MachineRepresentation=NULL
)
{
  // Create MachineRepresentation if not provided
  // Call mapping with mapping vertices == IDs from the input adapter.
  // Similar in spirit to previous case, but runs more slowly since current
  // task mapping is done in serial.
  // Mehmet has done experiments with Scotch, comparing case 3 with 4.
  // Case 3 is faster; case 4 has higher quality.


}


In general, the input Adapters' applyPartitioningSolution method should take an 
optional MappingSolution.

Should MappingSolution provide a re-numbered communicator reflecting the new mapping?
