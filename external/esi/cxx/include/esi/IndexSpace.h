#ifndef __ESI_IndexSpace_h
#define __ESI_IndexSpace_h

namespace esi {

/** ESI IndexSpace base class.

    \verbatim
    esi::IndexSpace interface exposes information about the finite dimensional 
    Hilbert space.
    
    Interface method naming conventions:

       {verb}{Set-Restriction}{Set}{Attribute}

       {Set-Restriction} = [ Global | Local ]
       {Set}             = <set_name>
  
       getLocalBlockSetSize   -- size of local( { Block } )
       getLocalBlockSizes     -- size of each member of local( { Block } )
  
    Given that the ESI IndexSpace is all about describing a basis then
    interface methods that access the basis may ommit 'Basis',
    for example:

       getGlobalSize       <=> getGlobalBasisSize
       getLocalSize        <=> getLocalBasisSize
       getLocalColors      <=> getLocalBasisColors
       getLocalIdentifiers <=> getLocalBasisIdentifiers
  
   Sets and attributes:

      Basis         -- set of basis vectors
  
      Color         -- attribute of a basis vector
      ColorSet      -- set of basis vector colors,
                       required to be [0, 1, 2, ..., #colors - 1]
  
      Identifier    -- attribute of a basis vector
      IdentifierSet -- set of basis vector identifiers
                       the element-type of this set is arbitrary
  
      Partition     -- set of basis vectors mapped to an execution context
      PartitionSet  -- set of partitions
  
      Block         -- set of basis vectors in a logical group
      BlockSet      -- set of blocks

    Change log:

      10/29/2001 RLC Added comments about the behavior of zero-length arrays
                 (pointers) which are returned as output args, per Ben Allan's
                 suggestion. Cleaned up the documentation for consistency.
    \endverbatim
*/
template<class Ordinal>
class IndexSpace : public virtual Object
{
  public:
  
  /** Default destructor. */
  virtual ~IndexSpace( void ) {};

  /** Get the global size.
      Query the number of basis vectors in the partitioned
      finite dimensional Hilbert space.
      \param globalSize The number of basis vectors in the partitioned
      finite dimensional Hilbert space.
  */
  virtual ErrorCode getGlobalSize( Ordinal & globalSize ) = 0;

  /** Get the local size.

      Query the number of basis vectors in the partition that
      is mapped to this instance of an ESI IndexSpace within a
      well-defined set of ESI IndexSpace instances.
      There is a consistent requirement s.t. the sum over all
      processors of (localSize) = globalSize.
  */
  virtual ErrorCode getLocalSize( Ordinal & localSize ) = 0;

  /** Get the global color set size.
      
      Each basis vector is associated with a particular FIELD
      from the physics domain (e.g., pressure, temperature).
      Basis vectors may also be classified according to the
      mapping of the discretization to the basis vectors
      (e.g., red-black nodes in a grid or nodes shared between subdomains).
      Such physics or structural information may be used in preconditioning
      or even in the solver.  For example, take the Schur complement the
      "green<->green" submatrix and then iteratively solve the "red<->red"
      submatrix. The attribute used to denote physics/structural
      classification of a basis vector is termed the COLOR of the
      basis vector.
      The color attribute should not be used to denote fine-grained
      blocking.  This attribute should be expressed through the
      esi::IndexSpaceBlock interface (that doesn't yet exist).

      The colors of the global basis vectors must be the set
      [0, 1, ... , colorSize - 1 ] of <ordinal_type>.
      A "colorless" index space may set colorSize == 0 in which case
      the 'getLocalBasisColors' is a no-op.
  */
  virtual ErrorCode getGlobalColorSetSize( Ordinal & colorSetSize ) = 0;

  /** Get the local basis vector colors. 

      Query the color of each local basis vector.

      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned 
      for the data type.
      
  */
  virtual ErrorCode getLocalColors( Ordinal * localBasisColors ) = 0;

  /** Get the local basis vector identifiers.

      Each basis vector may be given a globally unique identifier.
      Globally unique identifiers are of an "arbitrary" type as
      opposed to the ordinal type.  Globally unique identifiers
      may be used to translate data between two different maps
      that have an identical set of identifiers, e.g., a global
      permutation.
      
      <identifier_type> localIdentifiers[ localSize ]

      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned 
      for the data type.
      
  */
  virtual ErrorCode getLocalIdentifiers( Ordinal * localBasisIdentifiers ) = 0;

  /** Query the number of partitions.

      \verbatim
      This information must be consistent with the information produced by
      the following operations:
      1) query the ESI Object's run-time model,
      2) switch over the ESI Object's run-time model possibilities, and
      3) in a run-time model dependent call, query the number of
         execution contexts in the run-time model.
      This interface provides this information independently of
      the run-time model.
      \endverbatim
      \param numberOfPartitions The ordinal number of partitions.
  */
  virtual ErrorCode getGlobalPartitionSetSize( Ordinal & numberOfPartitions ) =
0;

  /** Query which partition is associated with the local execution context.

      \verbatim
      The set of local partitions of an ESI IndexSpace is defined over the interval
      [ 0 , 1 , ... , globalPartitionSetSize - 1 ].  The partition rank of
      the ESI IndexSpace instance may or may not be the same as the identifier for
      the associated execution context.
      \endverbatim
      \param localRank The local execution context (e.g., MPI rank).
  */
  virtual ErrorCode getLocalPartitionRank( Ordinal & localRank ) = 0;

  /** Query the local partition offset.

      \verbatim
      The ordinal offsets of the basis vectors form the set
      [ 0 , 1 , ... , globalSize - 1 ].  This set is partitioned into spans
      and mapped to execution contexts.  The first offset mapped to the
      local execution context is queried with this method.  The span of
      offsets mapped to the local execution context is
      [ globalOffset .. globalOffset + localSize - 1 ].  This information
      must be consistent with querying each localSize and prefixing these
      values in the execution-context rank order.
      \endverbatim
      \param globalOffset The global offset for the local execution
      context (e.g., MPI rank).
  */
  virtual ErrorCode getLocalPartitionOffset( Ordinal & globalOffset ) = 0;

  /** Get the global partition sizes.

      \verbatim
      Query the number of basis vectors in each of the partitions.  This
      information must be consistent with the information produced by the
      following operations:
      1) query the ESI Object's run-time model,
      2) switch over the ESI Object's run-time model possibilities,
      3) gather the localSize from each execution context /
         ESI IndexSpace instance to each execution context, and
      4) Order the gathered sizes by partition rank.
      \endverbatim

      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned 
      for the data type.
      
      \param partitionSizes The array of length 'number of partitions' with
      values of the partition sizes for each execution context (e.g.,
      MPI rank).
  */
  virtual ErrorCode getGlobalPartitionSizes( Ordinal * partitionSizes ) = 0;

  /** Get the global partition offsets .

      Query the offsets of the first basis vector in each partition.  This
      information must be consistent with applying the prefix operation to
      the partitionSize.

      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned 
      for the data type.
      
      \param partitionOffsets The array of length 'number of partitions' with
      values of the partition offsets for each execution context (e.g., MPI rank).
  */
  virtual ErrorCode getGlobalPartitionOffsets( Ordinal * partitionOffsets ) = 0;

};     // esi::IndexSpace class
};     // esi namespace
#endif // __ESI_IndexSpace_h
