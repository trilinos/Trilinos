#include "Cthulhu_ConfigDefs.hpp"

#ifdef HAVE_CTHULHU_EPETRA

#include "Cthulhu_EpetraMap.hpp"

namespace Cthulhu {
  namespace useEpetra {
  
    /** \brief Non-member function to create a contiguous Map with user-defined weights and a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates EpetraMap
    */
    Teuchos::RCP< const EpetraMap >
    createWeightedContigMapWithNode(int myWeight, global_size_t numElements, 
                                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Kokkos::DefaultNode::DefaultNodeType > &node) { 
      Teuchos::RCP< EpetraMap > map;
      int sumOfWeights, elemsLeft, localNumElements;
      const int numImages = comm->getSize(), 
        myImageID = comm->getRank();
      Teuchos::reduceAll<int>(*comm,Teuchos::REDUCE_SUM,myWeight,Teuchos::outArg(sumOfWeights));
      const double myShare = ((double)myWeight) / ((double)sumOfWeights);
      localNumElements = (int)std::floor( myShare * ((double)numElements) );
      // std::cout << "numElements: " << numElements << "  myWeight: " << myWeight << "  sumOfWeights: " << sumOfWeights << "  myShare: " << myShare << std::endl;
      Teuchos::reduceAll<int>(*comm,Teuchos::REDUCE_SUM,localNumElements,Teuchos::outArg(elemsLeft));
      elemsLeft = numElements - elemsLeft;
      // std::cout << "(before) localNumElements: " << localNumElements << "  elemsLeft: " << elemsLeft << std::endl;
      // i think this is true. just test it for now.
      TEST_FOR_EXCEPT(elemsLeft < -numImages || numImages < elemsLeft);
      if (elemsLeft < 0) {
        // last elemsLeft nodes lose an element
        if (myImageID >= numImages-elemsLeft) --localNumElements;
      }
      else if (elemsLeft > 0) {
        // first elemsLeft nodes gain an element
        if (myImageID < elemsLeft) ++localNumElements;
      }
      // std::cout << "(after) localNumElements: " << localNumElements << std::endl;
      return createContigMapWithNode(numElements,localNumElements,comm,node);
    }

} // useEpetra namespace

}

#endif
