#ifndef __Panzer_UtilityAlgs_hpp__
#define __Panzer_UtilityAlgs_hpp__

#include <vector>
#include <functional>

namespace panzer{

/**
 * \brief Using a functor, reorder an array using a order vector.
 *
 * \param[in,out] order Vector that describes the desired order. 
 *                      Note on output, this array will be sorted.
 * \param[in] swapper Functor that swaps entries in the array to be
 *                    sorted. Data being swapped must be captured
 *                    by the functor.
 */
void reorder(std::vector<int> & order,std::function<void(int,int)> swapper);

}

#endif
