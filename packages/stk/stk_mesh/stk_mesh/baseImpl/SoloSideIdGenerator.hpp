
#ifndef SOLOSIDEIDGENERATOR_HPP_
#define SOLOSIDEIDGENERATOR_HPP_

#include <stk_util/environment/ReportHandler.hpp>
#include <stk_mesh/base/Types.hpp>

namespace stk {
namespace mesh {
namespace impl {

class SoloSideIdGenerator
{
public:
    SoloSideIdGenerator(int nProcs, int proc, uint64_t maxSideId)
      : numProcs(nProcs), pseudoOrdinal(6)
    {
        ThrowRequire(maxSideId > 10);
        m_maxPseudoElement = (maxSideId - 10) / 10;

        pseudoElement = proc + 1;
    }

    stk::mesh::EntityId get_solo_side_id()
    {
        stk::mesh::EntityId id = get_solo_side_id_using_formula(pseudoElement, pseudoOrdinal);
        ThrowRequireMsg(pseudoElement <= max_pseudo_element(), "Exhausted solo side ids. Please report to sierra-help@sandia.gov");
        incrementPseudoElementAndOrdinal();
        return id;
    }

    inline uint64_t max_pseudo_element() const  { return m_maxPseudoElement; }

protected:
    inline stk::mesh::EntityId get_solo_side_id_using_formula(unsigned elementId, unsigned sideOrdinal)
    {
        ThrowRequireMsg(elementId <= m_maxPseudoElement, "Exhausted solo side ids for this processor. Please report to sierra-help@sandia.gov");

        //this is the side-id formula used by IO. the "+1" is because IO always uses one-based side ordinals
        return 10*elementId + sideOrdinal + 1;
    }

private:
    void incrementPseudoElementAndOrdinal()
    {
        pseudoOrdinal++;
        if(pseudoOrdinal%10==0)
        {
            pseudoOrdinal = 6;
            pseudoElement += numProcs;
        }
    }
    int numProcs;
    uint64_t m_maxPseudoElement;
    unsigned pseudoOrdinal;
    stk::mesh::EntityId pseudoElement;
};

}
}
}

#endif /* SOLOSIDEIDGENERATOR_HPP_ */
