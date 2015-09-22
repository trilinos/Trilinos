#ifndef STK_COMMLISTUPDATER_HPP
#define STK_COMMLISTUPDATER_HPP

#include <stk_mesh/base/EntityCommListInfo.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>

namespace stk {
namespace mesh {

class CommListUpdater  : public CommMapChangeListener {
public:
    CommListUpdater(EntityCommListInfoVector& comm_list) : m_comm_list(comm_list) {}
    virtual ~CommListUpdater(){}

    virtual void removedKey(const EntityKey& key) {
        EntityCommListInfoVector::iterator iter =
                std::lower_bound(m_comm_list.begin(), m_comm_list.end(), key);
        if (iter != m_comm_list.end() && iter->key == key) {
            iter->entity_comm = nullptr;
        }
    }

private:
  EntityCommListInfoVector& m_comm_list;
};

}
}

#endif
