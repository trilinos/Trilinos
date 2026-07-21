#ifndef stk_mesh_GhostingCache_hpp
#define stk_mesh_GhostingCache_hpp

#include <vector>
#include <string>
namespace stk::mesh {
class Ghosting;
class Part;
class MetaData;
class BulkData;
}
namespace stk::mesh::impl {
class GhostingCache
{
  public:
    explicit GhostingCache(BulkData& bulk_data);

    ~GhostingCache();

    GhostingCache(const GhostingCache&) = delete;

    GhostingCache& operator=(const GhostingCache&) = delete;

    Ghosting* create_ghosting(const std::string& name);

    void destroy_ghosting(Ghosting* ghosting);

    const std::vector<Ghosting*>& get_ghostings() const;

    Part* get_ghosting_part(const Ghosting* ghosting) const;

    bool is_ghost_part(Part* part) const;

  private:
    Ghosting* get_ghosting_from_cache(const std::string& name);

    Ghosting* create_new_ghosting(const std::string& name);

    void create_new_part(const std::string& ghosting_name);

    void check_name_debug(const std::string& name);


    MetaData& m_meta_data;
    BulkData& m_bulk_data;

    std::vector<Ghosting*> m_ghosting;
    std::vector<Ghosting*> m_deleted_ghostings;
    std::vector<Part*> m_ghost_parts;
};

}  // namespace

#endif