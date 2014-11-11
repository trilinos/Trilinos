public: //methods
  PairIterEntityComm entity_comm_map(const EntityKey & key) const { return m_entity_comm_map.comm(key); }
  PairIterEntityComm entity_comm_map_shared(const EntityKey & key) const { return m_entity_comm_map.shared_comm_info(key); }
  PairIterEntityComm entity_comm_map(const EntityKey & key, const Ghosting & sub ) const { return m_entity_comm_map.comm(key,sub); }

  int entity_comm_map_owner(const EntityKey & key) const;

  // Comm-related convenience methods

  bool in_shared(EntityKey key) const { return !entity_comm_map_shared(key).empty(); }
  bool in_shared(EntityKey key, int proc) const;
  bool in_receive_ghost( EntityKey key ) const;
  bool in_receive_ghost( const Ghosting & ghost , EntityKey entity ) const;
  bool in_send_ghost( EntityKey key) const;
  bool in_send_ghost( EntityKey key , int proc ) const;
  bool is_ghosted_onto_another_proc( EntityKey key ) const;
  bool is_ghosted_onto_proc( EntityKey key, int otherProc ) const;
  bool in_ghost( const Ghosting & ghost , EntityKey key , int proc ) const;
  void comm_procs( EntityKey key, std::vector<int> & procs ) const; //shared and ghosted entities
  void comm_shared_procs( EntityKey key, std::vector<int> & procs ) const; // shared entities
  void shared_procs_intersection( std::vector<EntityKey> & keys, std::vector<int> & procs ) const;
  void comm_procs( const Ghosting & ghost ,
                   EntityKey key, std::vector<int> & procs ) const;

  void internal_change_owner_in_comm_data(const EntityKey& key, int new_owner);
  void internal_sync_comm_list_owners();

protected: //methods
  void update_comm_list(const std::vector<stk::mesh::Entity>& shared_modified);
  bool entity_comm_map_insert(Entity entity, const EntityCommInfo & val) { return m_entity_comm_map.insert(entity_key(entity), val, parallel_owner_rank(entity)); }
  bool entity_comm_map_erase(  const EntityKey & key, const EntityCommInfo & val) { return m_entity_comm_map.erase(key,val); }
  bool entity_comm_map_erase(  const EntityKey & key, const Ghosting & ghost) { return m_entity_comm_map.erase(key,ghost); }
  void entity_comm_map_clear_ghosting(const EntityKey & key ) { m_entity_comm_map.comm_clear_ghosting(key); }
  void entity_comm_map_clear(const EntityKey & key) { m_entity_comm_map.comm_clear(key); }

  void internal_update_fast_comm_maps();

  void internal_print_comm_map( std::string title);

private: //methods
  //no private methods?

public: //data
  //no public data

protected: //data
  EntityCommDatabase m_entity_comm_map;

private: //data
  EntityCommListInfoVector m_entity_comm_list;
  VolatileFastSharedCommMap m_volatile_fast_shared_comm_map;
