#include <gtest/gtest.h>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/baseImpl/Visitors.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>

#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/WallTime.hpp>

#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>

#include "stk_io/StkMeshIoBroker.hpp"

namespace
{

double get_cpu_or_wall_time()
{
#if defined(_OPENMP)
  return stk::wall_time();
#else
  return stk::cpu_time();
#endif
}

class FieldVertexWeightSettingsWithSearchForParticles : public stk::balance::GraphCreationSettings
{
public:
  FieldVertexWeightSettingsWithSearchForParticles(stk::mesh::BulkData &stkMeshBulkData,
                                                  const stk::balance::DoubleFieldType &weightField,
                                                  const double defaultWeight,
                                                  bool incrementalRebalance)
    : m_stkMeshBulkData(stkMeshBulkData),
      m_incrementalRebalance(incrementalRebalance)
  {
    m_method = "parmetis";
    setVertexWeightMethod(stk::balance::VertexWeightMethod::FIELD);
    setVertexWeightFieldName(weightField.name());
    setDefaultFieldWeight(defaultWeight);
  }
  virtual ~FieldVertexWeightSettingsWithSearchForParticles() = default;

  using stk::balance::GraphCreationSettings::getToleranceForFaceSearch;

  virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const { return 1.0; }
  virtual bool includeSearchResultsInGraph() const { return true; }
  virtual bool getEdgesForParticlesUsingSearch() const { return true; }
  virtual bool setVertexWeightsBasedOnNumberAdjacencies() const { return false; }
  virtual double getToleranceForParticleSearch() const { return 1.5; }
  virtual double getToleranceForFaceSearch() const { return 0.005; }
  virtual int getGraphVertexWeight(stk::topology type) const { return 1; }
  virtual double getImbalanceTolerance() const { return 1.05; }
  virtual void setDecompMethod(const std::string& input_method) { m_method = input_method;}
  virtual std::string getDecompMethod() const { return m_method; }
  virtual bool incrementalRebalance() const { return m_incrementalRebalance; }

protected:
  FieldVertexWeightSettingsWithSearchForParticles() = delete;
  FieldVertexWeightSettingsWithSearchForParticles(const FieldVertexWeightSettingsWithSearchForParticles&) = delete;
  FieldVertexWeightSettingsWithSearchForParticles& operator=(const FieldVertexWeightSettingsWithSearchForParticles&) = delete;

  const stk::mesh::BulkData & m_stkMeshBulkData;
  bool m_incrementalRebalance;
};

class IncrementalRebalance : public stk::unit_test_util::MeshFixture
{
protected:
  void check_migration()
  {
    size_t num_elements_migrated = calculate_migrated_elements();

    std::ostringstream os;
    os << "P[" << get_bulk().parallel_rank() << "] elements migrated: " << num_elements_migrated << std::endl;
    std::cerr << os.str();
  }

  size_t calculate_migrated_elements()
  {
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, get_meta().locally_owned_part(), elements);

    size_t num_elements_migrated = 0;
    for(const stk::mesh::Entity& element : elements )
    {
      double* data = stk::mesh::field_data(*procOwner, element);
      if ( *data != get_bulk().parallel_rank())
      {
        num_elements_migrated++;
      }
    }
    return num_elements_migrated;
  }

  void set_proc_owner_on_field()
  {
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, get_meta().locally_owned_part(), elements);

    for(const stk::mesh::Entity &element : elements )
    {
      double* data = stk::mesh::field_data(*procOwner, element);
      *data = get_bulk().parallel_rank();
    }
  }

  void set_weights_on_elements(float weight = 9.0)
  {
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, get_meta().locally_owned_part(), elements);

    for(const stk::mesh::Entity &element : elements )
    {
      double* data = stk::mesh::field_data(*weight_field, element);
      *data = weight;
    }
  }

  void set_weights_on_elements(const std::vector<stk::mesh::EntityId>& ids, float newWeight)
  {
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, get_meta().locally_owned_part(), elements);

    for(stk::mesh::EntityId id : ids )
    {
      stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, id);
      if(get_bulk().is_valid(element) && get_bulk().bucket(element).owned())
      {
        double* data = stk::mesh::field_data(*weight_field, element);
        *data = newWeight;
      }
    }
  }

  void change_weights_on_elements()
  {
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, get_meta().locally_owned_part(), elements);

    double* data = stk::mesh::field_data(*weight_field, elements[0]);
    *data = 1;
  }

  void setup_fields()
  {
    double init_value = 1.0;
    weight_field = & get_meta().declare_field<double>(stk::topology::ELEM_RANK, "Weights", 1);
    stk::mesh::put_field_on_mesh(*weight_field, get_meta().universal_part(), &init_value);
    double init_proc = 0.0;
    procOwner = & get_meta().declare_field<double>(stk::topology::ELEM_RANK, "ProcOwner", 1);
    stk::mesh::put_field_on_mesh(*procOwner, get_meta().universal_part(), &init_proc);
  }

  size_t get_global_element_count()
  {
    if(numGlobalElements == 0)
    {
      size_t num_local_elements = stk::mesh::count_entities(get_bulk(), stk::topology::ELEM_RANK, get_meta().locally_owned_part());
      stk::all_reduce_sum(get_comm(), &num_local_elements, &numGlobalElements, 1);
    }
    return numGlobalElements;
  }

  stk::mesh::Field<double>* get_weight_field() { return weight_field; }

  void move_coordinates(const std::string& coordinate_field_name)
  {
    const stk::mesh::FieldBase * coord = get_meta().get_field(stk::topology::NODE_RANK, coordinate_field_name);

    stk::mesh::EntityVector nodes;
    stk::mesh::get_entities(get_bulk(), stk::topology::NODE_RANK, get_meta().locally_owned_part() | get_meta().globally_shared_part(), nodes);

    double rotation = 15 * M_PI / 180.0;
    double dx = 3;
    double dy = 3;
    double dz = 3;

    for(stk::mesh::Entity node : nodes)
    {
      double *xyz = static_cast<double *>(stk::mesh::field_data(*coord, node));
      double newx = xyz[0]*cos(rotation) + xyz[1]*sin(rotation) + dx;
      double newy = -xyz[0]*sin(rotation) + xyz[1]*cos(rotation) + dy;
      double newz = xyz[2] + dz;
      xyz[0] = newx;
      xyz[1] = newy;
      xyz[2] = newz;
    }
  }

  void print_locally_owned_elements(std::ostream& out)
  {
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, get_meta().locally_owned_part(), elements);

    out << "P" << get_bulk().parallel_rank() << " locally owned elements = ";
    for(const stk::mesh::Entity &element : elements )
    {
      out << get_bulk().identifier(element) << ", ";
    }
    out << std::endl;
  }

  void run_test(const std::string& method)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    setup_fields();
    stk::io::fill_mesh("generated:1x1x20", get_bulk());

    std::ofstream out("elements." + std::to_string(get_bulk().parallel_rank()), std::ios_base::app);
    out << "after mesh read:" << std::endl;
    print_locally_owned_elements(out);

    stk::mesh::Selector selector = get_meta().universal_part();

    set_proc_owner_on_field();
    set_weights_on_elements();

    const double defaultVertexWeight = 0.0;
    stk::balance::FieldVertexWeightSettings graphSettings(get_bulk(), *get_weight_field(), defaultVertexWeight);
    graphSettings.setDecompMethod(method);
    stk::balance::balanceStkMesh(graphSettings, get_bulk(), {selector});

    out << "after first balance with " << method << std::endl;
    print_locally_owned_elements(out);

    check_migration();
    set_proc_owner_on_field();
    change_weights_on_elements();

    move_coordinates(graphSettings.getCoordinateFieldName());

    stk::balance::balanceStkMesh(graphSettings, get_bulk(), {selector});

    out << "after second balance with " << method << std::endl;
    print_locally_owned_elements(out);
    out.close();
    check_migration();
  }

  void decompose1x1x10beamThenRebalanceWithLast2ElementChanges(const std::string& method)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    setup_fields();
    stk::io::fill_mesh("generated:1x1x10", get_bulk());

    std::ofstream out("elements." + std::to_string(get_bulk().parallel_rank()), std::ios_base::app);
    out << "after mesh read:" << std::endl;
    print_locally_owned_elements(out);

    set_weights_on_elements(1);
    const double defaultVertexWeight = 0.0;
    stk::balance::FieldVertexWeightSettings graphSettings(get_bulk(), *get_weight_field(), defaultVertexWeight);
    graphSettings.setDecompMethod(method);
    stk::mesh::Selector selector = get_meta().universal_part();
    stk::balance::balanceStkMesh(graphSettings, get_bulk(), {selector});

    out << "after first balance:" << std::endl;
    print_locally_owned_elements(out);

    set_weights_on_elements({9, 10}, 2.0);
    set_proc_owner_on_field();
    stk::balance::balanceStkMesh(graphSettings, get_bulk(), {selector});

    out << "after second balance:" << std::endl;
    print_locally_owned_elements(out);
    out.close();

    check_migration();
    size_t num_elements_migrated_to_me = calculate_migrated_elements();
    if(get_bulk().parallel_rank() == 0)
      EXPECT_EQ(1u, num_elements_migrated_to_me);
    else
      EXPECT_EQ(0u, num_elements_migrated_to_me);
  }

  void decompose1x9x9beamThenRebalanceWithLast2ElementChanges(const std::string& method1, const std::string& method2)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    setup_fields();
    stk::io::fill_mesh("generated:1x9x9", get_bulk());

    std::ofstream out("elements." + std::to_string(get_bulk().parallel_rank()));
    out << "after mesh read:" << std::endl;
    print_locally_owned_elements(out);

    set_weights_on_elements(1);
    set_proc_owner_on_field();
    const double defaultVertexWeight = 0.0;
    stk::balance::FieldVertexWeightSettings graphSettings(get_bulk(), *get_weight_field(), defaultVertexWeight);
    graphSettings.setDecompMethod(method1);
    if (get_bulk().parallel_rank()==0) std::cerr << "Decomposition method = " << method1 << std::endl;
    stk::mesh::Selector selector = get_meta().universal_part();
    stk::balance::balanceStkMesh(graphSettings, get_bulk(), {selector});

    out << "after first balance:" << std::endl;
    print_locally_owned_elements(out);
    check_migration();

    set_weights_on_elements({9, 10}, 2.0);
    set_proc_owner_on_field();
    graphSettings.setDecompMethod(method2);

    if (get_bulk().parallel_rank()==0) std::cerr << "Decomposition method = " << method2 << std::endl;
    stk::balance::balanceStkMesh(graphSettings, get_bulk(), {selector});

    out << "after second balance:" << std::endl;
    print_locally_owned_elements(out);
    check_migration();

    out.close();

    size_t num_elements_migrated_to_me = calculate_migrated_elements();
    if (get_bulk().parallel_rank() == 0) {
      EXPECT_EQ(24u, num_elements_migrated_to_me);
    }
    else {
      EXPECT_EQ(12u, num_elements_migrated_to_me);
    }
  }

  void decomposeWithRcbThenParmetisAndCheckMigration()
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    setup_fields();
    stk::io::fill_mesh("generated:10x10x10", get_bulk());

    std::ofstream out("elements." + std::to_string(get_bulk().parallel_rank()), std::ios_base::app);
    out << "after mesh read:" << std::endl;
    print_locally_owned_elements(out);

    set_weights_on_elements(1);
    const double defaultVertexWeight = 0.0;
    stk::balance::FieldVertexWeightSettings graphSettings(get_bulk(), *get_weight_field(), defaultVertexWeight);
    graphSettings.setDecompMethod("rcb");
    stk::mesh::Selector selector = get_meta().universal_part();
    stk::balance::balanceStkMesh(graphSettings, get_bulk(), {selector});

    out << "after first rebalance:" << std::endl;
    print_locally_owned_elements(out);

    set_proc_owner_on_field();
    graphSettings.setDecompMethod("parmetis");
    stk::balance::balanceStkMesh(graphSettings, get_bulk(), {selector});
    out << "after second rebalance:" << std::endl;
    print_locally_owned_elements(out);
    out.close();

    check_migration();

    size_t num_elements_migrated_to_me = calculate_migrated_elements();
    EXPECT_TRUE(0 < num_elements_migrated_to_me);
  }

  size_t count_global_non_particle_elements()
  {
    size_t numGlobalElements;
    size_t num_local_elements = stk::mesh::count_entities(get_bulk(), stk::topology::ELEM_RANK, get_meta().locally_owned_part() & !get_meta().get_topology_root_part(stk::topology::PARTICLE));
    stk::all_reduce_sum(get_comm(), &num_local_elements, &numGlobalElements, 1);
    return numGlobalElements;
  }

  double count_average_global_moved_elements(size_t localMovedElements)
  {
    size_t globalMovedElements;
    stk::all_reduce_sum(get_comm(), &localMovedElements, &globalMovedElements, 1);
    double averageMovedElements = globalMovedElements/get_bulk().parallel_size();
    return averageMovedElements;
  }

  size_t count_max_global_moved_elements(size_t localMovedElements)
  {
    size_t globalMaxMovedElements;
    stk::all_reduce_max(get_comm(), &localMovedElements, &globalMaxMovedElements, 1);
    return globalMaxMovedElements;
  }

  size_t count_local_elements()
  {
    return stk::mesh::count_entities(get_bulk(), stk::topology::ELEM_RANK, get_meta().locally_owned_part());
  }

  stk::mesh::EntityVector get_non_particle_elements_inside_sphere(double* center, double radius)
  {
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, get_meta().locally_owned_part() & !get_meta().get_topology_root_part(stk::topology::PARTICLE), elements);
    stk::mesh::EntityVector selectedElements;
    for(size_t i=0 ; i<elements.size() ; ++i)
    {
      double centroid[3] = {0,0,0};
      compute_centroid(elements[i],centroid);
      if (distance(center,centroid) <= radius)
      {
        selectedElements.push_back(elements[i]);
      }
    }
    return selectedElements;
  }

  double distance(double* A, double* B)
  {
    return sqrt( (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) + (A[2]-B[2])*(A[2]-B[2]) );
  }

  void compute_centroid(stk::mesh::Entity element, double *centroid)
  {
    const stk::mesh::Entity* nodes = get_bulk().begin_nodes(element);
    const unsigned numNodesThisEntity = get_bulk().num_nodes(element);
    if(get_bulk().is_valid(element))
    {
      for(unsigned j=0; j<numNodesThisEntity; j++)
      {
        if (get_bulk().is_valid(nodes[j]))
        {
          double *coordDataForNode = static_cast<double*>(stk::mesh::field_data(*get_meta().coordinate_field(),nodes[j]));
          centroid[0] += coordDataForNode[0];
          centroid[1] += coordDataForNode[1];
          centroid[2] += coordDataForNode[2];
        }
      }
      centroid[0] /= numNodesThisEntity;
      centroid[1] /= numNodesThisEntity;
      centroid[2] /= numNodesThisEntity;
    }
  }

  void set_weight(stk::mesh::Entity element, double weight)
  {
    double* data = stk::mesh::field_data(*weight_field, element);
    *data = weight;
  }

  void set_coords(stk::mesh::Entity node, double* coords)
  {
    double* data = static_cast<double*>(stk::mesh::field_data(*get_meta().coordinate_field(), node));
    data[0] = coords[0];
    data[1] = coords[1];
    data[2] = coords[2];
  }

  void balance_mesh_with_search_for_particles(std::string decompMethod, bool incrementalRebalance)
  {
    const double defaultVertexWeight = 0.0;
    FieldVertexWeightSettingsWithSearchForParticles graphSettings(get_bulk(), *get_weight_field(), defaultVertexWeight, incrementalRebalance);
    graphSettings.setDecompMethod(decompMethod);
    stk::mesh::Selector selector = get_meta().universal_part();
    stk::balance::balanceStkMesh(graphSettings, get_bulk(), {selector});
  }

  void destroy_element_and_lower(stk::mesh::Entity element)
  {
    stk::mesh::EntityVector entities;
    stk::mesh::impl::StoreInVector<stk::mesh::EntityVector> siv(entities);
    stk::mesh::impl::VisitClosure(get_bulk(),element,siv);
    for (stk::mesh::Entity entity : entities)
    {
      get_bulk().destroy_entity(entity);
    }
  }

  void convert_elements_to_particles(const stk::mesh::EntityVector& elements, int numParticlesPerElement, double particleWeight, size_t& entityId)
  {
    stk::mesh::Part& particleTopologyPart = get_meta().get_topology_root_part(stk::topology::PARTICLE);
    stk::mesh::Part& nodeTopologyPart = get_meta().get_topology_root_part(stk::topology::NODE);
    stk::mesh::Part* block2Part = get_meta().get_part("block_2");
    get_bulk().modification_begin();
    for (size_t i=0 ; i<elements.size() ; ++i)
    {
      double centroid[3] = {0,0,0};
      compute_centroid(elements[i], centroid);
      destroy_element_and_lower(elements[i]);
      for (int particleIndex=0 ; particleIndex < numParticlesPerElement ; ++particleIndex)
      {
        stk::mesh::Entity node = get_bulk().declare_node(++entityId, stk::mesh::ConstPartVector{&nodeTopologyPart});
        stk::mesh::Entity particle = get_bulk().declare_element(entityId, stk::mesh::ConstPartVector{&particleTopologyPart,block2Part});
        get_bulk().declare_relation(particle, node, 0);

        set_coords(node, centroid);
        set_weight(particle, particleWeight);
      }

    }
    get_bulk().modification_end();
  }

  void setup_io_parts()
  {
    stk::mesh::Part& block2 = get_meta().declare_part_with_topology("block_2", stk::topology::PARTICLE);
    stk::io::put_io_part_attribute(block2);
  }

  unsigned get_number_of_digits (unsigned i)
  {
    return i > 0 ? (int) log10 ((double) i) + 1 : 1;
  }


  void write_mesh_to_exodus(std::string fileString, int iteration, double simulationTime)
  {
    stk::io::StkMeshIoBroker writer(get_bulk().parallel());
    writer.set_bulk_data(get_bulk());
    std::string filename = fileString + ".e";
    if (iteration > 0)
    {
      filename += "-s";
      unsigned numZeros = 4-get_number_of_digits(iteration);
      for (unsigned i=0 ; i<numZeros ; ++i) { filename += "0"; }
      filename += std::to_string(iteration++);
    }
    //        if (get_bulk().parallel_rank() == 0) { std::cerr << "filename = " << filename << std::endl; }
    size_t outputHandle = writer.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    writer.add_field(outputHandle, *procOwner);
    writer.add_field(outputHandle, *weight_field);
    writer.begin_output_step(outputHandle, simulationTime);
    writer.write_defined_output_fields(outputHandle);
    writer.end_output_step(outputHandle);
  }

  void incrementally_convert_elements_to_particles(std::string meshFile, std::string decompMethod, double hexWeight, double center[3], double deltaR, int numParticlesPerElement, double particleWeight, int maxIters)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    setup_io_parts();
    setup_fields();
    stk::io::fill_mesh(meshFile, get_bulk());

    std::ofstream out("elements." + std::to_string(get_bulk().parallel_size()) + "." + std::to_string(get_bulk().parallel_rank()));
    out << "numProcs=" << get_bulk().parallel_size()
        << " meshFile=" << meshFile
        << " decompMethod=" << decompMethod
        << " center=(" << center[0] << "," << center[1] << "," << center[2]
        << ") deltaR=" << deltaR
        << " numParticlesPerElement=" << numParticlesPerElement
        << " particleWeight=" << particleWeight
        << " maxIters=" << maxIters
        << std::endl;

    set_weights_on_elements(1);
    bool incrementalRebalance = false;
    balance_mesh_with_search_for_particles(decompMethod, incrementalRebalance);
    incrementalRebalance = true;
    set_proc_owner_on_field();
    double simulationTime = 0.0;
    int iteration = 0;
    write_mesh_to_exodus(decompMethod, iteration, simulationTime);

    size_t numTotalElements = get_global_element_count();
    size_t entityId = (get_bulk().parallel_rank()+1)*numTotalElements*8*numParticlesPerElement;

    double radius = 0.0;
    double totalRebalanceTime = 0.0;
    while (count_global_non_particle_elements() > 0)
    {
      radius += deltaR;
      stk::mesh::EntityVector nonParticleElements = get_non_particle_elements_inside_sphere(center, radius);
      out << "iteration=" << ++iteration
          << " radius=" << radius
          << " elementsConverted=" << nonParticleElements.size();
      convert_elements_to_particles(nonParticleElements, numParticlesPerElement, particleWeight, entityId);
      size_t oldNumElements = count_local_elements();
      {
        const double timeStart = get_cpu_or_wall_time();
        // time this operation and report in output
        balance_mesh_with_search_for_particles(decompMethod, incrementalRebalance);
        const double balanceTime = get_cpu_or_wall_time()-timeStart;
        out << " balanceTime=" << balanceTime;
        totalRebalanceTime += balanceTime;
      }
      const size_t newNumElements = count_local_elements();
      const size_t numElementsGained = calculate_migrated_elements();
      const size_t numElementsLost = oldNumElements - (newNumElements - numElementsGained);
      out << " old=" << oldNumElements
          << " new=" << newNumElements
          << " gained=" << numElementsGained
          << " lost=" << numElementsLost;
      const double numAverageGlobalElementsMoved = count_average_global_moved_elements(numElementsGained+numElementsLost);
      const size_t numMaxGlobalElementsMoved = count_max_global_moved_elements(numElementsGained+numElementsLost);
      out << " maxGlobalMoved=" << numMaxGlobalElementsMoved;
      out << " averageGlobalMoved=" << numAverageGlobalElementsMoved;
      out << " max/avg=" << static_cast<double>(numMaxGlobalElementsMoved == 0 ? 0.0 : numMaxGlobalElementsMoved/numAverageGlobalElementsMoved);
      out << std::endl;
      set_proc_owner_on_field();

      simulationTime += 1.0;
      write_mesh_to_exodus(decompMethod, iteration, simulationTime);

      if (iteration >= maxIters) break;
    }
    out << "totalRebalanceTime=" << totalRebalanceTime << std::endl;
  }

private:
  stk::mesh::Field<double>* weight_field = nullptr;
  stk::mesh::Field<double>* procOwner = nullptr;
  size_t numGlobalElements = 0;
};

TEST_F(IncrementalRebalance, parmetis)
{
  if(stk::parallel_machine_size(get_comm())<=4)
    run_test("parmetis");
}

TEST_F(IncrementalRebalance, rcb)
{
  if(stk::parallel_machine_size(get_comm())<=4)
    run_test("rcb");
}

TEST_F(IncrementalRebalance, rib)
{
  if(stk::parallel_machine_size(get_comm())<=4)
    run_test("rib");
}

TEST_F(IncrementalRebalance, multijagged)
{
  if(stk::parallel_machine_size(get_comm())<=4)
    run_test("multijagged");
}

TEST_F(IncrementalRebalance, rcb_case1)
{
  if(stk::parallel_machine_size(get_comm())==2)
    decompose1x1x10beamThenRebalanceWithLast2ElementChanges("rcb");
}

TEST_F(IncrementalRebalance, parmetis_case1)
{
  if(stk::parallel_machine_size(get_comm())==2)
    decompose1x1x10beamThenRebalanceWithLast2ElementChanges("parmetis");
}

#if !defined(__APPLE__)
TEST_F(IncrementalRebalance, rib_then_parmetis_case2)
{
  if(stk::parallel_machine_size(get_comm())==3)
    decompose1x9x9beamThenRebalanceWithLast2ElementChanges("rib", "parmetis");
}
#endif

//TEST_F(IncrementalRebalance, hypergraph_case1)
//{
//    if(stk::parallel_machine_size(get_comm())==2)
//        decompose1x1x10beamThenRebalanceWithLast2ElementChanges("zoltan"); // hypergraph
//}

TEST_F(IncrementalRebalance, rcb_to_parmetis_check_migration)
{
  if(stk::parallel_machine_size(get_comm())==5)
    decomposeWithRcbThenParmetisAndCheckMigration();
}

TEST_F(IncrementalRebalance, element_to_particle_conversion_rcb)
{
  double hexWeight = 1.0;
  double center[3] = {0,0,0};
  double deltaR = 100.0; // change to 1.0 for a nice movie
  int numParticlesPerElement = 1;
  double particleWeight = 10.0;
  int maxIters = 1000;
  incrementally_convert_elements_to_particles("generated:1x100x100", "rcb", hexWeight, center, deltaR, numParticlesPerElement, particleWeight, maxIters);
}

TEST_F(IncrementalRebalance, element_to_particle_conversion_multijagged)
{
  double hexWeight = 1.0;
  double center[3] = {0,0,0};
  double deltaR = 100.0; // change to 1.0 for a nice movie
  int numParticlesPerElement = 1;
  double particleWeight = 10.0;
  int maxIters = 1000;
  incrementally_convert_elements_to_particles("generated:1x100x100", "multijagged", hexWeight, center, deltaR, numParticlesPerElement, particleWeight, maxIters);
}

TEST_F(IncrementalRebalance, element_to_particle_conversion_parmetis)
{
  double hexWeight = 1.0;
  double center[3] = {0,0,0};
  double deltaR = 100.0; // change to 1.0 for a nice movie
  int numParticlesPerElement = 1;
  double particleWeight = 10.0;
  int maxIters = 1000;
  incrementally_convert_elements_to_particles("generated:1x100x100", "parmetis", hexWeight, center, deltaR, numParticlesPerElement, particleWeight, maxIters);
}

}
