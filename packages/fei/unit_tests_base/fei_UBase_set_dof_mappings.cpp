
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <fei_DofMapper.hpp>
#include <fei_set_dof_mappings.hpp>


namespace {

template<class LocalOrdinal,class GlobalOrdinal,class DofOrder>
void fill_dof_mapper_1(fei::DofMapper<LocalOrdinal,GlobalOrdinal,DofOrder>& dof_mapper)
{
  LocalOrdinal rank = 0;
  const LocalOrdinal field = 0;
  const LocalOrdinal field_size = 2;
  GlobalOrdinal id = 0;

  dof_mapper.addDOF(rank, id, field);

  id = 1;
  dof_mapper.addDOF(rank, id, field);

  dof_mapper.setFieldSize(field, field_size);

  fei::set_dof_mappings(0, dof_mapper);
}

template<class LocalOrdinal,class GlobalOrdinal,class DofOrder>
void fill_dof_mapper_2(fei::DofMapper<LocalOrdinal,GlobalOrdinal,DofOrder>& dof_mapper)
{
  LocalOrdinal rank = 0;
  const LocalOrdinal field1 = 0;
  const LocalOrdinal field1_size = 2;
  const LocalOrdinal field2 = 1;
  const LocalOrdinal field2_size = 1;
  GlobalOrdinal id = 0;

  dof_mapper.addDOF(rank, id, field1);
  dof_mapper.addDOF(rank, id, field2);

  id = 1;
  dof_mapper.addDOF(rank, id, field1);
  dof_mapper.addDOF(rank, id, field2);

  dof_mapper.setFieldSize(field1, field1_size);
  dof_mapper.setFieldSize(field2, field2_size);

  fei::set_dof_mappings(0, dof_mapper);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(set_dof_mappings, test1, LocalOrdinal, GlobalOrdinal)
{
  fei::DofMapper<LocalOrdinal,GlobalOrdinal> dofmapper;

  TEUCHOS_TEST_EQUALITY(dofmapper.maps_are_valid(), false, out, success);

  fill_dof_mapper_1(dofmapper);

  LocalOrdinal rank = 0;
  GlobalOrdinal id = 1;
  LocalOrdinal field = 0;
  GlobalOrdinal idx = dofmapper.getGlobalIndex(rank, id, field);
  TEUCHOS_TEST_EQUALITY(idx, 2, out, success);

  idx = 1;
  std::pair<const fei::Dof<LocalOrdinal,GlobalOrdinal>*,GlobalOrdinal> result=
    dofmapper.getDof(idx);
  TEUCHOS_TEST_EQUALITY(result.first->id(), 0, out, success);
  TEUCHOS_TEST_EQUALITY(result.second, 1, out, success);

  idx = 2;
  result = dofmapper.getDof(idx);
  TEUCHOS_TEST_EQUALITY(result.first->id(), 1, out, success);
  TEUCHOS_TEST_EQUALITY(result.second, 0, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(set_dof_mappings, test2, LocalOrdinal, GlobalOrdinal)
{
  fei::DofMapper<LocalOrdinal,GlobalOrdinal> dofmapper;

  TEUCHOS_TEST_EQUALITY(dofmapper.maps_are_valid(), false, out, success);

  fill_dof_mapper_1(dofmapper);

  LocalOrdinal rank = 0;
  GlobalOrdinal id = 8;
  LocalOrdinal field = 0;
  TEUCHOS_TEST_THROW(dofmapper.getGlobalIndex(rank, id, field), std::runtime_error, out, success);

  GlobalOrdinal idx = 0;
  std::pair<const fei::Dof<LocalOrdinal,GlobalOrdinal>*,GlobalOrdinal> result=
    dofmapper.getDof(idx);
  TEUCHOS_TEST_EQUALITY(result.first->id(), 0, out, success);
  TEUCHOS_TEST_EQUALITY(result.second, 0, out, success);

  idx = 5;
  TEUCHOS_TEST_THROW(dofmapper.getDof(idx), std::runtime_error, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(set_dof_mappings, test3, LocalOrdinal, GlobalOrdinal)
{
  fei::DofMapper<LocalOrdinal,GlobalOrdinal,fei::less_field_rank_id<LocalOrdinal,GlobalOrdinal> > dofmapper;

  TEUCHOS_TEST_EQUALITY(dofmapper.maps_are_valid(), false, out, success);

  fill_dof_mapper_2(dofmapper);

  LocalOrdinal rank = 0;
  GlobalOrdinal id = 0;
  LocalOrdinal field1 = 0;
  LocalOrdinal field2 = 1;
  GlobalOrdinal idx = dofmapper.getGlobalIndex(rank, id, field2);
  TEUCHOS_TEST_EQUALITY(idx, 4, out, success);

  id = 1;
  idx = dofmapper.getGlobalIndex(rank, id, field2);
  TEUCHOS_TEST_EQUALITY(idx, 5, out, success);

  idx = 1;
  std::pair<const fei::Dof<LocalOrdinal,GlobalOrdinal>*,GlobalOrdinal> result=
    dofmapper.getDof(idx);
  TEUCHOS_TEST_EQUALITY(result.first->id(), 0, out, success);
  TEUCHOS_TEST_EQUALITY(result.first->field(), field1, out, success);
  TEUCHOS_TEST_EQUALITY(result.second, 1, out, success);

  idx = 2;
  result = dofmapper.getDof(idx);
  TEUCHOS_TEST_EQUALITY(result.first->id(), 1, out, success);
  TEUCHOS_TEST_EQUALITY(result.first->field(), field1, out, success);
  TEUCHOS_TEST_EQUALITY(result.second, 0, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(set_dof_mappings, test4, LocalOrdinal, GlobalOrdinal)
{
  fei::DofMapper<LocalOrdinal,GlobalOrdinal,fei::less_field_rank_id<LocalOrdinal,GlobalOrdinal> > dofmapper;

  TEUCHOS_TEST_EQUALITY(dofmapper.maps_are_valid(), false, out, success);

  fill_dof_mapper_2(dofmapper);

  LocalOrdinal rank = 0;
  GlobalOrdinal id = 8;
  LocalOrdinal field1 = 0;
//  LocalOrdinal field2 = 1;
  TEUCHOS_TEST_THROW(dofmapper.getGlobalIndex(rank, id, field1), std::runtime_error, out, success);

  GlobalOrdinal idx = 0;
  std::pair<const fei::Dof<LocalOrdinal,GlobalOrdinal>*,GlobalOrdinal> result=
    dofmapper.getDof(idx);
  TEUCHOS_TEST_EQUALITY(result.first->id(), 0, out, success);
  TEUCHOS_TEST_EQUALITY(result.first->field(), field1, out, success);
  TEUCHOS_TEST_EQUALITY(result.second, 0, out, success);

  idx = 6;
  TEUCHOS_TEST_THROW(dofmapper.getDof(idx), std::runtime_error, out, success);
}

#define UNIT_TEST_GROUP(LO,GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(set_dof_mappings, test1, LO, GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(set_dof_mappings, test2, LO, GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(set_dof_mappings, test3, LO, GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(set_dof_mappings, test4, LO, GO)

UNIT_TEST_GROUP(int,int)

}//namespace <anonymous>

