#include <cgns/Iocgns_StructuredZoneData.h>
#include <vector>

void cleanup(std::vector<Iocgns::StructuredZoneData *> &zones);
void check_split_assign(std::vector<Iocgns::StructuredZoneData *> &zones,
                        double load_balance_tolerance, size_t proc_count, double min_toler = 0.9,
                        double max_toler = 1.0);
