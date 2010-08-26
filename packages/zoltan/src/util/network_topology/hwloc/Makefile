#
# Build a apps that uses hwloc library
#
HWLOC_DIR=/home/lriesen/projects/hierarchical_partitioning/hwloc
HWLOC_INCLUDE=$(HWLOC_DIR)/include
HWLOC_LIB=$(HWLOC_DIR)/lib/libhwloc.so

all: node_topology zoltan_get_topology

CC=mpicc
CFLAGS+=-I$(HWLOC_INCLUDE)

%: %.c
	$(CC) $(CFLAGS) -o $@ $< $(HWLOC_LIB)

clean:
	rm -f node_topology zoltan_get_topology
