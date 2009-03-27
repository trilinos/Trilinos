.PHONY: clean cleandoc subdirs all doc

include ../make.include

all: subdirs $(PB_DIR)/lib/libpb.a

$(PB_DIR)/lib/libpb.a:
	ar rcs $(PB_DIR)/lib/libpb.a $(PB_DIR)/objs/*.o

subdirs:
	cp ../make.include .
	make -C src/
	ar rcs $(PB_DIR)/lib/libpb.a $(PB_DIR)/objs/*.o
	make -C tests/

linecount:
	find ./src -name \*.h -or -name \*.cpp -or -name *.hpp | xargs wc -l

doc:
	$(DOXYGEN) doxygen.cfg

clean:
	cp ../make.include .
	make -C tests/ clean
	make -C src/ clean
	-rm -f $(PB_DIR)/lib/libpb.a
	-rm -f make.include
	-rm -f objs/*.o
	-rm -f lib/*.a

cleandoc:
	-rm -fr doc/html
