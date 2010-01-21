.PHONY: clean cleandoc subdirs all doc

include ./make.include

all: subdirs $(Teko_DIR)/lib/libteko.a

$(Teko_DIR)/lib/libteko.a:
	mkdir -p lib
	ar rcs $(Teko_DIR)/lib/libteko.a $(Teko_DIR)/objs/*.o

subdirs:
	mkdir -p objs
	mkdir -p lib
	make -C src/
	ar rcs $(Teko_DIR)/lib/libteko.a $(Teko_DIR)/objs/*.o
	make -C tests/

linecount:
	find ./src -name \*.h -or -name \*.cpp -or -name *.hpp | xargs wc -l

total-linecount:
	find . -name \*.h -or -name \*.cpp -or -name *.hpp | xargs wc -l

doc:
	$(DOXYGEN) doxygen.cfg

clean:
	make -C tests/ clean
	make -C src/ clean
	-rm -f $(Teko_DIR)/lib/libteko.a
	-rm -f objs/*.o
	-rm -f lib/*.a

cleandoc:
	-rm -fr doc/html
