CC=g++

basker: basker.o basker_driver.o
	$(CC)  -o basker.exe basker_driver.o basker.o -lm
basker.o: basker.cpp basker.hpp basker_decl.hpp basker_def.hpp
	$(CC)  -c basker.cpp
basker_driver: basker_driver.cpp basker.hpp basker_decl.hpp basker_def.hpp
	$(CC)  -c basker_driver.cpp
matlab: baskermex.cpp basker.cpp
	mex -v -largeArrayDims basker.cpp baskermex.cpp 
matlab_complex: baskermex_complex.cpp basker.cpp
	mex -v -largeArrayDims basker.cpp baskermex_complex.cpp



