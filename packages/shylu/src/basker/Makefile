CC=g++

basker: basker_driver.o basker.o
	$(CC)  -o basker.exe basker_driver.o basker.o -lm
basker.o: basker.cpp basker.hpp basker_decl.hpp basker_def.hpp
	$(CC)  -c basker.cpp
basker_driver: basker_driver.cpp basker.hpp basker_decl.hpp basker_def.hpp
	$(CC)  -c basker_driver.cpp


