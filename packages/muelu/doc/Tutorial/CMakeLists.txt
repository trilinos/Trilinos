ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(tex)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_tutorial_build_files SOURCE_FILES src/Challenge.cpp src/laplace2d.cpp src/MLParameterList.cpp)
