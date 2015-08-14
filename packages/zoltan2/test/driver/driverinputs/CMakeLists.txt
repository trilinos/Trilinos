# command that copies files during build
TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_input_files_for_driver
    SOURCE_FILES
        input_071315.xml
        input_072315.xml
        rcbTest.xml
        multiJaggedTest.xml
        GeomGenWeight.txt
        chacoSimple3dTest.xml
        chacoSimpleTest.xml
        chacoVwgt2Test.xml    
        chacoEdgeWeightTest.xml
        Rectilinear.pmgen
        Rectilinear2D.pmgen
        Spherical.pmgen
        pamgenRectilinearTest.xml 
        pamgenSphericalTest.xml 
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

