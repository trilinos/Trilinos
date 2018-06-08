ADD_SUBDIRECTORY(coffey)
ADD_SUBDIRECTORY(glaze)
ADD_SUBDIRECTORY(glaze2)
ADD_SUBDIRECTORY(rectilinearRCB)
ADD_SUBDIRECTORY(ordering)

# command that copies files during build
TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_input_files_for_driver
    SOURCE_FILES
        input_template.xml
        rcbTest.xml
        multiJaggedTest.xml
        GeomGenWeight.txt
        chacoGraphMetricsTest.xml
        chacoSimple3dTest.xml
        chacoSimpleTest.xml
        chacoVwgtTest.xml    
        chacoVwgt2Test.xml    
        chacoEdgeWeightTest.xml
        geometricTest.xml
        geometricVwgtTest.xml
        geometricVwgt2Test.xml
        Rectilinear.pmgen
        Rectilinear2D.pmgen
        Spherical.pmgen
        Poisson.pmgen
        pamgenGraphMetricsTest.xml
        pamgenRectilinearTest2D.xml 
        pamgenRectilinearTest3D.xml 
        pamgenSphericalTest.xml 
        pamgenPoissonTest.xml 
        parmetisSubComm.xml
        galeriTest.xml 
        simple.pmgen
        simplePamgen.xml
        simple3d.pmgen
        simple3dPamgen.xml
        pulpTest.xml
        pulpVwgtTest.xml
        pulpTestNoDistrib.xml
        pulpVwgtTestNoDistrib.xml
        test578.xml
        zoltanGraphAdapter.xml
        zoltanGraphAdapterVwgt.xml
        zoltanGraphAdapterEwgt.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

