// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __IMAGE_SET_OPERATION_VALID_TESTS_H
#define __IMAGE_SET_OPERATION_VALID_TESTS_H

#include <string>

std::string ImagesetOperationThreshold1 = R"(
    begin catalyst
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        threshold = scalar VON_MISES keep between 28.5 30.5
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetOperationThreshold2 = R"(
    begin catalyst
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        threshold = vector magnitude displ keep between 0.003 0.008
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetOperationThreshold3 = R"(
    begin catalyst
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        threshold = vector component displ_x keep between 0.003 0.008
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetOperationThreshold4 = R"(
    begin catalyst
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        threshold = tensor component stress_xy keep between 1.0 1.7
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetOperationThreshold5 = R"(
    begin catalyst
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        threshold = scalar VON_MISES keep above 28.5
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetOperationThreshold6 = R"(
    begin catalyst
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        threshold = scalar VON_MISES keep below 28.5
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetOperationBoxClip1 = R"(
    begin catalyst
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        boxclip = center 4 1 -0.61 extents 6 1.5 2 keep inside
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetOperationClip1 = R"(
    begin catalyst
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        clip = point 0 1 2 normal 3 4 5
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetOperationSlice1 = R"(
    begin catalyst
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        slice = point 0 1 2 normal 3 4 5
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetOperationContour1 = R"(
    begin catalyst
      begin camera singletestcamera
        look direction = -1 -1 -1
        look at absolute point = 5.0 1.0 0.0
        look at absolute distance = 10.0
      end
      begin imageset fooImageset
        camera = singletestcamera
        contour = scalar VON_MISES value list 30.0
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetOperationContour2 = R"(
    begin catalyst
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        contour = vector magnitude displ value list 0.0001 0.0021 0.0041 0.0061 0.0081 0.011
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetOperationContour3 = R"(
    begin catalyst
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        image size = 800 450
        contour = vector component displ_x value list 0.0001 0.0021 0.0041 0.0061 0.0081 0.011
      end imageset
    end
    )";

std::string ImagesetOperationContour4 = R"(
    begin catalyst
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        contour = tensor component stress_xy value list 0.1 0.5 1.0 1.5 1.8
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetOperationContour5 = R"(
    begin catalyst
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        contour = vector component displ_x value list 0.0001 0.0011 0.0021 0.0031 0.0041 0.0051 0.0061 0.0071 0.0081 0.0091 0.0101
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetOperationContour6 = R"(
    begin catalyst
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        contour = vector component displ_x value sequence 0.0001 0.001 0.03
        image size = 800 450
      end imageset
    end
    )";

#endif /* __IMAGE_SET_OPERATION_VALID_TESTS_H */
