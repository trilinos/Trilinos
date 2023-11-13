// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __IMAGE_SET_VALID_TESTS_H
#define __IMAGE_SET_VALID_TESTS_H

#include <string>

std::string ImagesetCamera = R"(
     begin catalyst
      begin camera fooCamera
      end camera
      begin imageset fooImageset
        camera = fooCamera
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentation = R"(
    begin catalyst
      begin representation fooRepresentation
      end representation
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        representation = fooRepresentation
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetOperation = R"(
    begin catalyst
      begin clip fooOperation
      end clip
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        operation = fooOperation
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetBasename = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        image basename = foo
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetBasedirectory = R"(
     begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        image basedirectory = foo
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetFormat = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        image format = jpg
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetColorSettings = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        image size = 800 450
        axes color = 0.0 1.0 1.0
        edge color = 1.0 1.0 0.0
        text color = 1.0 0.0 1.0
        background color = 0.0 0.0 0.0
        show edges = true
        show axes = true
        show time annotation = true
        color by scalar = VON_MISES
      end imageset
    end
    )";

std::string ImagesetDigitCount = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        image digit count = 6
        image size = 800 450
      end imageset
    end
    )";

std::string ImageFilenameDatetimeSimtime = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset1
        image basename = "is_datetime."
        camera = singletestcamera
        image digit count = 6
        filename use datetime = true
        image size = 800 450
      end imageset
      begin imageset fooImageset2
        image basename = "is_iso_datetime."
        camera = singletestcamera
        filename use datetime = true
        filename datetime convert to integer = false
        image size = 800 450
      end imageset
      begin imageset fooImageset3
        image basename = "is_simulation_time."
        camera = singletestcamera
        image digit count = 6
        filename use simulation time = true
        image size = 800 450
      end imageset
      begin imageset fooImageset4
        image basename = "is_datetime_and_simulation_time."
        camera = singletestcamera
        image digit count = 6
        filename use simulation time = true
        filename use datetime = true
        image size = 800 450
      end imageset
      begin imageset fooImageset5
        image basename = "is_datetime_and_simulation_time_no_call_count."
        camera = singletestcamera
        filename use simulation time = true
        filename use datetime = true
        filename use call count = false
        image size = 800 450
      end imageset
      begin imageset fooImageset6
        filename use simulation time = true
        filename use datetime = true
        filename use call count = true
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetDefault = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
      end imageset
    end
    )";

#endif /* __IMAGE_SET_VALID_TESTS_H */
