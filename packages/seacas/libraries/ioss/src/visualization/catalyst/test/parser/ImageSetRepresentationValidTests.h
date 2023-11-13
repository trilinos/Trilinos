// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __IMAGE_SET_REPRESENTATION_VALID_TESTS_H
#define __IMAGE_SET_REPRESENTATION_VALID_TESTS_H

#include <string>

std::string ImagesetRepresentationSurfaces = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        show surfaces = false
        show edges = true
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationEdges = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        show edges = true
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationBoundingBox = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        show surfaces = false
        show bounding box = true
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationColorByScalar = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        color by scalar = VON_MISES
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationColorByVectorMagnitude = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        color by vector magnitude = displ
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationColorByVectorComponent = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        color by vector component = displ_x
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationColorByTensorComponent = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        color by tensor component = stress_xy
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationColorBySolidColor = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        color by solid color = 1.0 0.5 0.0
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationColorByBlockId = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        image size = 800 450
        color by blockid
      end imageset
    end
    )";

std::string ImagesetRepresentationTimeAnnotation = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        show time annotation = true
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationColorLegend = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        show color legend = true
        color by scalar = VON_MISES
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationColorLegendRange = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        color legend range = 22 32
        color by scalar = VON_MISES
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationColorLegendUseCurrentDataRange = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        color legend use current data range
        color by tensor component = stress_yz
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationColorLegendUseCumulativeDataRange = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        color legend use cumulative data range
        color by tensor component = stress_yz
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationColorLegendMinimumRange = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        color legend minimum range = -5 15
        color by scalar = VON_MISES
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationColorLegendMaximumRange = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        color legend maximum range = -5 15
        color by scalar = VON_MISES
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationColorLegendPosition = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        color legend position = right 0.5
        color by scalar = VON_MISES
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationPresetColorScale = R"(
    begin catalyst
      begin representation BlueToRedRainbowInvertedRep
        color by tensor component = stress_xy
        preset color scale = Blue_to_Red_Rainbow
        invert color scale = true
      end representation
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset BlueToRedRainbowInverted
        camera = singletestcamera
        representation = BlueToRedRainbowInvertedRep
        image size = 800 450
        image basename = "BlueToRedRainbow."
      end imageset
    end
    )";

std::string ImagesetRepresentationTimeAnnotationPosition = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        time annotation position = top left 0.5
        show time annotation = true
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetMultiColorBy = R"(
    begin catalyst
     begin camera singletestcamera
       look direction = -1 -1 -1
     end
     begin imageset is1
       image size = 800 450
       camera = singletestcamera
       color by scalar = VON_MISES
       image basename = is1.
     end imageset
     begin imageset is2
       image size = 800 450
       camera = singletestcamera
       color by tensor component = stress_xy
       image basename = is2.
     end imageset
     begin imageset is3
       image size = 800 450
       camera = singletestcamera
       color by vector magnitude = displ
       image basename = is3.
     end imageset
     begin imageset is4
       image size = 800 450
       camera = singletestcamera
       color by solid color = 0.1 0.9 0.95
       image basename = is4.
     end imageset
     begin imageset is5
       image size = 800 450
       camera = singletestcamera
       image basename = is5.
     end imageset
    end
    )";

std::string ImagesetAxesLegendOnOffMultiImageset = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset is1
        image size = 800 450
        camera = singletestcamera
        color by scalar = VON_MISES
        image basename = is1.
        show axes = true
        show orientation axes = true
        show time annotation = true
        show edges = true
      end imageset
      begin imageset is2
        image size = 800 450
        camera = singletestcamera
        color by tensor component = stress_xy
        image basename = is2.
        show color legend = false
        show orientation axes = false
      end imageset
      begin imageset is3
        image size = 800 450
        camera = singletestcamera
        color by scalar = VON_MISES
        image basename = is3.
        show color legend = true
        show axes = true
      end imageset
      begin imageset is4
        image size = 800 450
        camera = singletestcamera
        image basename = is4.
        show axes = true
        show edges = true
      end imageset
      begin imageset is5
        image size = 800 450
        camera = singletestcamera
        color by solid color = 1.0 0.5 0.0
        image basename = is5.
      end imageset
      begin imageset is6
        image size = 800 450
        camera = singletestcamera
        image basename = is6.
        show axes = true
        show edges = true
        show surfaces = false
      end imageset
      begin imageset is7
        image size = 800 450
        camera = singletestcamera
        color by solid color = 0.0 1.0 1.0
        image basename = is7.
        show axes = true
      end imageset
    end
    )";

std::string ImagesetRepresentationAxes = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        show axes = true
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationAxisLabelName = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        show axes = true
        x axis label name = xAxis
        y axis label name = yAxis
        z axis label name = zAxis
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationAxisLabel = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        show axes = true
        show x axis label = false
        show y axis label = false
        show z axis label = false
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationAxisTicMarks = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        show axes = true
        show x axis tic marks = false
        show y axis tic marks = false
        show z axis tic marks = false
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationAxisMinorTicMarks = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        show axes = true
        show x axis minor tic marks = false
        show y axis minor tic marks = false
        show z axis minor tic marks = false
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationAxisTicMarkLabels = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        show axes = true
        show x axis tic mark labels = false
        show y axis tic mark labels = false
        show z axis tic mark labels = false
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetRepresentationOrientationAxes = R"(
    begin catalyst
     begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        show orientation axes = true
        image size = 800 450
      end imageset
      begin imageset fooImageset2
        camera = singletestcamera
        show orientation axes = false
        image size = 800 450
      end imageset
    end
    )";

#endif /* __IMAGE_SET_REPRESENTATION_VALID_TESTS_H */
