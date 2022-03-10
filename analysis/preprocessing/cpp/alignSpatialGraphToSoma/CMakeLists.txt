cmake_minimum_required(VERSION 2.6)

PROJECT(Connectome)

# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

# Find VTK
INCLUDE(${CMAKE_ROOT}/Modules/FindVTK.cmake)
IF(USE_VTK_FILE)
  INCLUDE(${USE_VTK_FILE})
ELSE(USE_VTK_FILE)
  MESSAGE(FATAL_ERROR
          "Cannot build without VTK.  Please set VTK_DIR.")
ENDIF(USE_VTK_FILE)

#Find GSL
FIND_PACKAGE(GSL REQUIRED)

#CMake Optimizing
SET(CMAKE_BUILD_TYPE RELEASE )
SET(CMAKE_CXX_FLAGS_DEBUG "-g -wd1268")
SET(CMAKE_C_FLAGS_DEBUG "-g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -openmp")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -openmp")
SET(CMAKE_VERBOSE_MAKEFILE "TRUE")

#Sources
INCLUDE_DIRECTORIES(src ${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})

SET(HelperSources src/helper.h src/helper.cpp)
SET(AlignSpatialGraphToSomaSrc src/alignSpatialGraphToSoma.cpp)
SET(CommonSources src/typedefs.h src/basics.h src/basics.cpp src/amiraReader.h src/amiraReader.cpp)

ADD_EXECUTABLE(AlignSpatialGraphToSoma ${AlignSpatialGraphToSomaSrc} ${CommonSources} ${HelperSources})

TARGET_LINK_LIBRARIES(AlignSpatialGraphToSoma ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(AlignSpatialGraphToSoma PROPERTIES LINK_FLAGS "-openmp")