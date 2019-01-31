# ==================================================================
#  tubex-lib - cmake configuration file
# ==================================================================

# Settings:

  cmake_minimum_required(VERSION 2.8.7)

  set(CMAKE_CXX_STANDARD 11)
  set(PACKAGE_NAME       tubex-lib)
  set(PACKAGE_VERSION    "2.0")
  option(WITH_PYTHON     "build library with python bindings" ON)

  project(${PACKAGE_NAME})

  get_filename_component(TUBELIB_BIN_DIR      "./bin"     ABSOLUTE)
  get_filename_component(TUBELIB_LIB_DIR      "./lib"     ABSOLUTE)
  get_filename_component(TUBELIB_BUILD_DIR    "./make"    ABSOLUTE)
  get_filename_component(TUBELIB_SRC_DIR      "./build"   ABSOLUTE)
  get_filename_component(TUBELIB_SCRIPTS_DIR  "./scripts" ABSOLUTE)

  set(CMAKE_MODULE_PATH  "${TUBELIB_SCRIPTS_DIR}/cmake_modules/")

  #set(CMAKE_CXX_FLAGS ${CMAKE_CSS_FLAGS} "-DNDEBUG") # Comment for debug
  set(CMAKE_CXX_FLAGS ${CMAKE_CSS_FLAGS} "-fPIC")

# IBEX-lib package

  set(IbexLib_USE_STATIC 1) # force ibex static linking

  find_package(IbexLib)
    if(IBEX_FOUND)
      include_directories(${IBEX_INCLUDE_DIRS})
      set(LIBS ${LIBS} ${IBEX_LIBRARIES})
    else()
      message(FATAL_ERROR "Unable to find ibex-lib. You may need to set IBEX_ROOT ${IBEX_ROOT}")
  endif()

  ###set(LIBS ${LIBS} -libex -lm -lgaol -lgdtoa -lultim) # using library Gaol
  ##set(LIB S${LIBS} -lm -lprim -libex) # using library Filib

  #set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Sources

  include_directories(${TUBELIB_BUILD_DIR}/include/)
  add_subdirectory(src)

# Tests

  option(BUILD_TESTS "Build test" ON)
  if(BUILD_TESTS)
    enable_testing(true)
    add_subdirectory(tests)
    include(CTest)
  endif()

# Examples

  add_subdirectory(examples)

# Python bindings (pyibex)

  if(WITH_PYTHON)
    add_subdirectory(pybind11)
    add_subdirectory(python)
  endif(WITH_PYTHON)

# Uninstallation

  configure_file("${CMAKE_HOME_DIRECTORY}/scripts/cmake_modules/cmake_uninstall.cmake.in" "cmake_uninstall.cmake" IMMEDIATE @ONLY)
  add_custom_target(uninstall "${CMAKE_COMMAND}" -P "cmake_uninstall.cmake")