# ==================================================================
#  tubex-lib - cmake configuration file
# ==================================================================

# Installation process handled by setup.py
  
  set(PYTHON_PACKAGE_NAME    "pytubex")
  set(PYTHON_PACKAGE_DIR     "${CMAKE_CURRENT_BINARY_DIR}/python_package")
  file(MAKE_DIRECTORY        ${PYTHON_PACKAGE_DIR})

  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/${PYTHON_PACKAGE_NAME}/" "${PYTHON_PACKAGE_DIR}/${PYTHON_PACKAGE_NAME}")

  # find_package(pybind11)

  file(WRITE "${PYTHON_PACKAGE_DIR}/empty.c" "// empty source file")

  configure_file(
   ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in
   ${PYTHON_PACKAGE_DIR}/setup.py
  )

  include_directories(.)
  include_directories(../src/core/dynamics/tube/) # ?

  set(PYBIND11_CPP_STANDARD -std=c++11)
  pybind11_add_module(pytubex SHARED
    pyibex_tubex.cpp
  )

  target_link_libraries(pytubex PRIVATE
    tubex-core
  )

  target_include_directories(pytubex
    PUBLIC ${IBEX_INCLUDE_DIRS}
  )

  # Copie des fichiers .so dans le module python
  add_custom_command(TARGET pytubex POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:pytubex>" "${PYTHON_PACKAGE_DIR}/pytubex"
  )

  # ADD_TEST(NAME python_module
  #       COMMAND ${PYTHON_EXECUTABLE} -m unittest discover pytubex
  #       WORKING_DIRECTORY ${PYTHON_PACKAGE_DIR}
  # )

  install(DIRECTORY "${PYTHON_PACKAGE_DIR}/pytubex" DESTINATION "${PYTHON_INSTDIR}")

  add_custom_target(pip_package)
  add_custom_command(TARGET pip_package PRE_BUILD
    COMMAND ${PYTHON_EXECUTABLE} ARGS setup.py bdist_wheel -d ${CMAKE_BINARY_DIR}
    WORKING_DIRECTORY ${PYTHON_PACKAGE_DIR}
  )

  add_custom_target(install_python)
  add_custom_command(TARGET install_python PRE_BUILD
    COMMAND pip install --upgrade .
    WORKING_DIRECTORY ${PYTHON_PACKAGE_DIR}
  )

  add_dependencies(install_python pytubex)

  # --build-purelib

  # ------------------------------
  #add_custom_target(install_python)
  #add_custom_command(TARGET install_python PRE_BUILD
  #    COMMAND pip install --upgrade .
  #    WORKING_DIRECTORY ${PYTHON_PACKAGE_DIR}
  #)