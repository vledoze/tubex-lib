# ==================================================================
#  tubex-lib - cmake configuration file
# ==================================================================

# Installation process handled by setup.py

  set(PYTHON_PACKAGE_NAME    "pyibex_tubex")
  set(PYTHON_PACKAGE_DIR     "${CMAKE_CURRENT_BINARY_DIR}/python_package")
  file(MAKE_DIRECTORY        ${PYTHON_PACKAGE_DIR})

  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/${PYTHON_PACKAGE_NAME}/" "${PYTHON_PACKAGE_DIR}/${PYTHON_PACKAGE_NAME}")
  # execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/tests/" "${PYTHON_PACKAGE_DIR}/pyibex/tests")

  # find_package(pybind11 REQUIRED)

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in
    ${PYTHON_PACKAGE_DIR}/setup.py
  )

  set(pip_dir "${CMAKE_CURRENT_BINARY_DIR}/pip_package")

  #include_directories(.)

  pybind11_add_module(tubex SHARED
    pyibex_tubex.cpp
  )

  target_link_libraries(tubex PRIVATE ${LIBS} tubex-core) # tubex-graphics)

  add_custom_command(TARGET tubex POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:tubex>" "${PYTHON_PACKAGE_DIR}/${PYTHON_PACKAGE_NAME}"
  )

# Tests

  ######ADD_TEST(NAME python_module
  ######      COMMAND ${PYTHON_EXECUTABLE} -m unittest discover ${PYTHON_PACKAGE_NAME}.tests
  ######      WORKING_DIRECTORY ${PYTHON_PACKAGE_DIR}
  ######)

# Installation

  #install(DIRECTORY "${PYTHON_PACKAGE_DIR}/${PYTHON_PACKAGE_NAME}" DESTINATION "${PYTHON_INSTDIR}")

  #add_custom_target(pip_package)
  #add_custom_command(TARGET pip_package PRE_BUILD
  #    COMMAND ${PYTHON_EXECUTABLE} ARGS setup.py bdist_wheel -d ${CMAKE_BINARY_DIR}
  #    WORKING_DIRECTORY ${PYTHON_PACKAGE_DIR}
  #)

  add_custom_target(install_python)
  add_custom_command(TARGET install_python PRE_BUILD
      COMMAND pip install --upgrade .
      WORKING_DIRECTORY ${PYTHON_PACKAGE_DIR}
  )

  #add_dependencies(install_python pyibex)