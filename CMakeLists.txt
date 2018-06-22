cmake_minimum_required(VERSION 2.8)

project(meshDiff)
include_directories(src/lib)

# Blergh, I can't do the normal find_package(PythonLibs REQUIRED)
# b/c I have to support 32 and 64 bit builds with strange paths
# start annoying hard-coding
set(PY_VERSION 27_64 CACHE STRING "The python version as part of the path")
set(PY_PATH C:/Python${PY_VERSION})
string(REPLACE "_" ";" _PY_VLIST ${PY_VERSION})
list(GET _PY_VLIST 0 PY_SIMPLE_VERSION)
set(PYTHON_INCLUDE_DIRS ${PY_PATH}/include)
set(PYTHON_LIBRARIES ${PY_PATH}/libs/python${PY_SIMPLE_VERSION}.lib)
# End annoying hard-coding

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/../output/Python${PY_VERSION})
set(SOURCE_FILES
    "src/lib/mesh.cpp"
    "src/python/meshDiff.cpp"
)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(${PYTHON_INCLUDE_DIRS})
add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES})
set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".pyd")

target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES})
message("Python Libraries: ${PYTHON_LIBRARIES}")

install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION lib
    ARCHIVE DESTINATION lib
)

