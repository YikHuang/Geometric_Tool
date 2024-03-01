cmake_minimum_required(VERSION 3.1)
project(hw5)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_OPENGL             "Build target igl::opengl"         ON)
option(LIBIGL_GLFW               "Build target igl::glfw"           ON)
option(LIBIGL_IMGUI              "Build target igl::imgui"          ON)
option(LIBIGL_USE_STATIC_LIBRARY "Build LibIGL as a static library" OFF)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../libigl)
    message("libigl directory exists")
    set(FETCHCONTENT_SOURCE_DIR_LIBIGL "${CMAKE_CURRENT_SOURCE_DIR}/../libigl" CACHE PATH "Path to local libigl copy (to avoid re-cloning for each build directory)" FORCE)
else()
    message("libigl directory doesn't exist")
endif()
include(libigl)

# Build mode that compiles with optimizations but leaves in assertions.
set(CMAKE_CXX_FLAGS_RELWITHASSERT        "-O3" CACHE STRING "Flags used during RelWithAssert builds" FORCE)
set(CMAKE_C_FLAGS_RELWITHASSERT          "-O3" CACHE STRING "Flags used during RelWithAssert builds" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_RELWITHASSERT    "" CACHE STRING "Flags used during RelWithAssert builds" FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_RELWITHASSERT "" CACHE STRING "Flags used during RelWithAssert builds" FORCE)
set(CMAKE_MODULE_LINKER_FLAGS_RELWITHASSERT "" CACHE STRING "Flags used during RelWithAssert builds" FORCE)
set(CMAKE_STATIC_LINKER_FLAGS_RELWITHASSERT "" CACHE STRING "Flags used during RelWithAssert builds" FORCE)
mark_as_advanced(
    CMAKE_CXX_FLAGS_RELWITHASSERT
    CMAKE_C_FLAGS_RELWITHASSERT
    CMAKE_EXE_LINKER_FLAGS_RELWITHASSERT
    CMAKE_SHARED_LINKER_FLAGS_RELWITHASSERT
    CMAKE_MODULE_LINKER_FLAGS_RELWITHASSERT
    CMAKE_STATIC_LINKER_FLAGS_RELWITHASSERT
)

set(SHADER_PATH ${CMAKE_CURRENT_SOURCE_DIR}/shaders CACHE PATH "location of shader programs")
set(MESH_PATH ${CMAKE_CURRENT_SOURCE_DIR}/meshes CACHE PATH "location of example meshes")

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to RelWithAssert")
    set(CMAKE_BUILD_TYPE "RelWithAssert")
endif()

# Export compile flags (used for autocompletion of the C++ code)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Require C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Make sure warnings/errors are still colorized when using Ninja for building.
add_definitions(-fdiagnostics-color=always)

# Add M_PI
add_compile_definitions(_USE_MATH_DEFINES)

# Add your project files
if (NOT SOLUTION)
    add_library(ecs278 bezier.cc polyinterp.cc spline.cc spline_interpolation.cc tensor_product_surface.cc triangle_surf.cc)
else()
    add_library(ecs278 bezier_solution.cc polyinterp_solution.cc spline_solution.cc spline_interpolation_solution.cc tensor_product_surface_solution.cc triangle_surf_solution.cc)
endif()

target_link_libraries(ecs278 PUBLIC igl::core igl::glfw igl::imgui)
target_compile_definitions(ecs278 PUBLIC "-DMESH_PATH=\"${MESH_PATH}\"")

add_executable(${PROJECT_NAME} main.cc KnotViewer.cc UVDomainWidget.cc)
target_compile_definitions(${PROJECT_NAME} PRIVATE "-DSHADER_PATH=\"${SHADER_PATH}\"")
target_link_libraries(${PROJECT_NAME} ecs278)

add_executable(test_derivatives test_derivatives.cc)
target_link_libraries(test_derivatives ecs278)
