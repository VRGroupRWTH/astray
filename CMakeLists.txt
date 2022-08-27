##################################################    Project     ##################################################
cmake_minimum_required(VERSION 3.23 FATAL_ERROR)
project               (astray VERSION 1.0 LANGUAGES C CXX)
list                  (APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set_property          (GLOBAL PROPERTY USE_FOLDERS ON)
set                   (CMAKE_CXX_STANDARD 17)
set                   (CMAKE_BUILD_TYPE Release)

include               (set_max_warning_level)
set_max_warning_level ()

##################################################    Options     ##################################################
option      (ASTRAY_BUILD_APPLICATIONS "Build the example applications (default: OFF)."  OFF)
option      (ASTRAY_BUILD_TESTS        "Build the tests (default: OFF)."                 OFF)
option      (ASTRAY_USE_FFMPEG         "Use FFMPEG for video output (default: OFF)."     OFF)
option      (ASTRAY_USE_MPI            "Use MPI for distributed compute (default: OFF)." OFF)

set         (ASTRAY_DEVICE_SYSTEM OMP CACHE STRING "Device system (default: OMP).")
set_property(CACHE ASTRAY_DEVICE_SYSTEM PROPERTY STRINGS CPP CUDA OMP TBB)

##################################################    Sources     ##################################################
file(GLOB_RECURSE PROJECT_HEADERS include/*.h include/*.hpp)
file(GLOB_RECURSE PROJECT_CMAKE_UTILS cmake/*.cmake)
file(GLOB_RECURSE PROJECT_MISC *.md *.txt)
set (PROJECT_FILES ${PROJECT_HEADERS} ${PROJECT_CMAKE_UTILS} ${PROJECT_MISC})

include            (assign_source_group)
assign_source_group(${PROJECT_FILES})

##################################################  Dependencies  ##################################################
include(import_library)

if               (ASTRAY_BUILD_TESTS)
  find_package   (doctest CONFIG REQUIRED)
  list           (APPEND PROJECT_LIBRARIES doctest::doctest)
  set            (PROJECT_DEPENDENCY_DOCTEST "find_dependency(doctest)")
endif            ()

if               (ASTRAY_USE_FFMPEG)
  find_package   (FFMPEG REQUIRED)
  import_library (FFMPEG_INCLUDE_DIRS FFMPEG_LIBRARIES)
  list           (APPEND PROJECT_COMPILE_DEFINITIONS -DASTRAY_USE_FFMPEG)
  set            (PROJECT_DEPENDENCY_FFMPEG "find_dependency(FFMPEG)")
endif            ()

if               (ASTRAY_USE_MPI)
  find_package   (MPI REQUIRED COMPONENTS C)
  list           (APPEND PROJECT_LIBRARIES MPI::MPI_C)
  list           (APPEND PROJECT_COMPILE_DEFINITIONS -DASTRAY_USE_MPI)
  set            (PROJECT_DEPENDENCY_MPI "find_dependency(MPI)")
endif            ()

if               (${ASTRAY_DEVICE_SYSTEM} STREQUAL CUDA)
  enable_language(CUDA)
  find_package   (CUDA REQUIRED)
  set            (CMAKE_CUDA_ARCHITECTURES "all-major")
  list           (APPEND PROJECT_COMPILE_OPTIONS $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda --expt-relaxed-constexpr --use_fast_math>)
  import_library (CUDA_INCLUDE_DIRS CUDA_LIBRARIES)
  list           (APPEND PROJECT_COMPILE_DEFINITIONS -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA)
  set            (PROJECT_DEPENDENCY_DEVICE "find_dependency(CUDA)")
elseif           (${ASTRAY_DEVICE_SYSTEM} STREQUAL OMP)
  find_package   (OpenMP REQUIRED)
  list           (APPEND PROJECT_LIBRARIES OpenMP::OpenMP_CXX)
  list           (APPEND PROJECT_COMPILE_DEFINITIONS -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP )
  set            (PROJECT_DEPENDENCY_DEVICE "find_dependency(OpenMP)")
elseif           (${ASTRAY_DEVICE_SYSTEM} STREQUAL TBB)
  find_package   (TBB REQUIRED)
  list           (APPEND PROJECT_LIBRARIES TBB::tbb TBB::tbbmalloc)
  list           (APPEND PROJECT_COMPILE_DEFINITIONS -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_TBB )
  set            (PROJECT_DEPENDENCY_DEVICE "find_dependency(TBB)")
elseif           (${ASTRAY_DEVICE_SYSTEM} STREQUAL CPP)
  list           (APPEND PROJECT_COMPILE_DEFINITIONS -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP )
endif            ()

find_package     (Eigen3 REQUIRED NO_MODULE)
list             (APPEND PROJECT_LIBRARIES Eigen3::Eigen)

find_package     (thrust CONFIG REQUIRED)
list             (APPEND PROJECT_LIBRARIES thrust::thrust)

##################################################    Targets     ##################################################
add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)
target_include_directories(${PROJECT_NAME} INTERFACE ${PROJECT_INCLUDE_DIRS})
target_link_libraries     (${PROJECT_NAME} INTERFACE ${PROJECT_LIBRARIES})
target_compile_definitions(${PROJECT_NAME} INTERFACE ${PROJECT_COMPILE_DEFINITIONS})
target_compile_options    (${PROJECT_NAME} INTERFACE ${PROJECT_COMPILE_OPTIONS})

# Hack for header-only project to appear in the IDEs.
add_library(${PROJECT_NAME}_ STATIC ${PROJECT_FILES})
target_include_directories(${PROJECT_NAME}_ PUBLIC 
  ${CMAKE_CURRENT_SOURCE_DIR}/include 
  ${CMAKE_CURRENT_BINARY_DIR})     
target_include_directories(${PROJECT_NAME}_ PUBLIC ${PROJECT_INCLUDE_DIRS})
target_link_libraries     (${PROJECT_NAME}_ PUBLIC ${PROJECT_LIBRARIES})
target_compile_definitions(${PROJECT_NAME}_ PUBLIC ${PROJECT_COMPILE_DEFINITIONS})
target_compile_options    (${PROJECT_NAME}_ PUBLIC ${PROJECT_COMPILE_OPTIONS})
set_target_properties     (${PROJECT_NAME}_ PROPERTIES LINKER_LANGUAGE CXX)

##################################################  Applications  ##################################################
if(ASTRAY_BUILD_APPLICATIONS)
  file(GLOB PROJECT_APPLICATION_HPPS applications/*.hpp)
  file(GLOB PROJECT_APPLICATION_CPPS applications/*.cu )

  # Hack for compiling .cu files as .cpp.
  if   (NOT ${ASTRAY_DEVICE_SYSTEM} STREQUAL CUDA)
    set_source_files_properties(${PROJECT_APPLICATION_CPPS} PROPERTIES LANGUAGE CXX)
  endif()

  foreach(_SOURCE ${PROJECT_APPLICATION_CPPS})
    get_filename_component(_NAME ${_SOURCE} NAME_WE)
    add_executable        (${_NAME} ${PROJECT_APPLICATION_HPPS} ${_SOURCE})
    target_link_libraries (${_NAME} ${PROJECT_NAME})
    target_link_libraries (${_NAME} ${PROJECT_LIBRARIES})
    set_property          (TARGET ${_NAME} PROPERTY FOLDER applications)
    set_property          (TARGET ${_NAME} PROPERTY LINKER_LANGUAGE CXX)
    assign_source_group   (${PROJECT_APPLICATION_HPPS} ${_SOURCE})
  endforeach()
endif()

##################################################    Testing     ##################################################
if(ASTRAY_BUILD_TESTS)
  enable_testing       ()
  set                  (TEST_MAIN_NAME test_main)
  set                  (TEST_MAIN_SOURCES tests/internal/main.cpp)
  add_library          (${TEST_MAIN_NAME} OBJECT ${TEST_MAIN_SOURCES})
  target_link_libraries(${TEST_MAIN_NAME} ${PROJECT_NAME})
  set_property         (TARGET ${TEST_MAIN_NAME} PROPERTY FOLDER tests)
  assign_source_group  (${TEST_MAIN_SOURCES})

  file(GLOB_RECURSE PROJECT_TEST_CPPS tests/*.cu)

  # Hack for compiling .cu files as .cpp.
  if   (NOT ${ASTRAY_DEVICE_SYSTEM} STREQUAL CUDA)
    set_source_files_properties(${PROJECT_TEST_CPPS} PROPERTIES LANGUAGE CXX)
  endif()

  foreach(_SOURCE ${PROJECT_TEST_CPPS})
    get_filename_component(_NAME      ${_SOURCE}    NAME_WE  )
    get_filename_component(_DIRECTORY ${_SOURCE}    DIRECTORY)
    get_filename_component(_PARENT    ${_DIRECTORY} NAME     )
    add_executable        (${_NAME} ${_SOURCE} $<TARGET_OBJECTS:${TEST_MAIN_NAME}>)
    target_link_libraries (${_NAME} ${PROJECT_NAME})
    target_link_libraries (${_NAME} ${PROJECT_LIBRARIES})
    add_test              (${_NAME} ${_NAME})
    set_property          (TARGET ${_NAME} PROPERTY FOLDER tests/${_PARENT})
    assign_source_group   (${_SOURCE})
  endforeach()
endif()

##################################################  Installation  ##################################################
install(DIRECTORY include/ DESTINATION include)
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME})
install(EXPORT  ${PROJECT_NAME} FILE   ${PROJECT_NAME}-targets.cmake NAMESPACE ${PROJECT_NAME}:: DESTINATION cmake)
export (TARGETS ${PROJECT_NAME} FILE   ${PROJECT_NAME}-targets.cmake)

include(CMakePackageConfigHelpers)
configure_package_config_file   (
  ${PROJECT_NAME}-config.cmake.in  
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake INSTALL_DESTINATION cmake)
write_basic_package_version_file(
  ${PROJECT_NAME}-config-version.cmake 
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake
  DESTINATION cmake)
