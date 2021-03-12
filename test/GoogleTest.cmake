# Download and unpack googletest at configure time
configure_file(${PROJECT_SOURCE_DIR}/test/CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

# Include file lists
include(test/FileList.cmake)

# Now simply link against gtest or gtest_main as needed. Eg
add_executable(GOMC_Test ${sources} ${headers} ${libHeaders} ${libSources} ${TestHeaders} ${TestSources})
target_link_libraries(GOMC_Test gtest_main)
add_test(NAME BasicTypesTest COMMAND BasicTypesTest)
add_test(NAME CircuitTester COMMAND DialaTest)
add_test(NAME MolLookupTest COMMAND CheckConsensusBeta)
add_test(NAME PSFParserTest COMMAND CheckProtAndWaterTest)
add_test(NAME EndianTest COMMAND TestBitSwap)
add_test(NAME ParallelTemperingTest COMMAND ParallelTemperingTest)
if(MPI_FOUND)
  target_link_libraries(GOMC_Test ${MPI_LIBRARIES})
endif()