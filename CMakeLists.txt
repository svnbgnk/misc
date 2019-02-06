# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (misc CXX)
set (CUR_DIR "${CMAKE_CURRENT_LIST_DIR}/..")

if (NOT DEFINED ${CMAKE_BUILD_TYPE})
   set (CMAKE_BUILD_TYPE Release)
   message("*** CMAKE_BUILD_TYPE was not passed to CMake, therefore it was set to Release. ***")
endif()

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)
find_package (BZip2)

set (CMAKE_PREFIX_PATH "/home/sven/devel/seqan/util/cmake")
set (SEQAN_INCLUDE_PATH "/home/sven/devel/seqan/include")

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)
find_package (OpenMP REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
#set (CMAKE_CXX_FLAGS "-std=c++14 -O0 ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -I/home/sven/devel/sdsl-lite/include")
set (CMAKE_CXX_FLAGS "-g -std=c++17 -pthread -DSEQAN_DISABLE_VERSION_CHECK=YES -I ${CUR_DIR}/sdsl-lite/include ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
add_executable (sample sample.cpp)
target_link_libraries (sample ${SEQAN_LIBRARIES})
add_executable (modify modify.cpp)
target_link_libraries (modify ${SEQAN_LIBRARIES})
add_executable (split split.cpp)
target_link_libraries (split ${SEQAN_LIBRARIES})
add_executable (selectRandom selectRandom.cpp)
target_link_libraries (selectRandom ${SEQAN_LIBRARIES})
add_executable (merge merge.cpp)
target_link_libraries (merge ${SEQAN_LIBRARIES})
add_executable (testOSSAlignments testOSSAlignments.cpp)
target_link_libraries (testOSSAlignments ${SEQAN_LIBRARIES})



