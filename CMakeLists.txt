cmake_minimum_required (VERSION 3.4)
project (dnapermut CXX)
 
# add seqan3 to search path
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/seqan3/build_system")
 
# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)
 
# build app with seqan3
add_executable (dnapermut dnapermut.cpp)
target_link_libraries (dnapermut seqan3::seqan3)

# Check support for C++20
set_property(TARGET dnapermut PROPERTY CXX_STANDARD 20)
set_property(TARGET dnapermut PROPERTY CXX_STANDARD_REQUIRED ON)
