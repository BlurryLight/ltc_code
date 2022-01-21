cmake_minimum_required(VERSION 3.15)
project(prefilter)
add_executable(prefilter prefilterAreaLight.cpp)
target_include_directories(prefilter PRIVATE ../external/glm)
target_include_directories(prefilter PRIVATE ../external/)
target_include_directories(prefilter PRIVATE ../external/CImg)
target_compile_definitions(prefilter PRIVATE cimg_display=0)

