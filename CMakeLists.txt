cmake_minimum_required(VERSION 3.10)
project(asteroid)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  

file(GLOB SOFA_SOURCES "external/sofa/20231011/c/src/*.c")
file(GLOB SOFA_HEADERS "external/sofa/20231011/c/src/*.h")

add_executable(${PROJECT_NAME}
    src/data_reader.cpp
    src/ephemeris.cpp
    src/geometry.cpp
    src/rk4Integrator.cpp
    src/util.cpp
    src/observatories.cpp

    external/libephaccess/daf_reader.c
    external/libephaccess/daf_reader.h
    external/libephaccess/ephaccess.c
    external/libephaccess/ephaccess.h
    
    ${SOFA_SOURCES}
    ${SOFA_HEADERS}

    main.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC 
    include
    external/sofa/20231011/s/src
)
