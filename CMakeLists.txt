cmake_minimum_required(VERSION 3.10)
project(asteroid)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  

add_executable(${PROJECT_NAME}
    src/reader_obs.cpp

    main.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC include)
