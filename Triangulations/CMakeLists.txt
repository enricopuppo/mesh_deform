cmake_minimum_required(VERSION 3.2)
project(Triangulations)


add_executable(${PROJECT_NAME} main.cpp)

set_target_properties(Triangulations PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)
set(CINOLIB_USES_TRIANGLE ON)
set(cinolib_DIR /Users/enricopuppo/Documents/GitHub/lib/cinolib)
find_package(cinolib REQUIRED)
target_link_libraries(Triangulations cinolib)






