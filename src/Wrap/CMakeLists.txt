find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
find_package(Boost COMPONENTS python REQUIRED)

# Without this, any build libraries automatically have names "lib{x}.so"
set(CMAKE_SHARED_MODULE_PREFIX "")


python3_add_library(cpyshapeit MODULE cpyshapeit.cpp)
target_link_libraries(cpyshapeit PRIVATE shapeit_lib Boost::python)
set_target_properties(cpyshapeit PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../pyshapeit)
# add_pytest(pyDeprotect ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py)
