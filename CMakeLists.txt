cmake_minimum_required (VERSION 2.6) 

# Project
project (shape-it)

set (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
find_package (OpenBabel3 REQUIRED)
# Version
set (VERSION_MAJOR 1)
set (VERSION_MINOR 0)
set (VERSION_PATCH 1)
file (WRITE  ${PROJECT_SOURCE_DIR}/include/config.h "#define SHAPEIT_VERSION ${VERSION_MAJOR}\n")
file (APPEND ${PROJECT_SOURCE_DIR}/include/config.h "#define SHAPEIT_RELEASE ${VERSION_MINOR}\n")
file (APPEND ${PROJECT_SOURCE_DIR}/include/config.h "#define SHAPEIT_SUBRELEASE ${VERSION_PATCH}\n")

# Release type
set (CMAKE_BUILD_TYPE Release)

# Source files
file (GLOB_RECURSE sourcefiles "${PROJECT_SOURCE_DIR}/src/*.cpp")

# Include files
include_directories (${PROJECT_SOURCE_DIR}/include)
include_directories (${OPENBABEL3_INCLUDE_DIRS})

# Executables
add_executable (${PROJECT_NAME} ${sourcefiles})
target_link_libraries (${PROJECT_NAME} ${OPENBABEL3_LIBRARIES})
install (TARGETS ${PROJECT_NAME} DESTINATION bin)
message ("Installation directory of executable: ${CMAKE_INSTALL_PREFIX}/bin")

# Testing
#enable_testing()
#include(${CMAKE_MODULE_PATH}/IncludeAllTests.cmake)

# Packaging
set (CPACK_PACKAGE_NAME ${PROJECT_NAME})
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Shape-it - a program for shape-matching molecules")
set (CPACK_PACKAGE_VENDOR "Silicos-it, a division of Imacosi bvba")
set (CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR})
set (CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR})
set (CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
set (CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/COPYING.LESSER")
set (CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README")
set (CPACK_RESOURCE_FILE_INSTALL "${PROJECT_SOURCE_DIR}/INSTALL")
set (CPACK_SOURCE_GENERATOR "TGZ")
set (CPACK_SOURCE_IGNORE_FILES 
     ".DS_Store";
     ".git";
     "${PROJECT_NAME}$";
	 "${PROJECT_NAME}-1.0.0.tar.gz$";
     "${PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.tar.gz$";
     "CMakeCache.txt$";
     "CMakeFiles";
     "CPackConfig.cmake$";
     "CPackSourceConfig.cmake$";
     "CTestTestfile.cmake$";
     "Makefile$";
     "Testing";
     "_CPack_Packages";
     "cmake_install.cmake$";
     "install_manifest.txt$"
     )
set (CPACK_SOURCE_PACKAGE_FILE_NAME 
    "${PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
include (CPack)
