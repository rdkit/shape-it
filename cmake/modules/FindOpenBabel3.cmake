# Try to find OpenBabel2
# Once done this will define:
#
#    OPENBABEL3_FOUND - system has OpenBabel3
#    OPENBABEL3_INCLUDE_DIRS - the OpenBabel2 include directory
#    OPENBABEL3_INCLUDE_DIR - idem
#    OPENBABEL3_LIBRARIES - Link these to use OpenBabel3
#    OPENBABEL3_LIBRARY - idem
#
# Copyright (c) 2011-2012 Silicos-it, a division of Imacosi bvba

if (OPENBABEL3_INCLUDE_DIRS AND OPENBABEL3_LIBRARIES)
	set (OPENBABEL3_FOUND TRUE)
else ()
 	find_path (
		OPENBABEL3_INCLUDE_DIRS 
		NAMES
			openbabel3/openbabel/obconversion.h
		PATHS
			ENV BABEL_INCLUDEDIR
		)
	find_library (
		OPENBABEL3_LIBRARIES
		NAMES
			openbabel
		PATHS
			/usr/
			/usr/local/
			/usr/local/lib/
			/usr/local/openbabel/
			/usr/local/openbabel/lib/
			ENV BABEL_LIBDIR
		)
	if (OPENBABEL3_INCLUDE_DIRS AND OPENBABEL3_LIBRARIES)
		set (OPENBABEL3_INCLUDE_DIRS ${OPENBABEL3_INCLUDE_DIRS}/openbabel3)
		set (OPENBABEL3_FOUND TRUE)
	endif ()
endif ()

if (NOT OPENBABEL3_FOUND)
      message (FATAL_ERROR "Could NOT find 3penBabel 2.2 or later ")
else ()
	message ("OpenBabel include directories: ${OPENBABEL3_INCLUDE_DIRS}")
	message ("OpenBabel link library: ${OPENBABEL3_LIBRARIES}")
#	set (OPENBABEL3_INCLUDE_DIR ${OPENBABEL3_INCLUDE_DIRS})
#	set (OPENBABEL3_LIBRARY ${OPENBABEL3_LIBRARIES})
endif ()
