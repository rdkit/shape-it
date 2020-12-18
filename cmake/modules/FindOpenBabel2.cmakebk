# Try to find OpenBabel2
# Once done this will define:
#
#    OPENBABEL2_FOUND - system has OpenBabel2
#    OPENBABEL2_INCLUDE_DIRS - the OpenBabel2 include directory
#    OPENBABEL2_INCLUDE_DIR - idem
#    OPENBABEL2_LIBRARIES - Link these to use OpenBabel2
#    OPENBABEL2_LIBRARY - idem
#
# Copyright (c) 2011-2012 Silicos-it, a division of Imacosi bvba

if (OPENBABEL2_INCLUDE_DIRS AND OPENBABEL2_LIBRARIES)
	set (OPENBABEL2_FOUND TRUE)
else ()
 	find_path (
		OPENBABEL2_INCLUDE_DIRS 
		NAMES
			openbabel-2.0/openbabel/obconversion.h
		PATHS
			ENV BABEL_INCLUDEDIR
		)
	find_library (
		OPENBABEL2_LIBRARIES
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
	if (OPENBABEL2_INCLUDE_DIRS AND OPENBABEL2_LIBRARIES)
		set (OPENBABEL2_INCLUDE_DIRS ${OPENBABEL2_INCLUDE_DIRS}/openbabel-2.0)
		set (OPENBABEL2_FOUND TRUE)
	endif ()
endif ()

if (NOT OPENBABEL2_FOUND)
      message (FATAL_ERROR "Could NOT find OpenBabel 2.2 or later ")
else ()
	message ("OpenBabel include directories: ${OPENBABEL2_INCLUDE_DIRS}")
	message ("OpenBabel link library: ${OPENBABEL2_LIBRARIES}")
#	set (OPENBABEL2_INCLUDE_DIR ${OPENBABEL2_INCLUDE_DIRS})
#	set (OPENBABEL2_LIBRARY ${OPENBABEL2_LIBRARIES})
endif ()
