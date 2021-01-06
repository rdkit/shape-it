# FindRDKit.cmake
# Copyright (C) 2013-2017 NextMove Software
# Try to find RDKit headers and libraries
# Defines:
#
#  RDKIT_FOUND - system has RDKit
#  RDKIT_INCLUDE_DIR - the RDKit include directory
#  RDKIT_LIBRARIES - Link these to use RDKit

if(RDKIT_INCLUDE_DIR AND RDKIT_LIBRARIES)
  # in cache already or user-specified
  set(RDKIT_FOUND TRUE)

else()

  if(NOT RDKIT_INCLUDE_DIR)
    if(WIN32)
      find_path(RDKIT_INCLUDE_DIR GraphMol/RDKitBase.h
        PATHS
        ${RDKIT_DIR}\\Code
        ${RDKIT_DIR}\\include\\rdkit
        $ENV{RDKIT_INCLUDE_DIR}
        $ENV{RDKIT_INCLUDE_PATH}
        $ENV{RDKIT_BASE}\\Code
        $ENV{RDBASE}\\Code
        C:\\RDKit\\include
        C:\\RDKit\\Code
      )
    else()
      find_path(RDKIT_INCLUDE_DIR GraphMol/RDKitBase.h
        PATHS
          ${RDKIT_DIR}/Code
          ${RDKIT_DIR}/include/rdkit
          $ENV{RDKIT_INCLUDE_DIR}
          $ENV{RDKIT_INCLUDE_PATH}
          $ENV{RDKIT_BASE}/Code
          $ENV{RDBASE}/Code
          /usr/local/rdkit/include/Code
          /usr/local/rdkit/include
          /usr/local/rdkit/Code
          /usr/local/include/rdkit
          ~/rdkit/Code
      )
    endif()
    if(RDKIT_INCLUDE_DIR)
       message(STATUS "Found RDKit include files at ${RDKIT_INCLUDE_DIR}")
    endif()
  endif()

  if(NOT RDKIT_LIBRARIES)
    find_library(RDKIT_FILEPARSERS_LIB
      NAMES RDKitFileParsers_static RDKitFileParsers
            FileParsers_static FileParsers
      NO_DEFAULT_PATH
      PATHS
        ${RDKIT_DIR}/lib
        $ENV{RDKIT_LIB_DIR}
        $ENV{RDKIT_LIB_PATH}
        $ENV{RDKIT_LIBRARIES}
        $ENV{RDKIT_BASE}/lib
        $ENV{RDBASE}/lib
        /usr/local/rdkit/lib
        ~/rdkit/lib
        $ENV{LD_LIBRARY_PATH}
    )
    find_library(RDKIT_FILEPARSERS_LIB # repeat but use the default path
      NAMES RDKitFileParsers_static RDKitFileParsers
            FileParsers_static FileParsers
    )
    if(RDKIT_FILEPARSERS_LIB)
      GET_FILENAME_COMPONENT(RDKIT_LIBRARY_DIR ${RDKIT_FILEPARSERS_LIB} PATH)

      # Note that the order of the following libraries is significant!!
      find_library(SMILESPARSE_LIB NAMES RDKitSmilesParse_static
                                         RDKitSmilesParse
                                         SmilesParse_static
                                         SmilesParse
                                   HINTS ${RDKIT_LIBRARY_DIR} NO_DEFAULT_PATH)
      find_library(DEPICTOR_LIB NAMES RDKitDepictor_static
                                      RDKitDepictor
                                      Depictor_static
                                      Depictor
                                HINTS ${RDKIT_LIBRARY_DIR} NO_DEFAULT_PATH)
      find_library(SUBSTRUCTMATCH_LIB NAMES RDKitSubstructMatch_static
                                            RDKitSubstructMatch
                                            SubstructMatch_static
                                            SubstructMatch
                                HINTS ${RDKIT_LIBRARY_DIR} NO_DEFAULT_PATH)
      find_library(RINGDECOMPOSER_LIB NAMES RDKitRingDecomposerLib_static
                                      RDKitRingDecomposerLib
                                      RingDecomposerLib_static
                                      RingDecomposerLib
                                HINTS ${RDKIT_LIBRARY_DIR} NO_DEFAULT_PATH)
      if(NOT RINGDECOMPOSER_LIB)
         set(RINGDECOMPOSER_LIB "")
      endif()
      find_library(GRAPHMOL_LIB NAMES RDKitGraphMol_static
                                      RDKitGraphMol
                                      GraphMol_static
                                      GraphMol
                                HINTS ${RDKIT_LIBRARY_DIR} NO_DEFAULT_PATH)
      find_library(RDGEOMETRYLIB_LIB NAMES RDKitRDGeometryLib_static
                                           RDKitRDGeometryLib
                                           RDKitGeometryLib_static
                                           RDKitGeometryLib
                                           RDGeometryLib_static
                                           RDGeometryLib
                                     HINTS ${RDKIT_LIBRARY_DIR} NO_DEFAULT_PATH)
      find_library(RDGENERAL_LIB NAMES RDKitRDGeneral_static
                                       RDKitRDGeneral
                                       RDKitGeneral_static
                                       RDKitGeneral
                                       RDGeneral_static
                                       RDGeneral
                                 HINTS ${RDKIT_LIBRARY_DIR} NO_DEFAULT_PATH)
      find_library(COORDGEN_LIB NAMES RDKitcoordgen_static
                                      RDKitcoordgen
                                      RDKitcoordgenlib_static
                                      RDKitcoordgenlib
                                      coordgenlib_static
                                      coordgenlib
                                      coordgen_static
                                      coordgen
                                 HINTS ${RDKIT_LIBRARY_DIR} NO_DEFAULT_PATH)
      find_library(MAEPARSER_LIB NAMES RDKitmaeparser_static
                                       RDKitmaeparser
                                       maeparser_static
                                       maeparser
                                 HINTS ${RDKIT_LIBRARY_DIR} NO_DEFAULT_PATH)
      find_library(RDDATASTRUCTS_LIB NAMES RDKitDataStructs_static
                                           RDKitDataStructs
                                           DataStructs_static
                                           DataStructs
                                 HINTS ${RDKIT_LIBRARY_DIR} NO_DEFAULT_PATH)
      find_library(RDMOLTRANSFORMS_LIB NAMES RDKitMolTransforms_static
                                           RDKitMolTransforms
                                           MolTransforms_static
                                           MolTransforms
                                 HINTS ${RDKIT_LIBRARY_DIR} NO_DEFAULT_PATH)
      find_library(RDEIGENSOLVERS_LIB NAMES RDKitEigenSolvers_static
                                           RDKitEigenSolvers
                                           EigenSolvers_static
                                           EigenSolvers
                                 HINTS ${RDKIT_LIBRARY_DIR} NO_DEFAULT_PATH)
      set (RDKIT_LIBRARIES ${RDKIT_FILEPARSERS_LIB} ${SMILESPARSE_LIB}
        ${DEPICTOR_LIB} ${RDMOLTRANSFORMS_LIB} ${RDEIGENSOLVERS_LIB} ${SUBSTRUCTMATCH_LIB}  ${GRAPHMOL_LIB} ${RINGDECOMPOSER_LIB}
        ${RDDATASTRUCTS_LIB} ${RDGEOMETRYLIB_LIB} ${RDGENERAL_LIB})
      if(COORDGEN_LIB AND MAEPARSER_LIB)
        set(RDKIT_LIBRARIES ${RDKIT_LIBRARIES} ${COORDGEN_LIB} ${MAEPARSER_LIB})
      endif()
      message(STATUS "Found RDKit libraries at ${RDKIT_LIBRARY_DIR}")
    endif()
  endif()

  if(RDKIT_INCLUDE_DIR AND RDKIT_LIBRARIES)
    set(RDKIT_FOUND TRUE)
  endif()

  mark_as_advanced(RDKIT_INCLUDE_DIR RDKIT_LIBRARIES RDKIT_FILEPARSERS_LIB)
endif()

