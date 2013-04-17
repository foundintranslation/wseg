PROJECT (wseg)

SET(WSEG_VERSION_MAJOR "1")
SET(WSEG_VERSION_MINOR "0")
SET(WSEG_VERSION_PATCH "0")



# Some options

MATH(EXPR WSEG_INT_VERSION "(${WSEG_VERSION_MAJOR} * 10000) + (${WSEG_VERSION_MINOR} * 100) + (${WSEG_VERSION_PATCH} * 1)" )

SET(WSEG_VERSION "${WSEG_VERSION_MAJOR}.${WSEG_VERSION_MINOR}.${WSEG_VERSION_PATCH}")

#CMake 2.6+ is recommended to an improved Boost module
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0 FATAL_ERROR)

IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    SET(CMAKE_INSTALL_PREFIX
        "inst" CACHE PATH "WSEG install prefix" FORCE
        )
    MESSAGE(STATUS "You have not set the install dir, default to './inst', if
    you want to set it, use cmake -DCMAKE_INSTALL_PREFIX to do so")
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

#set various platform specific global options
IF(WIN32)
 SET(CMAKE_DEBUG_POSTFIX "d")
 ADD_DEFINITIONS("-DWIN32")
 ADD_DEFINITIONS("-DW32")
 ADD_DEFINITIONS("-D_WIN32")
ENDIF(WIN32)

INCLUDE(cmake/FindTR1.cmake)

# include specific modules
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
IF(WIN32)
	OPTION( ENABLE_MULTI_THREAD "Set to OFF to build single-thread executable" OFF )
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/vccomp)
ELSE()
	OPTION( ENABLE_MULTI_THREAD "Set to OFF to build single-thread executable" ON )
ENDIF()	
OPTION( WITH_SRILM "Set the path to SRILM" )
OPTION( SRILM_MACHINETYPE "SRILM Machine Type" "i686-m64")

ADD_DEFINITIONS("-D_FILE_OFFSET_BITS=64")
ADD_DEFINITIONS("-D_LARGE_FILES")


IF(ENABLE_MULTI_THREAD)
	IF(WIN32)
		ADD_DEFINITIONS("-DWITH_THREADS")
	ELSE(WIN32)
	    ADD_DEFINITIONS("-DWITH_THREADS")
		LINK_LIBRARIES(pthread)
	ENDIF(WIN32)
ENDIF()

IF (WITH_SRILM)
    ADD_DEFINITIONS("-DLM_SRI")
    IF (SRILM_MACHINETYPE)
        SET (SRILM_LIB_PATH ${WITH_SRILM}/lib/${SRILM_MACHINETYPE})
    ELSE()
        SET (SRILM_MACHINETYPE "i686-m64")
        SET (SRILM_LIB_PATH ${WITH_SRILM}/lib/${SRILM_MACHINETYPE})
    ENDIF()
    
    LINK_DIRECTORIES(${SRILM_LIB_PATH}) 
    LINK_LIBRARIES(oolm lattice flm dstruct misc )
ENDIF()


#IF (UNIX AND NOT APPLE)
#    LINK_LIBRARIES (rt)
#ENDIF()

# UNCOMMENT WHEN YOU WANT BOOST
# Find Boost library, specify
# BOOST_ROOT=/e/programs/boost_1_35_0
# BOOST_LIBRARYDIR=$BOOST_ROOT/stage/lib

#set(Boost_USE_STATIC_LIBS        ON)
#set(Boost_USE_MULTITHREADED      ON)
#set(Boost_USE_STATIC_RUNTIME    OFF)


#IF(WIN32)
#FIND_PACKAGE( Boost 1.46 COMPONENTS regex thread "program_options" filesystem
#    system unit_test_framework iostreams zlib )
#ELSE()
#FIND_PACKAGE( Boost 1.46 COMPONENTS regex thread "program_options" filesystem
#    system unit_test_framework iostreams )
#ENDIF(WIN32)

#IF(Boost_FOUND)
#    IF (NOT _boost_IN_CACHE)
#     MESSAGE( "Boost found" )
#     message(STATUS "Boost_INCLUDE_DIR    : ${Boost_INCLUDE_DIR}")
#    ENDIF (NOT _boost_IN_CACHE)#
#	LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS} )
#	INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )
#ELSE()
#    MESSAGE(FATAL_ERROR "Boost not found, please set the BOOST_ROOT and BOOST_LIBRARYDIR environment variables " )
#ENDIF()

# Configure and find ZLIB


#IF(WIN32)
    #FIND_PACKAGE(ZLIB)
    #IF(ZLIB_FOUND)
    #INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIRS})
	#don'e link ZLIB
    #SET (ZLIB_LIBRARIES )
    #ELSE()
    #MESSAGE(FATAL_ERROR "ZLIB not found" )
    #ENDIF (ZLIB_FOUND)
    #ELSE()
    #FIND_PACKAGE(ZLIB)

    #IF(ZLIB_FOUND)
    #    INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIRS})
    #LINK_DIRECTORIES( ${ZLIB_LIBRARY_DIRS})
    #ELSE()
    #MESSAGE(FATAL_ERROR "ZLIB not found" )
    #ENDIF (ZLIB_FOUND)
#ENDIF()

#FIND_PACKAGE(BZip2)

#IF( BZIP2_FOUND )
#    INCLUDE_DIRECTORIES(${BZIP2_INCLUDE_DIRS})
#    LINK_DIRECTORIES( ${BZIP2_LIBRARY_DIRS} )
#ELSE()
#    MESSAGE(FATAL_ERROR "BZIP2 not found" )
#ENDIF(BZIP2_FOUND)

#IF(WIN32)
#ADD_SUBDIRECTORY (vccomp)
#ENDIF(WIN32)

SET( WSEGSRC Array.cpp     Counts.h   File.h     LM.cpp       Map2.cpp  MemStats.h  NgramLM.cpp     option.h  stdafx.cpp    targetver.h   Trie.h     wseg.cpp
Array.h       Debug.cpp  GBmap.cpp  LM.h         Map2.h    MLSeg.cpp   NgramStats.cpp  Prob.cpp  stdafx.h      tls.h         tserror.h  XCount.h
Boolean.h     Debug.h    LHash.cpp  LMStats.cpp  Map.cpp   MLSeg.h     NgramStats.h    Prob.h    SubVocab.cpp  TLSWrapper.h  Vocab.cpp  zio.cpp
cfuncproto.h  File.cpp   LHash.h    LMStats.h    Map.h     Ngram.h     option.cpp      SArray.h  SubVocab.h    Trie.cpp      Vocab.h    zio.h
)

ADD_EXECUTABLE(wseg ${WSEGSRC} )
ADD_DEFINITIONS( -DINSTANTIATE_TEMPLATES -DNO_TLS )
TARGET_LINK_LIBRARIES (
    wseg
    m
    )

#set properties on the libraries
SET_TARGET_PROPERTIES(wseg PROPERTIES
    VERSION ${WSEG_VERSION}
    COMPILE_DEFINITIONS_DEBUG _DEBUG
)


INSTALL(TARGETS wseg 
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS "ENABLE_MULTI_THREAD = ${ENABLE_MULTI_THREAD}" )
MESSAGE( STATUS "WITH_SRILM = ${WITH_SRILM}" )
MESSAGE( STATUS "SRILM_MACHINETYPE = ${SRILM_MACHINETYPE}" )
MESSAGE( STATUS "Change a value with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )


SET( ENABLE_MULTI_THREAD "${ENABLE_MULTI_THREAD}" CACHE BOOL
    "Set to OFF to build single-thread executable" FORCE )
SET( WITH_SRILM "${WITH_SRILM}" CACHE PATH
    "Set to path to SRILM Build" FORCE )
SET( SRILM_MACHINETYPE "${SRILM_MACHINETYPE}" CACHE STRING
    "Set the build type of SRILM" FORCE )

