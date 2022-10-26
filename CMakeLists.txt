
cmake_minimum_required(VERSION          3.10)
cmake_policy(SET CMP0048                NEW) 

option(BUILD_EXAMPLE "Build demo"        ON)
set(GK_VERSION_MAJOR                      3)
set(GK_VERSION_MINOR                      0)
set(GK_VERSION_PATCH                      0)

project(lib_opengjk 
    LANGUAGES C
	DESCRIPTION "openGJK library community edition (CE)"
    HOMEPAGE_URL "www.mattiamontanari.com/opengjk/"
    VERSION ${GK_VERSION_MAJOR}.${GK_VERSION_MINOR}.${GK_VERSION_PATCH}
)

set( C_STANDARD                          11)
set( CMAKE_CXX_VISIBILITY_PRESET         hidden)
set( CMAKE_VISIBILITY_INLINES_HIDDEN     1)
set( CMAKE_POLICY_DEFAULT_CMP0079        NEW)
set( CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS    ON)
set( CMAKE_INCLUDE_CURRENT_DIR           TRUE)

# Set compile flags
message(STATUS "Compiler in use: ${CMAKE_C_COMPILER_ID}")
set( CMAKE_C_FLAGS_DEBUG                 "-O0 -g -Wall")
set( CMAKE_C_FLAGS_RELEASE               "-O3 -Werror")




message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_CURRENT_SOURCE_DIR: ${PROJECT_CURRENT_SOURCE_DIR}")
message(STATUS "PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "PROJECT_CURRENT_BINARY_DIR: ${PROJECT_CURRENT_BINARY_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")



if (CMAKE_C_COMPILER_ID STREQUAL "Clang")
    # To use Clang, pass to cmake: -DCMAKE_C_COMPILER="clang"
elseif (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    # To use GNU, pass to cmake: -DCMAKE_C_COMPILER="gcc"
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pedantic -Wunused-macros -Wno-unused-command-line-argument")
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
elseif (CMAKE_C_COMPILER_ID STREQUAL "MSVC")
endif()

set( GK_PUBLIC_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(${PROJECT_NAME} 
    STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/openGJK.c 
    ${GK_PUBLIC_HEADER_DIR}/openGJK/openGJK.h
)

target_link_libraries(${PROJECT_NAME} PRIVATE m)

target_include_directories( 
    ${PROJECT_NAME} 
    PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}
)

set_target_properties(${PROJECT_NAME} 
    PROPERTIES 
        PUBLIC_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/include/openGJK/openGJK.h
)

if(BUILD_EXAMPLE)
    add_subdirectory(examples/c)
endif(BUILD_EXAMPLE)

# if (UNIX)
#     install(TARGETS ${PROJECT_NAME} PERMISSIONS WORLD_WRITE )
#     set(DESTDIR "/usr")
#     INSTALL(TARGETS ${PROJECT_NAME} 
#         LIBRARY DESTINATION "${DESTDIR}/lib"
#         PUBLIC_HEADER DESTINATION "${DESTDIR}/include"
#     )
# endif (UNIX)

# Wrap up feedback on setup
message(STATUS "Version     : " ${CMAKE_PROJECT_VERSION}  )
message(STATUS "Build type  : " ${CMAKE_BUILD_TYPE} )
