set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    set(COMPILE_DEFINITIONS -DDEBUG_SDPCOLORING -DDEBUG_ILPCOLORING  -DDEBUG_LPCOLORING)
else()
    set(COMPILE_DEFINITIONS "")
endif()

#set(COMPILE_DEFINITIONS -DDEBUG_LIWEI ${COMPILE_DEFINITIONS})

set(LIBS programoptions gdsdb gdsparser gzstream ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})
if(UNIX AND NOT APPLE)
    set(LIBS ${LIBS} rt)
endif(UNIX AND NOT APPLE)

if(GUROBI_FOUND)
    include_directories(
        ${GUROBI_INCLUDE_DIRS}
        )
    set(LIBS ${GUROBI_LIBRARIES} ${LIBS})
    set(COMPILE_DEFINITIONS "GUROBI=1" ${COMPILE_DEFINITIONS})
else(GUROBI_FOUND)
    set(COMPILE_DEFINITIONS GUROBI=0 ${COMPILE_DEFINITIONS})
endif(GUROBI_FOUND)
if(ENABLE_LEMONCBC)
    set(COMPILE_DEFINITIONS LEMONCBC=1 ${COMPILE_DEFINITIONS})
else(ENABLE_LEMONCBC)
    set(COMPILE_DEFINITIONS LEMONCBC=0 ${COMPILE_DEFINITIONS})
endif(ENABLE_LEMONCBC)
set(COMPILE_DEFINITIONS CSDP=1 ${COMPILE_DEFINITIONS})
set(LIBS ${LIBS} sdp openblas m ${CMAKE_THREAD_LIBS_INIT} gfortran)

file(GLOB SOURCES 
    DL_MPL.cpp
    DL_struct.cpp
    Enums.cpp
    GdsiiIO.cpp
    LayoutDB.cpp
    LayoutDBPolygon.cpp 
    LayoutDBRect.cpp
    Msg.cpp 
    Params.cpp
    RecoverHiddenVertex.cpp 
    SimpleMPL.cpp
    )

include_directories(${LIMBO_SOURCE_DIR} ${Boost_INCLUDE_DIRS})
add_library(${PROJECT_NAME}_lib STATIC ${SOURCES})
add_executable(${PROJECT_NAME} main.cpp)
if(OpenMP_CXX_FOUND)
    target_compile_options(${PROJECT_NAME}_lib PRIVATE ${OpenMP_CXX_FLAGS})
    target_compile_options(${PROJECT_NAME} PRIVATE ${OpenMP_CXX_FLAGS})
    set(LIBS ${LIBS} ${OpenMP_CXX_LIB_NAMES})
endif(OpenMP_CXX_FOUND)
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${PROJECT_NAME}_lib ${LIBS})
message("-- COMPILE_DEFINITIONS = ${COMPILE_DEFINITIONS}")
if(COMPILE_DEFINITIONS)
    target_compile_definitions(${PROJECT_NAME}_lib PRIVATE ${COMPILE_DEFINITIONS})
    target_compile_definitions(${PROJECT_NAME} PRIVATE ${COMPILE_DEFINITIONS})
endif(COMPILE_DEFINITIONS)
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
