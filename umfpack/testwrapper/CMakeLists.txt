SET(SRCS1 test_di_demo.c umfpack_di_demo.c)
set_source_files_properties(umfpack_di_demo.c PROPERTIES COMPILE_DEFINITIONS "main=umfpack_di_demo")

ADD_EXECUTABLE(test1 ${SRCS1})
target_include_directories(test1 PRIVATE ../AMD/Include ../UMFPACK/Include ../UFconfig ../blaswrapper)
target_link_libraries(test1 PRIVATE umfpack)

#In production
#set( CMAKE_SKIP_BUILD_RPATH true )
SET(SRCS2 test_blas.c)
ADD_EXECUTABLE(test2 ${SRCS2})
target_include_directories(test2 PRIVATE ../blaswrapper)
target_link_libraries(test2 PRIVATE umfpack)

