

SET(SRCS1 test_di_demo.c umfpack_di_demo.c)
set_source_files_properties(umfpack_di_demo.c PROPERTIES COMPILE_DEFINITIONS "main=umfpack_di_demo")
ADD_EXECUTABLE(test1 ${SRCS1})
target_include_directories(test1 PRIVATE ../AMD/Include ../UMFPACK/Include ../UFconfig ../blaswrapper)
target_link_libraries(test1 PRIVATE umfpack_lgpl )

SET(SRCS2 test_zi_demo.c umfpack_zi_demo.c)
set_source_files_properties(umfpack_zi_demo.c PROPERTIES COMPILE_DEFINITIONS "main=umfpack_zi_demo")
ADD_EXECUTABLE(test2 ${SRCS2})
target_include_directories(test2 PRIVATE ../AMD/Include ../UMFPACK/Include ../UFconfig ../blaswrapper)
target_link_libraries(test2 PRIVATE umfpack_lgpl )

SET(SRCS3 test_blas.c)
ADD_EXECUTABLE(test3 ${SRCS3})
target_include_directories(test3 PRIVATE ../blaswrapper)
target_link_libraries(test3 PRIVATE umfpack_lgpl )

FOREACH (I test1 test2 test3)
ADD_CUSTOM_COMMAND(TARGET ${I}
          POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${I}> ${CMAKE_BINARY_DIR}
)

IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
SET_TARGET_PROPERTIES(${I} PROPERTIES INSTALL_RPATH "\$ORIGIN/")
ENDIF (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
IF (APPLE)
SET_TARGET_PROPERTIES(${I} PROPERTIES INSTALL_RPATH "@executable_path/")
#SET_TARGET_PROPERTIES(${I} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup -flat_namespace")
ENDIF (APPLE)
ENDFOREACH ()
