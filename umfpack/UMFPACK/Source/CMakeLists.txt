INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/AMD/Include ${CMAKE_SOURCE_DIR}/UFconfig ../Include)
SET (GLOBAL_SRCS 
	umfpack_global.c
	umfpack_tictoc.c
	umfpack_timer.c
)

#need DINT and DZINT
SET (GROUP1_SRCS
umf_2by2.c
umf_analyze.c
umf_apply_order.c
umf_assemble.c
umf_blas3_update.c
umf_build_tuples.c
umf_colamd.c
umf_create_element.c
umf_dump.c
umf_extend_front.c
umf_free.c
umf_fsize.c
umf_garbage_collection.c
umf_get_memory.c
umf_grow_front.c
umf_init_front.c
umf_is_permutation.c
umf_kernel.c
umf_kernel_init.c
umf_kernel_wrapup.c
umf_local_search.c
umf_lsolve.c
umf_ltsolve.c
umf_malloc.c
umf_mem_alloc_element.c
umf_mem_alloc_head_block.c
umf_mem_alloc_tail_block.c
umf_mem_free_tail_block.c
umf_mem_init_memoryspace.c
umfpack_col_to_triplet.c
umfpack_defaults.c
umfpack_free_numeric.c
umfpack_free_symbolic.c
umfpack_get_determinant.c
umfpack_get_lunz.c
umfpack_get_numeric.c
umfpack_get_symbolic.c
umfpack_load_numeric.c
umfpack_load_symbolic.c
umfpack_numeric.c
umfpack_qsymbolic.c
umfpack_report_control.c
umfpack_report_info.c
umfpack_report_matrix.c
umfpack_report_numeric.c
umfpack_report_perm.c
umfpack_report_status.c
umfpack_report_symbolic.c
umfpack_report_triplet.c
umfpack_report_vector.c
umfpack_save_numeric.c
umfpack_save_symbolic.c
umfpack_scale.c
umfpack_solve.c
umfpack_symbolic.c
umfpack_transpose.c
umfpack_triplet_to_col.c
umf_realloc.c
umf_report_perm.c
umf_report_vector.c
umf_row_search.c
umf_scale_column.c
umf_scale.c
umf_set_stats.c
umf_singletons.c
umf_solve.c
umf_start_front.c
umf_store_lu.c
umf_symbolic_usage.c
umf_transpose.c
umf_triplet.c
umf_tuple_lengths.c
umf_usolve.c
umf_utsolve.c
umf_valid_numeric.c
umf_valid_symbolic.c
) # DINT only


ADD_LIBRARY(GROUP1_OBJS OBJECT ${GROUP1_SRCS})
target_compile_definitions(GROUP1_OBJS PRIVATE DINT)

SET (GROUP2_SRCS
umf_2by2.c
umf_assemble.c
umf_blas3_update.c
umf_build_tuples.c
umf_create_element.c
umf_dump.c
umf_extend_front.c
umf_garbage_collection.c
umf_get_memory.c
umf_grow_front.c
umf_init_front.c
umf_kernel.c
umf_kernel_init.c
umf_kernel_wrapup.c
umf_local_search.c
umf_lsolve.c
umf_ltsolve.c
umf_mem_alloc_element.c
umf_mem_alloc_head_block.c
umf_mem_alloc_tail_block.c
umf_mem_free_tail_block.c
umf_mem_init_memoryspace.c
umfpack_col_to_triplet.c
umfpack_defaults.c
umfpack_free_numeric.c
umfpack_free_symbolic.c
umfpack_get_determinant.c
umfpack_get_lunz.c
umfpack_get_numeric.c
umfpack_get_symbolic.c
umfpack_load_numeric.c
umfpack_load_symbolic.c
umfpack_numeric.c
umfpack_qsymbolic.c
umfpack_report_control.c
umfpack_report_info.c
umfpack_report_matrix.c
umfpack_report_numeric.c
umfpack_report_perm.c
umfpack_report_status.c
umfpack_report_symbolic.c
umfpack_report_triplet.c
umfpack_report_vector.c
umfpack_save_numeric.c
umfpack_save_symbolic.c
umfpack_scale.c
umfpack_solve.c
umfpack_symbolic.c
umfpack_transpose.c
umfpack_triplet_to_col.c
umf_report_vector.c
umf_row_search.c
umf_scale_column.c
umf_scale.c
umf_set_stats.c
umf_solve.c
umf_start_front.c
umf_store_lu.c
umf_symbolic_usage.c
umf_transpose.c
umf_triplet.c
umf_tuple_lengths.c
umf_usolve.c
umf_utsolve.c
umf_valid_numeric.c
umf_valid_symbolic.c
)
ADD_LIBRARY(GROUP2_OBJS OBJECT ${GROUP2_SRCS})
target_compile_definitions(GROUP2_OBJS PRIVATE ZINT)

SET (GROUP3_SRCS
umf_ltsolve.c
umf_utsolve.c
)
ADD_LIBRARY(GROUP3_OBJS OBJECT ${GROUP3_SRCS})
target_compile_definitions(GROUP3_OBJS PRIVATE DINT CONJUGATE_SOLVE)

SET(GROUP4_SRCS
  umf_triplet.c
)
ADD_LIBRARY(GROUP4_OBJS OBJECT ${GROUP4_SRCS})
target_compile_definitions(GROUP4_OBJS PRIVATE DINT DO_MAP)


SET(GROUP5_SRCS
umf_triplet.c
)
ADD_LIBRARY(GROUP5_OBJS OBJECT ${GROUP5_SRCS})
target_compile_definitions(GROUP5_OBJS PRIVATE DINT DO_VALUES)

SET(GROUP5A_SRCS
umf_triplet.c
)
ADD_LIBRARY(GROUP5A_OBJS OBJECT ${GROUP5A_SRCS})
target_compile_definitions(GROUP5A_OBJS PRIVATE DINT DO_MAP DO_VALUES)

SET(GROUP6_SRCS
umf_store_lu.c
)
ADD_LIBRARY(GROUP6_OBJS OBJECT ${GROUP6_SRCS})
target_compile_definitions(GROUP6_OBJS PRIVATE DINT DROP)

SET(GROUP7_SRCS
umf_assemble.c
)
ADD_LIBRARY(GROUP7_OBJS OBJECT ${GROUP7_SRCS})
target_compile_definitions(GROUP7_OBJS PRIVATE DINT FIXQ)

SET(GROUP8_SRCS
umfpack_solve.c
)
ADD_LIBRARY(GROUP8_OBJS OBJECT ${GROUP8_SRCS})
target_compile_definitions(GROUP8_OBJS PRIVATE DINT WSOLVE)

SET(GROUP9_SRCS umf_ltsolve.c)
ADD_LIBRARY(GROUP9_OBJS OBJECT ${GROUP9_SRCS})
target_compile_definitions(GROUP9_OBJS PRIVATE ZINT CONJUGATE_SOLVE)

SET(GROUP10_SRCS umf_utsolve.c)
ADD_LIBRARY(GROUP10_OBJS OBJECT ${GROUP10_SRCS})
target_compile_definitions(GROUP10_OBJS PRIVATE ZINT CONJUGATE_SOLVE)

SET(GROUP11_SRCS umf_triplet.c)
ADD_LIBRARY(GROUP11_OBJS OBJECT ${GROUP11_SRCS})
target_compile_definitions(GROUP11_OBJS PRIVATE ZINT DO_MAP)

SET(GROUP12_SRCS umf_triplet.c)
ADD_LIBRARY(GROUP12_OBJS OBJECT ${GROUP12_SRCS})
target_compile_definitions(GROUP12_OBJS PRIVATE ZINT DO_MAP DO_VALUES)

SET(GROUP13_SRCS umf_triplet.c)
ADD_LIBRARY(GROUP13_OBJS OBJECT ${GROUP13_SRCS})
target_compile_definitions(GROUP13_OBJS PRIVATE ZINT DO_VALUES)

SET(GROUP14_SRCS umf_store_lu.c)
ADD_LIBRARY(GROUP14_OBJS OBJECT ${GROUP14_SRCS})
target_compile_definitions(GROUP14_OBJS PRIVATE ZINT DROP)

SET(GROUP15_SRCS umf_assemble.c)
ADD_LIBRARY(GROUP15_OBJS OBJECT ${GROUP15_SRCS})
target_compile_definitions(GROUP15_OBJS PRIVATE ZINT FIXQ)

SET(GROUP16_SRCS umfpack_solve.c)
ADD_LIBRARY(GROUP16_OBJS OBJECT ${GROUP16_SRCS})
target_compile_definitions(GROUP16_OBJS PRIVATE ZINT WSOLVE)

ADD_LIBRARY(umfpack SHARED
${GLOBAL_SRCS}
$<TARGET_OBJECTS:GROUP1_OBJS>
$<TARGET_OBJECTS:GROUP2_OBJS>
$<TARGET_OBJECTS:GROUP3_OBJS>
$<TARGET_OBJECTS:GROUP4_OBJS>
$<TARGET_OBJECTS:GROUP5_OBJS>
$<TARGET_OBJECTS:GROUP5A_OBJS>
$<TARGET_OBJECTS:GROUP6_OBJS>
$<TARGET_OBJECTS:GROUP7_OBJS>
$<TARGET_OBJECTS:GROUP8_OBJS>
$<TARGET_OBJECTS:GROUP9_OBJS>
$<TARGET_OBJECTS:GROUP10_OBJS>
$<TARGET_OBJECTS:GROUP11_OBJS>
$<TARGET_OBJECTS:GROUP12_OBJS>
$<TARGET_OBJECTS:GROUP13_OBJS>
$<TARGET_OBJECTS:GROUP14_OBJS>
$<TARGET_OBJECTS:GROUP15_OBJS>
$<TARGET_OBJECTS:GROUP16_OBJS>
)
TARGET_LINK_LIBRARIES(umfpack PRIVATE  amd blaswrapper -lm)

