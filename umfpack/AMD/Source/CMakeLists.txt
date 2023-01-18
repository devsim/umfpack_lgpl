
INCLUDE_DIRECTORIES(../Include ${CMAKE_SOURCE_DIR}/UFconfig)
SET (AMD_SRCS
amd_aat.c
amd_1.c
amd_2.c
amd_dump.c
amd_postorder.c
amd_post_tree.c
amd_defaults.c
amd_order.c
amd_control.c
amd_info.c
amd_valid.c
amd_preprocess.c
amd_global.c
)
ADD_LIBRARY(amd STATIC ${AMD_SRCS})
TARGET_COMPILE_DEFINITIONS(amd PRIVATE DINT)
