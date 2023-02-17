from ctypes import *
import platform

#/* used in all UMFPACK_report_* routines: */
UMFPACK_PRL = 0                  # /* print level */

#/* used in UMFPACK_*symbolic only: */
UMFPACK_DENSE_ROW = 1            #/* dense row parameter */
UMFPACK_DENSE_COL = 2            #/* dense col parameter */
UMFPACK_BLOCK_SIZE = 4           #/* BLAS-3 block size */
UMFPACK_STRATEGY = 5             #/* auto, symmetric, unsym., or 2by2 */
UMFPACK_2BY2_TOLERANCE = 12      #/* 2-by-2 pivot tolerance */
UMFPACK_FIXQ = 13                #/* -1: no fixQ, 0: default, 1: fixQ */
UMFPACK_AMD_DENSE = 14           #/* for AMD ordering */
UMFPACK_AGGRESSIVE = 19          #/* whether or not to use aggressive

#/* default values of Control may change in future versions of UMFPACK. */

#/* -------------------------------------------------------------------------- */
#/* status codes */
#/* -------------------------------------------------------------------------- */
#
UMFPACK_OK = 0
#
#/* status > 0 means a warning, but the method was successful anyway. */
#/* A Symbolic or Numeric object was still created. */
##define UMFPACK_WARNING_singular_matrix (1)
#
#/* The following warnings were added in umfpack_*_get_determinant */
##define UMFPACK_WARNING_determinant_underflow (2)
##define UMFPACK_WARNING_determinant_overflow (3)
#
#/* status < 0 means an error, and the method was not successful. */
#/* No Symbolic of Numeric object was created. */
##define UMFPACK_ERROR_out_of_memory (-1)
##define UMFPACK_ERROR_invalid_Numeric_object (-3)
##define UMFPACK_ERROR_invalid_Symbolic_object (-4)
##define UMFPACK_ERROR_argument_missing (-5)
##define UMFPACK_ERROR_n_nonpositive (-6)
##define UMFPACK_ERROR_invalid_matrix (-8)
##define UMFPACK_ERROR_different_pattern (-11)
##define UMFPACK_ERROR_invalid_system (-13)
##define UMFPACK_ERROR_invalid_permutation (-15)
##define UMFPACK_ERROR_internal_error (-911) /* yes, call me if you get this! */
##define UMFPACK_ERROR_file_IO (-17)

#/* -------------------------------------------------------------------------- */
#/* solve codes */
#/* -------------------------------------------------------------------------- */
#
#/* Solve the system ( )x=b, where ( ) is defined below.  "t" refers to the */
#/* linear algebraic transpose (complex conjugate if A is complex), or the (') */
#/* operator in MATLAB.  "at" refers to the array transpose, or the (.') */
#/* operator in MATLAB. */
#
UMFPACK_A = (0)     # /* Ax=b    */
UMFPACK_At = (1) #    /* A'x=b   */
##define UMFPACK_Aat    (2)     /* A.'x=b  */
#
UMFPACK_Pt_L = (3)  #   /* P'Lx=b  */
##define UMFPACK_L      (4)     /* Lx=b    */
##define UMFPACK_Lt_P   (5)     /* L'Px=b  */
##define UMFPACK_Lat_P  (6)     /* L.'Px=b */
##define UMFPACK_Lt     (7)     /* L'x=b   */
##define UMFPACK_Lat    (8)     /* L.'x=b  */
#
UMFPACK_U_Qt = (9) #    /* UQ'x=b  */
##define UMFPACK_U      (10)    /* Ux=b    */
##define UMFPACK_Q_Ut   (11)    /* QU'x=b  */
##define UMFPACK_Q_Uat  (12)    /* QU.'x=b */
##define UMFPACK_Ut     (13)    /* U'x=b   */
##define UMFPACK_Uat    (14)    /* U.'x=b  */
#
#/* -------------------------------------------------------------------------- */
#
#/* Integer constants are used for status and solve codes instead of enum */
#/* to make it easier for a Fortran code to call UMFPACK. */

# /* -------------------------------------------------------------------------- */
# /* size of Info and Control arrays */
# /* -------------------------------------------------------------------------- */
#
# /* These might be larger in future versions, since there are only 3 unused
#  * entries in Info, and no unused entries in Control. */

UMFPACK_INFO = 90
UMFPACK_CONTROL = 20

def get_dll_naming():
    systems = {
      'Linux' : {"prefix" : "lib", "suffix" : ".so"},
      'Darwin' : {"prefix" : "lib", "suffix" : ".dylib"},
      'Windows' : {"prefix" : "", "suffix" : ".dll"},
    }
    return systems[platform.system()]

def get_umfpack_name():
    return "./%(prefix)sumfpack_lgpl%(suffix)s" % get_dll_naming()

def load_umfpack_dll():
    dll = cdll.LoadLibrary(get_umfpack_name())
    if not dll:
        raise RuntimeError("Cannot find UMFPACK dll")
    return dll

def get_blas_name():
    return "%(prefix)sopenblas%(suffix)s" % get_dll_naming()

def load_blas_dll(dll):
    #print(dll.blasw_load_dll)
    libname = c_char_p(get_blas_name().encode('utf-8'))
    #libname = c_char_p(b'mkl_rt.2.dll')
    msg = c_char_p()
    dll.blasw_load_dll.argtypes = [c_char_p, c_void_p]
    dll.blasw_load_dll.restype = c_void_p
    h = dll.blasw_load_dll(libname, byref(msg))
    if not h or msg:
      if msg:
          print(string_at(msg))
      raise RuntimeError("NO BLAS DLL LOADED")
    #print(h)
    return h

def load_blas_functions(dll, h):
    dll.blasw_load_functions.restype = c_int
    dll.blasw_load_functions.argtypes = [c_void_p]
    i = dll.blasw_load_functions(h)
    return i

def myprintcb(msg):
    pmsg = msg.decode("utf-8")
    pmsg = pmsg.replace("\n", "\nUMF: ")
    print(pmsg, end="")

global_callback = None

def set_python_print_callback(dll):
    global global_callback
    CALLBACK = CFUNCTYPE(None, c_char_p)
    dll.blasw_set_printer_callback.argtypes = [CALLBACK]
    dll.blasw_set_printer_callback.restype = None
    dcb = CALLBACK(myprintcb)
    dll.blasw_set_printer_callback(dcb)
    global_callback = dcb
    return dcb

def get_info():
    return (c_int * uml.UMFPACK_INFO)();

