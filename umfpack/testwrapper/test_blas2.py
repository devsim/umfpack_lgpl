
import umfpack_loader as umf
from ctypes import *


dll = umf.load_umfpack_dll()
h = umf.load_blas_dll(dll)
i = umf.load_blas_functions(dll, h)
#print(i)


def myprintcb(msg):
  pmsg = msg.decode("utf-8")
  pmsg = pmsg.replace("\n", "\nUMF: ")
  print(pmsg, end="")

# dcb handle needed
umf.set_python_print_callback(dll)




n = 5
nz = 12
Arow = (c_int * nz)(0,  4,  1,  1,   2,   2,  0,  1,  2,  3,  4,  4)
Acol = (c_int * nz)(0,  4,  0,  2,   1,   2,  1,  4,  3,  2,  1,  2)
Aval = (c_double * nz)(2., 1., 3., 4., -1., -3., 3., 6., 2., 1., 4., 2.)
b = (c_double * n)(8., 45., -3., 3., 19.)
x = (c_double * n)()
r = (c_double * n)()

#/* -------------------------------------------------------------------------- */
#/* resid: compute the residual, r = Ax-b or r = A'x=b and return maxnorm (r) */
#/* -------------------------------------------------------------------------- */
def resid ( transpose, Ap, Ai, Ax):
  n = len(b)
  for i in range(n): 
    r[i] = -b[i]
  if (transpose):
    for j in range(n): 
      for p in range(Ap[j], Ap[j+1]):
        i = Ai [p]
        r [j] += Ax [p] * x [i]
  else:
    for j in range(n):
      for p in range(Ap[j], Ap[j+1]):
        i = Ai [p]
        r [i] += Ax [p] * x [j]
  norm = 0.0
  print(len(r))
  for i in range(n):
    if abs(r[i]) > norm:
      norm = abs(r[i])
  return norm

Info = (c_int * umf.UMFPACK_INFO)()
nr = c_int()
nc = c_int()
n1 = c_int()
anz = c_int()
nfr = c_int()
Symbolic = c_void_p()
Numeric = c_void_p()

#
#    /* ---------------------------------------------------------------------- */
#    /* initializations */
#    /* ---------------------------------------------------------------------- */
#
umfclass = umf.umf_control(dll)

umfclass.tic()

umfclass.set_defaults()
umfclass.init_verbose()

#
#    /* ---------------------------------------------------------------------- */
#    /* print A and b, and convert A to column-form */
#    /* ---------------------------------------------------------------------- */
#
umfclass.print_vector(b, 'b')
#
#    /* print the triplet form of the matrix */

print ("\nA: ")
umfclass.print_triplet(Arow, Acol, Aval)

#    /* convert to column form */
nz1 = max(nz,1) #; /* ensure arrays are not of size zero. */
Ap = (c_int * (n+1))()
Ai = (c_int * (nz1))()
Ax = (c_double * (nz1))()
umfclass.triplet_to_col(Arow, Acol, Aval, Ap, Ai, Ax)

#    /* print the column-form of A */
umfclass.print_matrix(Ap, Ai, Ax)

#    /* ---------------------------------------------------------------------- */
#    /* symbolic factorization */
#    /* ---------------------------------------------------------------------- */
#
umfclass.symbolic(Ap, Ai, Ax, Symbolic, Info)
#
#
#    /* print the symbolic factorization */
#
umfclass.print_symbolic(Symbolic)

#
#    /* ---------------------------------------------------------------------- */
#    /* numeric factorization */
#    /* ---------------------------------------------------------------------- */
#
umfclass.numeric(Ap, Ai, Ax, Symbolic, Numeric, Info)


#    /* print the numeric factorization */
umfclass.print_numeric(Numeric, Info)

#
#    /* ---------------------------------------------------------------------- */
#    /* solve Ax=b */
#    /* ---------------------------------------------------------------------- */
#
status = umfclass.solve(Ap, Ai, Ax, x, b, Numeric, Info)
umfclass.print_info(Info)
umfclass.print_status(status)

print("\nx (solution of Ax=b): ")
umfclass.print_vector(x, 'x')
rnorm = resid(False, Ap, Ai, Ax)
print ("maxnorm of residual: %g\n\n" % rnorm)

#   /* ---------------------------------------------------------------------- */
#   /* compute the determinant */
#   /* ---------------------------------------------------------------------- */

umfclass.determinant(x, r, Numeric, Info)
umfclass.print_determinant(x, r)

#dll.umfpack_di_free_symbolic (byref(Symbolic))
#dll.umfpack_di_free_numeric (byref(Numeric))
umfclass.toc()
#
