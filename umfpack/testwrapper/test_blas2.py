
import umfpack_loader as umf
import array


dll = umf.load_umfpack_dll()
h = umf.load_blas_dll(dll)
i = umf.load_blas_functions(dll, h)
#print(i)

# dcb handle needed
umf.set_python_print_callback(dll)

n = 5
nz = 12
Arow = array.array('i', (0,  4,  1,  1,   2,   2,  0,  1,  2,  3,  4,  4))
Acol = array.array('i', (0,  4,  0,  2,   1,   2,  1,  4,  3,  2,  1,  2))
Aval = array.array('d', (2., 1., 3., 4., -1., -3., 3., 6., 2., 1., 4., 2.))
b = array.array('d', (8., 45., -3., 3., 19.))
x = array.array('d', [0.0]*n)
r = array.array('d', [0.0]*n)

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
[Ap, Ai, Ax] = umfclass.triplet_to_col(Arow, Acol, Aval)

#    /* print the column-form of A */
umfclass.print_matrix(Ap, Ai, Ax)

#    /* ---------------------------------------------------------------------- */
#    /* symbolic factorization */
#    /* ---------------------------------------------------------------------- */
#
Symbolic = umfclass.symbolic(Ap, Ai, Ax)
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
Numeric = umfclass.numeric(Ap, Ai, Ax, Symbolic)


#    /* print the numeric factorization */
umfclass.print_numeric(Numeric)

#
#    /* ---------------------------------------------------------------------- */
#    /* solve Ax=b */
#    /* ---------------------------------------------------------------------- */
#
status = umfclass.solve(Ap, Ai, Ax, x, b, Numeric)
umfclass.print_info()
umfclass.print_status()

print("\nx (solution of Ax=b): ")
umfclass.print_vector(x, 'x')
rnorm = resid(False, Ap, Ai, Ax)
print ("maxnorm of residual: %g\n\n" % rnorm)

#   /* ---------------------------------------------------------------------- */
#   /* compute the determinant */
#   /* ---------------------------------------------------------------------- */

umfclass.determinant(x, r, Numeric)
umfclass.print_determinant(x, r)

#dll.umfpack_di_free_symbolic (byref(Symbolic))
#dll.umfpack_di_free_numeric (byref(Numeric))
umfclass.toc()
#
