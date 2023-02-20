
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

#
#    /* ---------------------------------------------------------------------- */
#    /* initializations */
#    /* ---------------------------------------------------------------------- */
#
umfclass = umf.umf_control(dll, matrix_type='real')
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
triplet_matrix = umf.di_triplet(umfclass, Arow, Acol, Aval)
umfclass.print_triplet(triplet_matrix)

#    /* convert to column form */
matrix = umfclass.triplet_to_col(triplet_matrix)

#    /* print the column-form of A */
umfclass.print_matrix(matrix)

#    /* ---------------------------------------------------------------------- */
#    /* symbolic factorization */
#    /* ---------------------------------------------------------------------- */
#
Symbolic = umfclass.symbolic(matrix)
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
Numeric = umfclass.numeric(matrix, Symbolic)


#    /* print the numeric factorization */
umfclass.print_numeric(Numeric)

#
#    /* ---------------------------------------------------------------------- */
#    /* solve Ax=b */
#    /* ---------------------------------------------------------------------- */
#
status = umfclass.solve(matrix, x, b, Numeric, False)
umfclass.print_info()
umfclass.print_status()

print("\nx (solution of Ax=b): ")
umfclass.print_vector(x, 'x')
rnorm = umf.resid(transpose=False, matrix=matrix, x=x, r=r, b=b)
print ("maxnorm of residual: %g\n\n" % rnorm)

#   /* ---------------------------------------------------------------------- */
#   /* compute the determinant */
#   /* ---------------------------------------------------------------------- */

umfclass.determinant(x, r, Numeric)
umfclass.print_determinant(x, r)

umfclass.toc()

