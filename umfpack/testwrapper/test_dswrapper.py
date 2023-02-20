
import array
import dswrapper
import umfpack_loader as umf


is_complex = False
if is_complex:
    umfclass = dswrapper.solver(action="init", blas=umf.get_blas_name(), matrix_type='complex')
else:
    umfclass = dswrapper.solver(action="init", blas=umf.get_blas_name(), matrix_type='real')

n = 5
nz = 12
Arow = array.array('i', (0,  4,  1,  1,   2,   2,  0,  1,  2,  3,  4,  4))
Acol = array.array('i', (0,  4,  0,  2,   1,   2,  1,  4,  3,  2,  1,  2))
if is_complex:
    Aval = array.array('d', (2., 0., 1., 0., 3., 0., 4., 0., -1., 0., -3., 0., 3., 0., 6., 0., 2., 0., 1., 0., 4., 0., 2., 0.))
    b = array.array('d', (8., 0., 45., 0., -3., 0., 3., 0., 19., 0.))
    x = array.array('d', [0.0]*2*n)
    r = array.array('d', [0.0]*2*n)
else:
    Aval = array.array('d', (2., 1., 3., 4., -1., -3., 3., 6., 2., 1., 4., 2.))
    b = array.array('d', (8., 45., -3., 3., 19.))
    x = array.array('d', [0.0]*n)
    r = array.array('d', [0.0]*n)

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

matrix = dswrapper.solver(action="matrix", math_object=umfclass, matrix_format="csc", Ap=matrix.Ap, Ai=matrix.Ai, Ax=matrix.Ax)

#    /* print the column-form of A */
umfclass.print_matrix(matrix)

#    /* ---------------------------------------------------------------------- */
#    /* symbolic factorization */
#    /* ---------------------------------------------------------------------- */

Symbolic = dswrapper.solver(action="symbolic", math_object=umfclass, matrix=matrix)

#    /* print the symbolic factorization */

umfclass.print_symbolic(Symbolic)


#    /* ---------------------------------------------------------------------- */
#    /* numeric factorization */
#    /* ---------------------------------------------------------------------- */

Numeric = dswrapper.solver(action="numeric", math_object=umfclass, matrix=matrix, symbolic=Symbolic)

#    /* print the numeric factorization */
umfclass.print_numeric(Numeric)


#    /* ---------------------------------------------------------------------- */
#    /* solve Ax=b */
#    /* ---------------------------------------------------------------------- */

dswrapper.solver(action="solve", math_object=umfclass, matrix=matrix, numeric=Numeric, b=b, x=x, transpose=False)
umfclass.print_info()
umfclass.print_status()

print("\nx (solution of Ax=b): ")
umfclass.print_vector(x, 'x')
rnorm = umf.resid(transpose=False, is_complex=is_complex, matrix=matrix, x=x, r=r, b=b)
print ("maxnorm of residual: %g\n\n" % rnorm)

