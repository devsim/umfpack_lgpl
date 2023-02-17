
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

###    /* ---------------------------------------------------------------------- */
###    /* solve Ax=b, broken down into steps */
###    /* ---------------------------------------------------------------------- */
###
###    /* Rb = R*b */
#Rb = (c_double * (n))()
#y  = (c_double * (n))()
#
#status = dll.umfpack_di_scale (Rb, b, Numeric)
#if (status < 0):
#  raise RuntimeError ("umfpack_di_scale failed")
#
###    /* solve Ly = P*(Rb) */
#status = dll.umfpack_di_solve (umf.UMFPACK_Pt_L, Ap, Ai, Ax, y, Rb, Numeric, Control, Info)
#if (status < 0):
#  raise RuntimeError ("umfpack_di_solve failed")
#
###    /* solve UQ'x=y */
#status = dll.umfpack_di_solve (umf.UMFPACK_U_Qt, Ap, Ai, Ax, x, y, Numeric, Control, Info)
#if (status < 0):
#  raise RuntimeError ("umfpack_di_solve failed")
#
#print ("\nx (solution of Ax=b, solve is split into 3 steps): ")
#dll.umfpack_di_report_vector (n, x, Control)
#rnorm = resid (False, Ap, Ai, Ax)
#print ("maxnorm of residual: %g\n\n" % rnorm)
#
###    /* ---------------------------------------------------------------------- */
###    /* solve A'x=b */
###    /* ---------------------------------------------------------------------- */
###
#status = dll.umfpack_di_solve (umf.UMFPACK_At, Ap, Ai, Ax, x, b, Numeric, Control, Info)
#dll.umfpack_di_report_info (Control, Info)
#if (status < 0):
#  raise RuntimeError ("umfpack_di_solve failed")
#
#print ("\nx (solution of A'x=b): ")
#dll.umfpack_di_report_vector (n, x, Control)
#rnorm = resid (True, Ap, Ai, Ax)
#print ("maxnorm of residual: %g\n\n" % rnorm)
###
###    /* ---------------------------------------------------------------------- */
###    /* modify one numerical value in the column-form of A */
###    /* ---------------------------------------------------------------------- */
###
###    /* change A (1,4), look for row index 1 in column 4. */
#row = 1
#col = 4
#for p in range(Ap [col], Ap[col+1]):
#   if (row == Ai [p]):
#       print ("\nchanging A (%d,%d) to zero\n" % (row, col))
#       Ax [p] = 0.0
#       break
#print ("\nmodified A: ")
#dll.umfpack_di_report_matrix (n, n, Ap, Ai, Ax, 1, Control)
###
###    /* ---------------------------------------------------------------------- */
###    /* redo the numeric factorization */
###    /* ---------------------------------------------------------------------- */
###
###    /* The pattern (Ap and Ai) hasn't changed, so the symbolic factorization */
###    /* doesn't have to be redone, no matter how much we change Ax. */
###
###    /* We don't need the Numeric object any more, so free it. */
#dll.umfpack_di_free_numeric (byref(Numeric))
###
###    /* Note that a memory leak would have occurred if the old Numeric */
###    /* had not been free'd with umfpack_di_free_numeric above. */
#status = dll.umfpack_di_numeric (Ap, Ai, Ax, Symbolic, byref(Numeric), Control, Info)
#if status < 0:
#  dll.umfpack_di_report_info (Control, Info)
#  dll.umfpack_di_report_status (Control, status)
#  raise RuntimeError ("umfpack_di_numeric failed")
#
#print ("\nNumeric factorization of modified A: ")
#dll.umfpack_di_report_numeric (Numeric, Control)
###
###    /* ---------------------------------------------------------------------- */
###    /* solve Ax=b, with the modified A */
###    /* ---------------------------------------------------------------------- */
###
#status = dll.umfpack_di_solve (umf.UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, Control, Info)
#dll.umfpack_di_report_info (Control, Info)
#if status < 0:
#  dll.umfpack_di_report_status (Control, status)
#  raise RuntimeError ("umfpack_di_solve failed")
#
#print ("\nx (with modified A): ")
#dll.umfpack_di_report_vector (n, x, Control)
#rnorm = resid (False, Ap, Ai, Ax)
#print ("maxnorm of residual: %g\n\n" % rnorm)
###
###    /* ---------------------------------------------------------------------- */
###    /* modify all of the numerical values of A, but not the pattern */
###    /* ---------------------------------------------------------------------- */
###
#for col in range(n):
#  for p in range(Ap[col], Ap[col+1]):
#    row = Ai [p]
#    print ("changing A (%d,%d) from %g" % (row, col, Ax [p]), end="")
#    Ax [p] = Ax [p] + col*10 - row
#    print (" to %g\n" % Ax[p])
#print ("\ncompletely modified A (same pattern): ")
#dll.umfpack_di_report_matrix (n, n, Ap, Ai, Ax, 1, Control)
###
###    /* ---------------------------------------------------------------------- */
###    /* save the Symbolic object to file, free it, and load it back in */
###    /* ---------------------------------------------------------------------- */
###
###    /* use the default filename, "symbolic.umf" */
#print ("\nSaving symbolic object:\n")
#status = dll.umfpack_di_save_symbolic (Symbolic, NULL)
#if (status < 0):
#  dll.umfpack_di_report_status (Control, status)
#  raise RuntimeError ("umfpack_di_save_symbolic failed")
#
#print ("\nFreeing symbolic object:\n")
#dll.umfpack_di_free_symbolic (byref(Symbolic))
#print ("\nLoading symbolic object:\n")
#status = dll.umfpack_di_load_symbolic (byref(Symbolic), NULL)
#if (status < 0):
#  dll.umfpack_di_report_status (Control, status)
#  raise RuntimeError ("umfpack_di_load_symbolic failed")
#print ("\nDone loading symbolic object\n")
#
###
###    /* ---------------------------------------------------------------------- */
###    /* redo the numeric factorization */
###    /* ---------------------------------------------------------------------- */
###
#dll.umfpack_di_free_numeric (byref(Numeric))
#status = dll.umfpack_di_numeric (Ap, Ai, Ax, Symbolic, byref(Numeric), Control, Info)
#if (status < 0):
#  dll.umfpack_di_report_info (Control, Info)
#  dll.umfpack_di_report_status (Control, status)
#  raise RuntimeError ("umfpack_di_numeric failed")
#print ("\nNumeric factorization of completely modified A: ")
#dll.umfpack_di_report_numeric (Numeric, Control)
###
###    /* ---------------------------------------------------------------------- */
###    /* solve Ax=b, with the modified A */
###    /* ---------------------------------------------------------------------- */
###
#status = dll.umfpack_di_solve (umf.UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, Control, Info)
#dll.umfpack_di_report_info (Control, Info)
#if (status < 0):
#  dll.umfpack_di_report_status (Control, status)
#  raise RuntimeError ("umfpack_di_solve failed")
#print ("\nx (with completely modified A): ")
#dll.umfpack_di_report_vector (n, x, Control)
#rnorm = resid (False, Ap, Ai, Ax)
#print ("maxnorm of residual: %g\n\n" % rnorm)
#
###    /* ---------------------------------------------------------------------- */
###    /* free the symbolic and numeric factorization */
###    /* ---------------------------------------------------------------------- */
###
#dll.umfpack_di_free_symbolic (byref(Symbolic))
#dll.umfpack_di_free_numeric (byref(Numeric))
###
###    /* ---------------------------------------------------------------------- */
###    /* C = transpose of A */
###    /* ---------------------------------------------------------------------- */
###
#Cp = (c_int * (n + 1))()
#Ci = (c_int * (nz1 + 1))()
#Cx = (c_double * (nz1 + 1))()
#
#status = dll.umfpack_di_transpose (n, n, Ap, Ai, Ax, NULL, NULL, Cp, Ci, Cx)
#if (status < 0):
#  dll.umfpack_di_report_status (Control, status)
#  raise RuntimeError ("umfpack_di_transpose failed: ")
#
#print ("\nC (transpose of A): ")
#dll.umfpack_di_report_matrix (n, n, Cp, Ci, Cx, 1, Control)
#
###    /* ---------------------------------------------------------------------- */
###    /* symbolic factorization of C */
###    /* ---------------------------------------------------------------------- */
###
#status = dll.umfpack_di_symbolic (n, n, Cp, Ci, Cx, byref(Symbolic), Control, Info)
#if (status < 0):
#  dll.umfpack_di_report_info (Control, Info)
#  dll.umfpack_di_report_status (Control, status)
#  raise RuntimeError ("umfpack_di_symbolic failed")
#print ("\nSymbolic factorization of C: ")
#dll.umfpack_di_report_symbolic (Symbolic, Control)
#
###    /* ---------------------------------------------------------------------- */
###    /* copy the contents of Symbolic into user arrays print them */
###    /* ---------------------------------------------------------------------- */
###
#print ("\nGet the contents of the Symbolic object for C:\n")
#print ("(compare with dll.umfpack_di_report_symbolic output, above)\n")
#Pinit = (c_int * (n+1))()
#Qinit = (c_int * (n+1))()
#Front_npivcol = (c_int * (n+1))()
#Front_1strow = (c_int * (n+1))()
#Front_leftmostdesc = (c_int * (n+1))()
#Front_parent = (c_int * (n+1))()
#Chain_start = (c_int * (n+1))()
#Chain_maxrows = (c_int * (n+1))()
#Chain_maxcols = (c_int * (n+1))()
#nchains = c_int()
#
###
#status = dll.umfpack_di_get_symbolic (byref(nr), byref(nc), byref(n1), byref(anz), byref(nfr), byref(nchains), Pinit, Qinit, Front_npivcol, Front_parent, Front_1strow, Front_leftmostdesc, Chain_start, Chain_maxrows, Chain_maxcols, Symbolic)
###
#if (status < 0):
#  raise RuntimeError ("symbolic factorization invalid")
###
#print ("From the Symbolic object, C is of dimension %d-by-%d\n" % (nr.value, nc.value))
#print ("   with nz = %d, number of fronts = %d,\n" % (nz, nfr.value))
#print ("   number of frontal matrix chains = %d\n" % nchains.value)
###
#print ("\nPivot columns in each front, and parent of each front:\n")
#k = 0
#for i in range(nfr.value):
#  fnpiv = Front_npivcol [i]
#  print ("    Front %d: parent front: %d number of pivot cols: %d\n" %
#           (i, Front_parent [i], fnpiv))
#  for j in range(fnpiv):
#    col = Qinit [k]
#    print (
#    "        %d-th pivot column is column %d in original matrix\n" %
#        (k, col))
#    k+=1
#
#print ("\nNote that the column ordering, above, will be refined\n")
#print ("in the numeric factorization below.  The assignment of pivot\n")
#print ("columns to frontal matrices will always remain unchanged.\n")
#
#print ("\nTotal number of pivot columns in frontal matrices: %d\n" % k)
#
#print ("\nFrontal matrix chains:\n")
#for j in range(nchains.value):
#  print ("   Frontal matrices %d to %d are factorized in a single\n" %
#     (Chain_start [j], Chain_start [j+1] - 1))
#  print ("        working array of size %d-by-%d\n" %
#     (Chain_maxrows [j], Chain_maxcols [j]))
#
###
###    /* ---------------------------------------------------------------------- */
###    /* numeric factorization of C */
###    /* ---------------------------------------------------------------------- */
###
#status = dll.umfpack_di_numeric (Cp, Ci, Cx, Symbolic, byref(Numeric), Control, Info)
#if (status < 0):
#  raise RuntimeError ("umfpack_di_numeric failed")
#print ("\nNumeric factorization of C: ")
#dll.umfpack_di_report_numeric (Numeric, Control)
###
###    /* ---------------------------------------------------------------------- */
###    /* extract the LU factors of C and print them */
###    /* ---------------------------------------------------------------------- */
###
#lnz = c_int()
#unz = c_int()
#nzud = c_int()
#if (dll.umfpack_di_get_lunz (byref(lnz), byref(unz), byref(nr), byref(nc), byref(nzud), Numeric) < 0):
#  raise RuntimeError ("umfpack_di_get_lunz failed")
###    /* ensure arrays are not of zero size */
#lnz1 = max(lnz.value,1)
#unz1 = max(unz.value,1)
#Lp = (c_int * (n + 1))()
#Lj = (c_int * lnz1)()
#Lx = (c_double * lnz1)()
#Up = (c_int * (n + 1))()
#Ui = (c_int * unz1)()
#Ux = (c_double * unz1)()
#P = (c_int * n)()
#Q = (c_int * n)()
#Dx = NULL
#Rs  = (c_double * n)()
#do_recip = c_int()
#status = dll.umfpack_di_get_numeric (Lp, Lj, Lx, Up, Ui, Ux, P, Q, Dx, byref(do_recip), Rs, Numeric)
#if (status < 0):
#  raise RuntimeError ("umfpack_di_get_numeric failed")
#
#print ("\nL (lower triangular factor of C): ")
#dll.umfpack_di_report_matrix (n, n, Lp, Lj, Lx, 0, Control)
#print ("\nU (upper triangular factor of C): ")
#dll.umfpack_di_report_matrix (n, n, Up, Ui, Ux, 1, Control)
#print ("\nP: ")
#dll.umfpack_di_report_perm (n, P, Control)
#print ("\nQ: ")
#dll.umfpack_di_report_perm (n, Q, Control)
#print ("\nScale factors: row i of A is to be ")
#if (do_recip):
#  print ("multiplied by the ith scale factor\n")
#else:
#  print ("divided by the ith scale factor\n")
#
#for i in range(n):
#  print ("%d: %g\n" % (i, Rs [i]))
###
###    /* ---------------------------------------------------------------------- */
###    /* convert L to triplet form and print it */
###    /* ---------------------------------------------------------------------- */
###
###    /* Note that L is in row-form, so it is the row indices that are created */
###    /* by dll.umfpack_di_col_to_triplet. */
###
#print ("\nConverting L to triplet form, and printing it:\n")
#Li = (c_int * lnz1)()
#
#if (dll.umfpack_di_col_to_triplet (n, Lp, Li) < 0):
#  raise RuntimeError ("umfpack_di_col_to_triplet failed")
#print ("\nL, in triplet form: ")
#dll.umfpack_di_report_triplet (n, n, lnz, Li, Lj, Lx, Control)
###
###    /* ---------------------------------------------------------------------- */
###    /* save the Numeric object to file, free it, and load it back in */
###    /* ---------------------------------------------------------------------- */
###
###    /* use the default filename, "numeric.umf" */
#print ("\nSaving numeric object:\n")
#status = dll.umfpack_di_save_numeric (Numeric, NULL)
#if (status < 0):
#  dll.umfpack_di_report_status (Control, status)
#  raise RuntimeError ("umfpack_di_save_numeric failed")
#print ("\nFreeing numeric object:\n")
#dll.umfpack_di_free_numeric (byref(Numeric))
#print ("\nLoading numeric object:\n")
#status = dll.umfpack_di_load_numeric (byref(Numeric), NULL)
#if (status < 0):
#  dll.umfpack_di_report_status (Control, status)
#  raise RuntimeError ("umfpack_di_load_numeric failed")
#print ("\nDone loading numeric object\n")
###
###    /* ---------------------------------------------------------------------- */
###    /* solve C'x=b */
###    /* ---------------------------------------------------------------------- */
###
#status = dll.umfpack_di_solve (umf.UMFPACK_At, Cp, Ci, Cx, x, b, Numeric, Control, Info)
#dll.umfpack_di_report_info (Control, Info)
#if (status < 0):
#  dll.umfpack_di_report_status (Control, status)
#  raise RuntimeError ("umfpack_di_solve failed")
#
#print ("\nx (solution of C'x=b): ")
#dll.umfpack_di_report_vector (n, x, Control)
#rnorm = resid (True, Cp, Ci, Cx)
#print ("maxnorm of residual: %g\n\n" % rnorm)
###
###    /* ---------------------------------------------------------------------- */
###    /* solve C'x=b again, using dll.umfpack_di_wsolve instead */
###    /* ---------------------------------------------------------------------- */
###
#print ("\nSolving C'x=b again, using dll.umfpack_di_wsolve instead:\n")
#Wi = (c_int * n)()
#W = (c_double * (5*n))()
#
#status = dll.umfpack_di_wsolve (umf.UMFPACK_At, Cp, Ci, Cx, x, b, Numeric, Control, Info, Wi, W)
#dll.umfpack_di_report_info (Control, Info)
#if (status < 0):
#  dll.umfpack_di_report_status (Control, status)
#  raise RuntimeError ("umfpack_di_wsolve failed")
#print ("\nx (solution of C'x=b): ")
#dll.umfpack_di_report_vector (n, x, Control)
#rnorm = resid (True, Cp, Ci, Cx)
#print ("maxnorm of residual: %g\n\n" % rnorm)
###
#dll.umfpack_di_free_symbolic (byref(Symbolic))
#dll.umfpack_di_free_numeric (byref(Numeric))
###
###    /* ---------------------------------------------------------------------- */
###    /* print the total time spent in this demo */
###    /* ---------------------------------------------------------------------- */
###
#umfclass.toc()
#
