/* ========================================================================== */
/* === umfpack_col_to_triplet =============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

UMFPACK_PUBLIC
int umfpack_di_col_to_triplet
(
    int n_col,
    const int Ap [ ],
    int Tj [ ]
) ;

UF_long umfpack_dl_col_to_triplet
(
    UF_long n_col,
    const UF_long Ap [ ],
    UF_long Tj [ ]
) ;

UMFPACK_PUBLIC
int umfpack_zi_col_to_triplet
(
    int n_col,
    const int Ap [ ],
    int Tj [ ]
) ;

UF_long umfpack_zl_col_to_triplet
(
    UF_long n_col,
    const UF_long Ap [ ],
    UF_long Tj [ ]
) ;

/*
double int Syntax:

    #include "umfpack.h"
    int n_col, *Tj, *Ap, status ;
    status = umfpack_di_col_to_triplet (n_col, Ap, Tj) ;

double UF_long Syntax:

    #include "umfpack.h"
    UF_long n_col, *Tj, *Ap, status ;
    status = umfpack_dl_col_to_triplet (n_col, Ap, Tj) ;

complex int Syntax:

    #include "umfpack.h"
    int n_col, *Tj, *Ap, status ;
    status = umfpack_zi_col_to_triplet (n_col, Ap, Tj) ;

complex UF_long Syntax:

    #include "umfpack.h"
    UF_long n_col, *Tj, *Ap, status ;
    status = umfpack_zl_col_to_triplet (n_col, Ap, Tj) ;

Purpose:

    Converts a column-oriented matrix to a triplet form.  Only the column
    pointers, Ap, are required, and only the column indices of the triplet form
    are constructed.   This routine is the opposite of umfpack_*_triplet_to_col.
    The matrix may be singular and/or rectangular.  Analogous to [i, Tj, x] =
    find (A) in MATLAB, except that zero entries present in the column-form of
    A are present in the output, and i and x are not created (those are just Ai
    and Ax+Az*1i, respectively, for a column-form matrix A).

Returns:

    UMFPACK_OK if successful
    UMFPACK_ERROR_argument_missing if Ap or Tj is missing
    UMFPACK_ERROR_n_nonpositive if n_col <= 0
    UMFPACK_ERROR_invalid_matrix if Ap [n_col] < 0, Ap [0] != 0, or
	Ap [j] > Ap [j+1] for any j in the range 0 to n-1.
    Unsorted columns and duplicate entries do not cause an error (these would
    only be evident by examining Ai).  Empty rows and columns are OK.

Arguments:

    Int n_col ;		Input argument, not modified.

	A is an n_row-by-n_col matrix.  Restriction: n_col > 0.
	(n_row is not required)

    Int Ap [n_col+1] ;	Input argument, not modified.

	The column pointers of the column-oriented form of the matrix.  See
	umfpack_*_*symbolic for a description.  The number of entries in
	the matrix is nz = Ap [n_col].  Restrictions on Ap are the same as those
	for umfpack_*_transpose.  Ap [0] must be zero, nz must be >= 0, and
	Ap [j] <= Ap [j+1] and Ap [j] <= Ap [n_col] must be true for all j in
	the range 0 to n_col-1.  Empty columns are OK (that is, Ap [j] may equal
	Ap [j+1] for any j in the range 0 to n_col-1).

    Int Tj [nz] ;	Output argument.

	Tj is an integer array of size nz on input, where nz = Ap [n_col].
	Suppose the column-form of the matrix is held in Ap, Ai, Ax, and Az
	(see umfpack_*_*symbolic for a description).  Then on output, the
	triplet form of the same matrix is held in Ai (row indices), Tj (column
	indices), and Ax (numerical values).  Note, however, that this routine
	does not require Ai and Ax (or Az for the complex version) in order to
	do the conversion.
*/
