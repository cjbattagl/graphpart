#include <matrix.h>
#include <mex.h>   

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

//declare variables
    mxArray *a_in_m, *b_out_m;
    const mwSize *dims;
    double *a, *b;
    int dimx, dimy, numdims;
    int i,j;
    idx_t nparts;
    double gamma;

    if (nrhs != 4 || nlhs > 1) {
        mexErrMsgTxt ("Wrong # of arguments");
    }
    if (!mxIsSparse(prhs[0]) || mxGetN(prhs[0])!=mxGetM(prhs[0])) {
        mexErrMsgTxt ("First parameter must be a symmetric sparse matrix");
    }
    if (!mxIsMatrix(prhs[3])) {
        mexErrMsgTxt ("Fourth parameter must be a vorder");
    }

//associate inputs
    a_in_m = mxDuplicateArray(prhs[0]);
    v_in_m = mxDuplicateArray(prhs[3]);
    nparts = (idx_t) mxGetScalar (prhs[1]);
    gamma = (double) mxGetScalar (prhs[2]);

    if (nparts < 2) {
      mexErrMsgTxt("nparts must be at least 2");
    }
//figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    dimy = (int)dims[0]; dimx = (int)dims[1];

//associate outputs
    b_out_m = plhs[0] = mxCreateDoubleMatrix(dimy,1,mxREAL);

//associate pointers
    a = mxGetPr(a_in_m);
    v = mxGetPr(v_in_m);
    b = mxGetPr(b_out_m);

    mwIndex n, nnz, *jc, *ir;
    double *pr;

    n = mxGetN(a_in_m);
    jc = mxGetJc(a_in_m);
    nnz = jc[n];
    ir = mxGetIr(a_in_m);
    pr = mxGetPr(a_in_m);
    
    

//outer loop: all vertices
    for(i=0;i<dimx;i++) {
      int curr_v = v[i];
          for (i = 1; i <= n; i++) {
        for (j = jc[i-1]; j < jc[i]; j++) {
    }


    for(i=0;i<dimx;i++)
    {
        for(j=0;j<dimy;j++)
        {
            mexPrintf("element[%d][%d] = %f\n",j,i,a[i*dimy+j]);
           // c[i*dimy+j] = a[i*dimy+j]+5; //adds 5 to every element in a
           // d[i*dimy+j] = b[i*dimy+j]*b[i*dimy+j]; //squares b
        }
    }

    return;
}
