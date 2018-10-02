/*
 * Copyright (c) 2016, Anmol Popli <anmol.ap020@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * @file young.cpp
 * @brief Fast IIR approximation of gaussian filter using Young's approach
 *
 * @author PRAVIN NAIR <anmol.ap020@gmail.com>
 **/

#include <iostream>
using namespace std;
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>
#include <complex.h>
#include <sys/time.h>
#include <math.h>
static double bf[3], bb[3] , B;
static int w;
/**
 * \brief Dynamically allocate 2D array of doubles
 * \param rows      Number of rows
 * \param columns   Number of columns
 * \return pointer to 2D array
 *
 * This routine allocates memory in heap for a 2D
 * array of dimensions rows x columns and datatype
 * double.
 */
double **alloc_array(int rows, int columns)
{
    int i;
    int j;
    /* Allocate an array of pointers with size equal to number of rows */
    double** twoDary = new double*[rows] ;
    double* currentrow;

    /* For each row, allocate an array with size equal to number of columns */
    for ( i = 0; i < rows; i++ )
        twoDary[i] = new double[columns];

    /* Initialize the 2D array with zeros */
    for (j = 0; j < rows; j++) {
        currentrow = *(twoDary + j);
        for ( i = 0; i < columns; i++ ) {
            *(currentrow + i) = 0;
        }
    }
    return twoDary;
}


/**
 * \brief Deallocate dynamically allocated 2D array of doubles
 * \param arr       Pointer to 2D array
 * \param m         Number of rows
 *
 * This routine deallocates heap memory allocated for
 * 2D array of rows m and datatype double.
 */
void dealloc_array(double **arr,int m)
{
    int k;
    /* Free memory corresponding to each row */
    for(k=0;k<m;k++)
    {
        delete [] arr[k];
    }
    /* Free memory corresponding to the array of pointers to rows */
    delete [] arr;
}

void calculate_parameters(int sigma){
    double q;
    if (sigma < 2.5)
        q = 3.97156 - 4.14554*sqrt(1-0.26891*sigma);
    else
        q = 0.98711*sigma - 0.9633;

    /** \brief Filter parameters b0, b1, b2, b3 */
    double b0 = 1.57825 + 2.44413*q + 1.4281*q*q + 0.422205*q*q*q;
    double b1 = 2.44413*q + 2.85619*q*q + 1.26661*q*q*q;
    double b2 = -(1.4281*q*q + 1.26661*q*q*q);
    double b3 = 0.422205*q*q*q;

    /** \brief Filter parameters bf, bb, B */
    bf[0] = b3/b0; bf[1] = b2/b0; bf[2] = b1/b0;
    bb[0] = b1/b0; bb[1] = b2/b0; bb[2] = b3/b0;
    B = 1 - (b1+b2+b3)/b0;
    w = 3*sigma;
}

/**
 * \brief Convolve input array with 1D Gaussian filter
 *        (Young and van Vliet's algorithm) 
 * \param in        Pointer to input array
 * \param datasize  Input array size
 *
 * This routine performs constant time convolution of the
 * 1D input array of complex doubles with 1D Gaussian filter
 * using Young and van Vliet's algorithm. The input array is
 * first convolved with 1D Causal filter, the result of
 * which is convolved with 1D AntiCausal filter.
 */
void convolve_young1D(double* in, int datasize) {
    int i, j;
    in[0] = B*in[0];
    in[1] = B*in[1] + bf[2]*in[0];
    in[2] = B*in[2] + (bf[1]*in[0]+bf[2]*in[1]);
    for (i=3; i<datasize; ++i) 
	in[i] = B*in[i] + bf[0]*in[i-3] + bf[1]*in[i-2] + bf[2]*in[i-1];

    in[datasize-1] = B*in[datasize-1];
    in[datasize-2] = B*in[datasize-2] + bb[0]*in[datasize-1];
    in[datasize-3] = B*in[datasize-3] + (bb[0]*in[datasize-2]+bb[1]*in[datasize-1]);
    for (i=datasize-4; i>=w; --i)
	in[i] = B*in[i] + bb[0]*in[i+1] + bb[1]*in[i+2] + bb[2]*in[i+3];    
}

/**
 * \brief Apply 2D Gaussian filter to input image
 *        (Young and van Vliet's algorithm) 
 * \param ip_padded Pointer to input Matrix
 *
 * This routine applies 2D Gaussian filter of s.d.
 * sigma to input image ip_padded of dimensions
 * rows x columns and computes output image op_padded.
 * 1D filter is first convolved along rows and then
 * along columns. The 1D convolution is performed using
 * Young and van Vliet's fast recursive algorithm.
 */
double **convolve_young2D(double** ip_in,int rows, int columns, double sigma){
    double q;
    if (sigma < 2.5)
        q = 3.97156 - 4.14554*sqrt(1-0.26891*sigma);
    else
        q = 0.98711*sigma - 0.9633;

    /** \brief Filter parameters b0, b1, b2, b3 */
    double b0 = 1.57825 + 2.44413*q + 1.4281*q*q + 0.422205*q*q*q;
    double b1 = 2.44413*q + 2.85619*q*q + 1.26661*q*q*q;
    double b2 = -(1.4281*q*q + 1.26661*q*q*q);
    double b3 = 0.422205*q*q*q;
    /** \brief Filter parameters bf, bb, B */
    bf[0] = b3/b0; bf[1] = b2/b0; bf[2] = b1/b0;
    bb[0] = b1/b0; bb[1] = b2/b0; bb[2] = b3/b0;
    B = 1 - (b1+b2+b3)/b0;
    w = 3*sigma;
    double **ip_padded= alloc_array(rows+2*w,columns+2*w);
    int i,j;
    for (i=0; i<rows; ++i) {
        for (j=0; j<columns; ++j) {
            ip_padded[i+w][j+w]=ip_in[i][j];
        }
    }
    /* Convolve each row with 1D Gaussian filter */
    for (i=0; i<rows+2*w; ++i) 
        convolve_young1D(ip_padded[i], columns+2*w);
    double* intemp=new double[rows+(2*w)];
    double **op_in=alloc_array(rows,columns);
    for (j=w; j<columns+w; j++) 
    {
        /* Convolve each column with 1D Gaussian filter */
        for (i=0;i<rows+(2*w);i++)
            intemp[i]=ip_padded[i][j];      
        convolve_young1D(intemp, rows+2*w);
        /* Store the convolved column in row of output matrix*/
        for (i=w;i<rows+w;i++)
            op_in[i-w][j-w]=intemp[i];
    }
    dealloc_array(ip_padded,rows+2*w);
    delete [] intemp;
    return op_in;
}

#ifdef MATLAB_MEX_FILE  /* Only used if compiling as a MATLAB MEX function */
#include "mex.h"
#define IMAGE_IN        prhs[0]
#define SIGMA_IN        prhs[1]
#define IMAGE_OUT       plhs[0]

void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
    double sigma;
    const mwSize *size;
    long k, K, numpixels;
    int numsteps;

    if(nrhs < 2)
        mexErrMsgTxt("Two input arguments required.");
    else if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");
    if(mexCallMATLAB(1, &IMAGE_OUT, 1, (mxArray **)&IMAGE_IN, "single"))
        mexErrMsgTxt("First argument must be a numeric array.");
    if(!mxIsNumeric(SIGMA_IN) || mxGetNumberOfElements(SIGMA_IN) != 1
        || (sigma = (double)mxGetScalar(SIGMA_IN)) <= 0)
        mexErrMsgTxt("Second argument must be a positive scalar.");

    size = mxGetDimensions(IMAGE_IN);
    int rows=(int)size[0];
    int columns=(int)size[1];
    double *indata=new double[rows*columns];
    indata =mxGetPr(IMAGE_IN);
    IMAGE_OUT = mxCreateDoubleMatrix(rows, columns, mxREAL);
    double *outdata=new double[rows*columns];
    outdata =mxGetPr(IMAGE_OUT);
    double **inmat=alloc_array(rows,columns);
    for(int i=0;i<rows;i++){
	for(int j=0;j<columns;j++)
		inmat[i][j]=indata[j*rows+i];
    }
    double **outmat=convolve_young2D(inmat,rows,columns,sigma);
    dealloc_array(inmat,rows);
    for(int i=0;i<rows;i++){
        for(int j=0;j<columns;j++)
		outdata[j*rows+i]=outmat[i][j];
    }
    dealloc_array(outmat,rows);
    return;
}

#else
int main()
{
    return 0;
}
#endif
