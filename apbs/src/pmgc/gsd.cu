/**
 *  @ingroup PMGC
 *  @author  Tucker Beck [fortran ->c translation], Michael Holst [original]
 *  @brief
 *  @version $Id:
 *
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (nathan.baker@pnl.gov)
 * Pacific Northwest National Laboratory
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2010-2014 Battelle Memorial Institute. Developed at the Pacific Northwest National Laboratory, operated by Battelle Memorial Institute, Pacific Northwest Division for the U.S. Department Energy.  Portions Copyright (c) 2002-2010, Washington University in St. Louis.  Portions Copyright (c) 2002-2010, Nathan A. Baker.  Portions Copyright (c) 1999-2002, The Regents of the University of California. Portions Copyright (c) 1995, Michael Holst.
 * All rights reserved.
 *
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * -  Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * - Neither the name of Washington University in St. Louis nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */

#include "gsd.hu"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#define HANDLE_ERROR(x){									\
	cudaError_t _err = x;									\
	if(_err != cudaSuccess){								\
		printf("  (%s:%d)Cuda error: %s\n", __FILE__, __LINE__, cudaGetErrorString(_err));	\
		exit(-1);											\
	}														\
}

__global__ void cuTest(double *x, double *x2, double *fc, double *cc, double *oC, double *uC, double *oE, double *oN, int N, int dx, int dy, int dz){

	int ind = blockDim.x*blockIdx.x + threadIdx.x;
	
	int lb = (dx*dy) + dx + 1;
	int ub = (dz-2)*(dx*dy)+(dy-2)*dx+(dx-2);
	
	if(ind >= lb && ind <= ub && oC[ind] != 0){
		x2[ind] = (2.0/3.0)*( (fc[ind] 
			+ oN[ind]	* x[ind+dx] 
			+ oN[ind-dx]	* x[ind-dx] 
			+ oE[ind]	* x[ind+1] 
			+ oE[ind-1]	* x[ind-1]
			+ uC[ind]	* x[ind+dx*dy] 
			+ uC[ind-dx*dy]	* x[ind-dx*dy] ) 
			/ oC[ind] ) 
			+ cc[ind] + (1.0/3.0)*x[ind];
	}
}

VPUBLIC void Vgsrb(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *ac, double *cc, double *fc,
        double *x, double *w1, double *w2, double *r,
        int *itmax, int *iters,
        double *errtol, double *omega,
        int *iresid, int *iadjoint, int *gpu) {

    int numdia; /// @todo: doc

    MAT2(ac, *nx * *ny * *nz, 1);

    // Do in one step ***
    numdia = VAT(ipc, 11);
    if (numdia == 7) {
       if(*gpu == 0){
		Vgsrb7x(nx, ny, nz,
				ipc, rpc,
				RAT2(ac, 1,1), cc, fc,
				RAT2(ac, 1,2), RAT2(ac, 1,3), RAT2(ac, 1,4),
				x, w1, w2, r,
				itmax, iters, errtol, omega, iresid, iadjoint);
       }
       else
       {
	 clock_t start, stop;
	 start = clock();
    		Vgsrb7xGpu(nx, ny, nz,
				ipc, rpc,
				RAT2(ac, 1,1), cc, fc,
				RAT2(ac, 1,2), RAT2(ac, 1,3), RAT2(ac, 1,4),
				x, w1, w2, r,
				itmax, iters, errtol, omega, iresid, iadjoint);
		stop = clock();
		printf("**One call vo Vgsrb7zGpu: %f\n", ((stop-start)/(double)CLOCKS_PER_SEC)*1000);
       }
    } else if (numdia == 27) {
        Vgsrb27x(nx, ny, nz,
                 ipc, rpc,
                 RAT2(ac, 1, 1), cc, fc,
                 RAT2(ac, 1, 2), RAT2(ac, 1, 3), RAT2(ac, 1, 4),
                 RAT2(ac, 1, 5), RAT2(ac, 1, 6),
                 RAT2(ac, 1, 7), RAT2(ac, 1, 8), RAT2(ac, 1, 9), RAT2(ac, 1,10),
                 RAT2(ac, 1,11), RAT2(ac, 1,12), RAT2(ac, 1,13), RAT2(ac, 1,14),
                 x, w1, w2, r,
                 itmax, iters, errtol, omega, iresid, iadjoint);
    } else {
        Vnm_print(2, "GSRB: invalid stencil type given...\n");
    }
}

VPUBLIC void Vgsrb7xGpu(int *nx,int *ny,int *nz,
        int *ipc, double *rpc,
        double *oC, double *cc, double *fc,
        double *oE, double *oN, double *uC,
        double *x, double *w1, double *w2, double *r,
        int *itmax, int *iters,
        double *errtol, double *omega,
        int *iresid, int *iadjoint){
	
	int i, j, k;
	int sz = *nx * *ny * *nz; 					//<--- grid dimensions
	int threads = 256;							//<--- number of cuda threads per block (max 512) 
	int blocks = (int)ceil(sz/(float)threads); 	//<--- number of block of size threads needed to cover sz grid points

	//this macro creates variable dx_<arr>, dy_<arr>, and dz_<arr> with values nx, ny, and nz respectively
	MAT3(cc, *nx, *ny, *nz);
	MAT3(fc, *nx, *ny, *nz);
	MAT3( x, *nx, *ny, *nz);
	MAT3(w1, *nx, *ny, *nz);
	MAT3(w2, *nx, *ny, *nz);
	MAT3( r, *nx, *ny, *nz);
	MAT3(oE, *nx, *ny, *nz);
	MAT3(oN, *nx, *ny, *nz);
	MAT3(uC, *nx, *ny, *nz);
	MAT3(oC, *nx, *ny, *nz);

	//intialize cuda arrays and allocate the device memory
	double *d_x;  HANDLE_ERROR(cudaMalloc((void**)&d_x,  sizeof(double) * sz));
	double *d_x2; HANDLE_ERROR(cudaMalloc((void**)&d_x2, sizeof(double) * sz));
	double *d_cc; HANDLE_ERROR(cudaMalloc((void**)&d_cc, sizeof(double) * sz));
	double *d_fc; HANDLE_ERROR(cudaMalloc((void**)&d_fc, sizeof(double) * sz));
	double *d_oC; HANDLE_ERROR(cudaMalloc((void**)&d_oC, sizeof(double) * sz));
	double *d_uC; HANDLE_ERROR(cudaMalloc((void**)&d_uC, sizeof(double) * sz));
	double *d_oN; HANDLE_ERROR(cudaMalloc((void**)&d_oN, sizeof(double) * sz));
	double *d_oE; HANDLE_ERROR(cudaMalloc((void**)&d_oE, sizeof(double) * sz));

	//copy data from host to device
	HANDLE_ERROR(cudaMemcpy(d_x,  x,  sizeof(double)*sz, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_x2, x,  sizeof(double)*sz, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_cc, cc, sizeof(double)*sz, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_fc, fc, sizeof(double)*sz, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_oC, oC, sizeof(double)*sz, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_uC, uC, sizeof(double)*sz, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_oN, oN, sizeof(double)*sz, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_oE, oE, sizeof(double)*sz, cudaMemcpyHostToDevice));
	    
	for (*iters=1; *iters<=*itmax; (*iters)++) {

		double *temp;
		
		cuTest<<<blocks, threads>>>(d_x, d_x2, d_fc, d_cc, d_oC, d_uC, d_oE, d_oN, sz, *nx, *ny, *nz);
		HANDLE_ERROR(cudaGetLastError());
		temp = d_x;
		d_x = d_x2;
		d_x2 = temp;
		
		cuTest<<<blocks, threads>>>(d_x, d_x2, d_fc, d_cc, d_oC, d_uC, d_oE, d_oN, sz, *nx, *ny, *nz);
		HANDLE_ERROR(cudaGetLastError());
		temp = d_x;
		d_x = d_x2;
		d_x2 = temp;
		
		cuTest<<<blocks, threads>>>(d_x, d_x2, d_fc, d_cc, d_oC, d_uC, d_oE, d_oN, sz, *nx, *ny, *nz);
		HANDLE_ERROR(cudaGetLastError());
		temp = d_x;
		d_x = d_x2;
		d_x2 = temp;

	}
	HANDLE_ERROR(cudaThreadSynchronize());
	
	//copy data from device to host
	HANDLE_ERROR(cudaMemcpy(x,   d_x, sizeof(double)*sz, cudaMemcpyDeviceToHost));
	//these arrays shouldn't need to be copied back
	//HANDLE_ERROR(cudaMemcpy(fc, d_fc, sizeof(double)*sz, cudaMemcpyDeviceToHost));
	//HANDLE_ERROR(cudaMemcpy(cc, d_cc, sizeof(double)*sz, cudaMemcpyDeviceToHost));
	//HANDLE_ERROR(cudaMemcpy(oC, d_oC, sizeof(double)*sz, cudaMemcpyDeviceToHost));
	//HANDLE_ERROR(cudaMemcpy(uC, d_uC, sizeof(double)*sz, cudaMemcpyDeviceToHost));
	//HANDLE_ERROR(cudaMemcpy(oN, d_oN, sizeof(double)*sz, cudaMemcpyDeviceToHost));
	//HANDLE_ERROR(cudaMemcpy(oE, d_oE, sizeof(double)*sz, cudaMemcpyDeviceToHost));

	//release cuda memory
	HANDLE_ERROR(cudaFree(d_x));  HANDLE_ERROR(cudaFree(d_cc)); HANDLE_ERROR(cudaFree(d_fc));
	HANDLE_ERROR(cudaFree(d_oC)); HANDLE_ERROR(cudaFree(d_uC)); HANDLE_ERROR(cudaFree(d_oE));
	HANDLE_ERROR(cudaFree(d_oN)); HANDLE_ERROR(cudaFree(d_x2)); 
	
	//reset the cuda device
	HANDLE_ERROR(cudaDeviceReset());
	
	if (*iresid == 1){
		Vmresid7_1s(nx, ny, nz, ipc, rpc, oC, cc, fc, oE, oN, uC, x, r);
	}

}

VPUBLIC void Vgsrb7x(int *nx,int *ny,int *nz,
        int *ipc, double *rpc,
        double *oC, double *cc, double *fc,
        double *oE, double *oN, double *uC,
        double *x, double *w1, double *w2, double *r,
        int *itmax, int *iters,
        double *errtol, double *omega,
        int *iresid, int *iadjoint) {

    int i, j, k, ioff;    
    
    //this macro creates variable dx_<arr>, dy_<arr>, and dz_<arr> with values nx, ny, and nz respectively
    MAT3(cc, *nx, *ny, *nz);
    MAT3(fc, *nx, *ny, *nz);
    MAT3( x, *nx, *ny, *nz);
    MAT3(w1, *nx, *ny, *nz);
    MAT3(w2, *nx, *ny, *nz);
    MAT3( r, *nx, *ny, *nz);
    MAT3(oE, *nx, *ny, *nz);
    MAT3(oN, *nx, *ny, *nz);
    MAT3(uC, *nx, *ny, *nz);
    MAT3(oC, *nx, *ny, *nz);

    
    for (*iters=1; *iters<=*itmax; (*iters)++) {
    		
        // Do the red points ***
        #pragma omp parallel for private(i, j, k, ioff)
        for (k=2; k<=*nz-1; k++) {
            for (j=2; j<=*ny-1; j++) {
                ioff = (1 - *iadjoint) * (    (j + k + 2) % 2)
                     + (    *iadjoint) * (1 - (j + k + 2) % 2);
                for (i=2+ioff; i<=*nx-1; i+=2) {
                    VAT3(x, i, j, k) = (
                            VAT3(fc,   i,  j,  k)
                         +  VAT3(oN,   i,   j,   k) * VAT3(x,   i, j+1,   k)
                         +  VAT3(oN,   i, j-1,   k) * VAT3(x,   i, j-1,   k)
                         +  VAT3(oE,   i,   j,   k) * VAT3(x, i+1,   j,   k)
                         +  VAT3(oE, i-1,   j,   k) * VAT3(x, i-1,   j,   k)
                         + VAT3( uC,   i,   j, k-1) * VAT3(x,   i,   j, k-1)
                         + VAT3( uC,   i,   j,   k) * VAT3(x,   i,   j, k+1)
                         ) / (VAT3(oC, i, j, k) + VAT3(cc, i, j, k));
                }
            }
        }

        // Do the black points
        #pragma omp parallel for private(i, j, k, ioff)
        for (k=2; k<=*nz-1; k++) {
            for (j=2; j<=*ny-1; j++) {
                ioff =   (    *iadjoint) * (    (j + k + 2) % 2 )
                       + (1 - *iadjoint) * (1 - (j + k + 2) % 2 );
                for (i=2+ioff;i<=*nx-1; i+=2) {
                    VAT3(x, i, j, k) = (
                            VAT3(fc,   i,   j,   k)
                         +  VAT3(oN,   i,   j,   k) * VAT3(x,   i,j+1,  k)
                         +  VAT3(oN,   i, j-1,   k) * VAT3(x,   i,j-1,  k)
                         +  VAT3(oE,   i,   j,   k) * VAT3(x, i+1,  j,  k)
                         +  VAT3(oE, i-1,   j,   k) * VAT3(x, i-1,  j,  k)
                         + VAT3( uC,   i,   j, k-1) * VAT3(x,   i,  j,k-1)
                         + VAT3( uC,   i,   j,   k) * VAT3(x,   i,  j,k+1)
                         ) / (VAT3(oC, i, j, k) + VAT3(cc, i, j, k));
                }
            }
        }
    }
 
    if (*iresid == 1){
        Vmresid7_1s(nx, ny, nz, ipc, rpc, oC, cc, fc, oE, oN, uC, x, r);
    }
}



VPUBLIC void Vgsrb27x(int *nx,int *ny,int *nz,
        int *ipc, double *rpc,
        double  *oC, double  *cc, double  *fc,
        double  *oE, double  *oN, double  *uC, double *oNE, double *oNW,
        double  *uE, double  *uW, double  *uN, double  *uS,
        double *uNE, double *uNW, double *uSE, double *uSW,
        double *x, double *w1, double *w2, double *r,
        int *itmax, int *iters,
        double *errtol, double *omega,
        int *iresid, int *iadjoint) {

    int  i,  j,  k;
    int i1, j1, k1;
    int i2, j2, k2;
    int ioff;
    int istep;

    double tmpO, tmpU, tmpD;

    MAT3( cc, *nx, *ny, *nz);
    MAT3(fc, *nx, *ny, *nz);
    MAT3( x, *nx, *ny, *nz);
    MAT3(w1, *nx, *ny, *nz);
    MAT3(w2, *nx, *ny, *nz);
    MAT3( r, *nx, *ny, *nz);

    MAT3(oE, *nx, *ny, *nz);
    MAT3(oN, *nx, *ny, *nz);
    MAT3(uC, *nx, *ny, *nz);
    MAT3(oC, *nx, *ny, *nz);

    MAT3(oNE, *nx, *ny, *nz);
    MAT3(oNW, *nx, *ny, *nz);

    MAT3( uE, *nx, *ny, *nz);
    MAT3( uW, *nx, *ny, *nz);
    MAT3( uN, *nx, *ny, *nz);
    MAT3( uS, *nx, *ny, *nz);
    MAT3(uNE, *nx, *ny, *nz);
    MAT3(uNW, *nx, *ny, *nz);
    MAT3(uSE, *nx, *ny, *nz);
    MAT3(uSW, *nx, *ny, *nz);

    // Do the gauss-seidel iteration itmax times

    /*
    i1    = (1 - *iadjoint) *   2  + (    *iadjoint) * (*nx - 1);
    i2    = (    *iadjoint) *   2  + (1 - *iadjoint) * (*nx - 1);
    j1    = (1 - *iadjoint) *   2  + (    *iadjoint) * (*ny - 1);
    j2    = (    *iadjoint) *   2  + (1 - *iadjoint) * (*ny - 1);
    k1    = (1 - *iadjoint) *   2  + (    *iadjoint) * (*nz - 1);
    k2    = (    *iadjoint) *   2  + (1 - *iadjoint) * (*nz - 1);
    istep = (    *iadjoint) * (-1) + (1 - *iadjoint) * (1);
    */

    i1 = (1-*iadjoint) * 2 + *iadjoint     * (*nx-1);
    i2 = *iadjoint     * 2 + (1-*iadjoint) * (*nx-1);
    j1 = (1-*iadjoint) * 2 + *iadjoint     * (*ny-1);
    j2 = *iadjoint     * 2 + (1-*iadjoint) * (*ny-1);
    k1 = (1-*iadjoint) * 2 + *iadjoint     * (*nz-1);
    k2 = *iadjoint     * 2 + (1-*iadjoint) * (*nz-1);
    istep = *iadjoint*(-1) + (1-*iadjoint)*(1);

    for (*iters=1; *iters<=*itmax; (*iters)++) {

        //#pragma omp parallel for private(i, j, k, ioff, tmpO, tmpU, tmpD)
        for (k=2; k<=*nz-1; k++) {

            for (j=2; j<=*ny-1; j++) {

                ioff = (1 - *iadjoint) * (    (j + k + 2) % 2)
                     + (    *iadjoint) * (1 - (j + k + 2) % 2);

                for (i=2+ioff; i<=*nx-1; i+=2) {

                    tmpO =
                         + VAT3(  oN,   i,   j,   k) * VAT3(x,   i, j+1,   k)
                         + VAT3(  oN,   i, j-1,   k) * VAT3(x,   i, j-1,   k)
                         + VAT3(  oE,   i,   j,   k) * VAT3(x, i+1,   j,   k)
                         + VAT3(  oE, i-1,   j,   k) * VAT3(x, i-1,   j,   k)
                         + VAT3( oNE,   i,   j,   k) * VAT3(x, i+1, j+1,   k)
                         + VAT3( oNW,   i,   j,   k) * VAT3(x, i-1, j+1,   k)
                         + VAT3( oNW, i+1, j-1,   k) * VAT3(x, i+1, j-1,   k)
                         + VAT3( oNE, i-1, j-1,   k) * VAT3(x, i-1, j-1,   k);

                   tmpU =
                         + VAT3(  uC,   i,   j,   k) * VAT3(x,   i,   j, k+1)
                         + VAT3(  uN,   i,   j,   k) * VAT3(x,   i, j+1, k+1)
                         + VAT3(  uS,   i,   j,   k) * VAT3(x,   i, j-1, k+1)
                         + VAT3(  uE,   i,   j,   k) * VAT3(x, i+1,   j, k+1)
                         + VAT3(  uW,   i,   j,   k) * VAT3(x, i-1,   j, k+1)
                         + VAT3( uNE,   i,   j,   k) * VAT3(x, i+1, j+1, k+1)
                         + VAT3( uNW,   i,   j,   k) * VAT3(x, i-1, j+1, k+1)
                         + VAT3( uSE,   i,   j,   k) * VAT3(x, i+1, j-1, k+1)
                         + VAT3( uSW,   i,   j,   k) * VAT3(x, i-1, j-1, k+1);

                   tmpD =
                         + VAT3(  uC,   i,   j, k-1) * VAT3(x,   i,   j, k-1)
                         + VAT3(  uS,   i, j+1, k-1) * VAT3(x,   i, j+1, k-1)
                         + VAT3(  uN,   i, j-1, k-1) * VAT3(x,   i, j-1, k-1)
                         + VAT3(  uW, i+1,   j, k-1) * VAT3(x, i+1,   j, k-1)
                         + VAT3(  uE, i-1,   j, k-1) * VAT3(x, i-1,   j, k-1)
                         + VAT3( uSW, i+1, j+1, k-1) * VAT3(x, i+1, j+1, k-1)
                         + VAT3( uSE, i-1, j+1, k-1) * VAT3(x, i-1, j+1, k-1)
                         + VAT3( uNW, i+1, j-1, k-1) * VAT3(x, i+1, j-1, k-1)
                         + VAT3( uNE, i-1, j-1, k-1) * VAT3(x, i-1, j-1, k-1);

                       VAT3(x, i,j,k) = (VAT3(fc, i, j, k) + (tmpO + tmpU + tmpD))
                                / (VAT3(oC, i, j, k) + VAT3(cc, i, j, k));

                }
            }
        }

        //#pragma omp parallel for private(i, j, k, ioff, tmpO, tmpU, tmpD)
        for (k=2; k<=*nz-1; k++) {

            for (j=2; j<=*ny-1; j++) {

                ioff = (    *iadjoint) * (    (j + k + 2) % 2)
                     + (1 - *iadjoint) * (1 - (j + k + 2) % 2);

                for (i=2+ioff; i<=*nx-1; i+=2) {

                    tmpO =
                         + VAT3(  oN,   i,   j,   k) * VAT3(x,   i, j+1,   k)
                         + VAT3(  oN,   i, j-1,   k) * VAT3(x,   i, j-1,   k)
                         + VAT3(  oE,   i,   j,   k) * VAT3(x, i+1,   j,   k)
                         + VAT3(  oE, i-1,   j,   k) * VAT3(x, i-1,   j,   k)
                         + VAT3( oNE,   i,   j,   k) * VAT3(x, i+1, j+1,   k)
                         + VAT3( oNW,   i,   j,   k) * VAT3(x, i-1, j+1,   k)
                         + VAT3( oNW, i+1, j-1,   k) * VAT3(x, i+1, j-1,   k)
                         + VAT3( oNE, i-1, j-1,   k) * VAT3(x, i-1, j-1,   k);

                    tmpU =
                         + VAT3(  uC,   i,   j,   k) * VAT3(x,   i,   j, k+1)
                         + VAT3(  uN,   i,   j,   k) * VAT3(x,   i, j+1, k+1)
                         + VAT3(  uS,   i,   j,   k) * VAT3(x,   i, j-1, k+1)
                         + VAT3(  uE,   i,   j,   k) * VAT3(x, i+1,   j, k+1)
                         + VAT3(  uW,   i,   j,   k) * VAT3(x, i-1,   j, k+1)
                         + VAT3( uNE,   i,   j,   k) * VAT3(x, i+1, j+1, k+1)
                         + VAT3( uNW,   i,   j,   k) * VAT3(x, i-1, j+1, k+1)
                         + VAT3( uSE,   i,   j,   k) * VAT3(x, i+1, j-1, k+1)
                         + VAT3( uSW,   i,   j,   k) * VAT3(x, i-1, j-1, k+1);

                   tmpD =
                         + VAT3(  uC,   i,   j, k-1) * VAT3(x,   i,   j, k-1)
                         + VAT3(  uS,   i, j+1, k-1) * VAT3(x,   i, j+1, k-1)
                         + VAT3(  uN,   i, j-1, k-1) * VAT3(x,   i, j-1, k-1)
                         + VAT3(  uW, i+1,   j, k-1) * VAT3(x, i+1,   j, k-1)
                         + VAT3(  uE, i-1,   j, k-1) * VAT3(x, i-1,   j, k-1)
                         + VAT3( uSW, i+1, j+1, k-1) * VAT3(x, i+1, j+1, k-1)
                         + VAT3( uSE, i-1, j+1, k-1) * VAT3(x, i-1, j+1, k-1)
                         + VAT3( uNW, i+1, j-1, k-1) * VAT3(x, i+1, j-1, k-1)
                         + VAT3( uNE, i-1, j-1, k-1) * VAT3(x, i-1, j-1, k-1);

                   VAT3(x, i,j,k) = (VAT3(fc, i, j, k) + (tmpO + tmpU + tmpD))
                            / (VAT3(oC, i, j, k) + VAT3(cc, i, j, k));
                }
            }
        }
    }

    // If specified, return the new residual as well
    if (*iresid == 1)
        Vmresid27_1s(nx, ny, nz,
                     ipc, rpc,
                      oC,  cc,  fc,
                      oE,  oN,  uC,
                     oNE, oNW,
                     uE,   uW,  uN,  uS,
                     uNE, uNW, uSE, uSW,
                       x,   r);
}
