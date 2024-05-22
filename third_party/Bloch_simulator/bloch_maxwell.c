
#include "mex.h" 
#include <stdio.h>

#define _USE_MATH_DEFINES /* To use MATH constants */
#include <math.h>

#define TWOPI  2*M_PI
#define gamma  42577.46778     /* gyromagnetic ratio [Hz/mT] */
#define GAMMA  gamma/10*TWOPI  /* gyromagnetic ratio [rad/sec/G] */

//#define GAMMA   26753.0 /* [rad/sec/G] */
//#define TWOPI   6.283185

//#define DEBUG


/* Multiply 3x3 matrix by 3x1 vector */
void multmatvec(double *mat, double *vec, double *matvec)
{
    *matvec++ = mat[0] * vec[0] + mat[3] * vec[1] + mat[6] * vec[2];
    *matvec++ = mat[1] * vec[0] + mat[4] * vec[1] + mat[7] * vec[2];
    *matvec++ = mat[2] * vec[0] + mat[5] * vec[1] + mat[8] * vec[2];
}


/* Add two 3x1 vectors */
void addvecs(double *vec1, double *vec2, double *vecsum)
{
    *vecsum++ = *vec1++ + *vec2++;
    *vecsum++ = *vec1++ + *vec2++;
    *vecsum++ = *vec1++ + *vec2++;
}


/* ======== Adjoint of a 3x3 matrix ========= */
void adjmat(double *mat, double *adj)
{
    *adj++ =  (mat[4] * mat[8] - mat[7] * mat[5]);	
    *adj++ = -(mat[1] * mat[8] - mat[7] * mat[2]);
    *adj++ =  (mat[1] * mat[5] - mat[4] * mat[2]);
    *adj++ = -(mat[3] * mat[8] - mat[6] * mat[5]);
    *adj++ =  (mat[0] * mat[8] - mat[6] * mat[2]);
    *adj++ = -(mat[0] * mat[5] - mat[3] * mat[2]);
    *adj++ =  (mat[3] * mat[7] - mat[6] * mat[4]);
    *adj++ = -(mat[0]  *mat[7] - mat[6] * mat[1]);
    *adj++ =  (mat[0]  *mat[4] - mat[3] * mat[1]);
}


/* ====== Set a 3x3 matrix to all zeros	======= */
void zeromat(double *mat)
{
    *mat++ = 0;
    *mat++ = 0;
    *mat++ = 0;
    *mat++ = 0;
    *mat++ = 0;
    *mat++ = 0;
    *mat++ = 0;
    *mat++ = 0;
    *mat++ = 0;
}


/* ======== Return 3x3 Identity Matrix  ========= */
void eyemat(double *mat)
{
    zeromat(mat);
    mat[0] = 1;
    mat[4] = 1;
    mat[8] = 1;
}


/* ======== Determinant of a 3x3 matrix ======== */
double detmat(double *mat)
{
    double det;

    det = mat[0] * mat[4] * mat[8];
    det+= mat[3] * mat[7] * mat[2];
    det+= mat[6] * mat[1] * mat[5];
    det-= mat[0] * mat[7] * mat[5];
    det-= mat[3] * mat[1] * mat[8];
    det-= mat[6] * mat[4] * mat[2];

    return det;
}


/* ======== multiply a matrix by a scalar ========= */
void scalemat(double *mat, double scalar)
{
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
}


/* ======== Inverse of a 3x3 matrix ========= */
/*	DO NOT MAKE THE OUTPUT THE SAME AS ONE OF THE INPUTS!! */
void invmat(double *mat, double *imat)
{
    double det;
    int count;

    det = detmat(mat);  /* Determinant */
    adjmat(mat, imat);  /* Adjoint */

    for (count = 0; count < 9; count++)
    {
        *imat = *imat++ / det;
    }
}


/* ====== Add two 3x3 matrices	====== */
void addmats(double *mat1, double *mat2, double *matsum)
{
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
}


/* ======= Multiply two 3x3 matrices. ====== */
/*	DO NOT MAKE THE OUTPUT THE SAME AS ONE OF THE INPUTS!! */
void multmats(double *mat1, double *mat2, double *matproduct)
{
    *matproduct++ = mat1[0] * mat2[0] + mat1[3] * mat2[1] + mat1[6] * mat2[2];
    *matproduct++ = mat1[1] * mat2[0] + mat1[4] * mat2[1] + mat1[7] * mat2[2];
    *matproduct++ = mat1[2] * mat2[0] + mat1[5] * mat2[1] + mat1[8] * mat2[2];
    *matproduct++ = mat1[0] * mat2[3] + mat1[3] * mat2[4] + mat1[6] * mat2[5];
    *matproduct++ = mat1[1] * mat2[3] + mat1[4] * mat2[4] + mat1[7] * mat2[5];
    *matproduct++ = mat1[2] * mat2[3] + mat1[5] * mat2[4] + mat1[8] * mat2[5];
    *matproduct++ = mat1[0] * mat2[6] + mat1[3] * mat2[7] + mat1[6] * mat2[8];
    *matproduct++ = mat1[1] * mat2[6] + mat1[4] * mat2[7] + mat1[7] * mat2[8];
    *matproduct++ = mat1[2] * mat2[6] + mat1[5] * mat2[7] + mat1[8] * mat2[8];
}


/* Find the rotation matrix that rotates |n| radians about */
/* the vector given by nx,ny,nz                            */
void calcrotmat(double nx, double ny, double nz, double *rmat)
{
    double ar, ai, br, bi, hp, cp, sp;
    double arar, aiai, arai2, brbr, bibi, brbi2, arbi2, aibr2, arbr2, aibi2;
    double phi;

    phi = sqrt(nx*nx + ny*ny + nz*nz);
    /* printf("calcrotmat(%6.3f,%6.3f,%6.3f) -> phi = %6.3f\n",nx,ny,nz,phi); */

    if (phi == 0.0)
	{
        *rmat++ = 1;
        *rmat++	= 0;
        *rmat++ = 0;
        *rmat++ = 0;
        *rmat++ = 1;
        *rmat++	= 0;
        *rmat++ = 0;
        *rmat++ = 0;
        *rmat++ = 1;
    } else
	{
        /* First define Cayley-Klein parameters */
        hp = phi/2;
        cp = cos(hp);
        sp = sin(hp)/phi; /* /phi because n is unit length in defs. */
        ar = cp;
        ai = -nz*sp;
        br = ny*sp;
        bi = -nx*sp;

        /* Make auxiliary variables to speed this up */
        arar = ar*ar;
        aiai = ai*ai;
        arai2 = 2*ar*ai;
        brbr = br*br;
        bibi = bi*bi;
        brbi2 = 2*br*bi;
        arbi2 = 2*ar*bi;
        aibr2 = 2*ai*br;
        arbr2 = 2*ar*br;
        aibi2 = 2*ai*bi;
        
        /* Make rotation matrix */
        *rmat++ = arar - aiai - brbr + bibi;
        *rmat++ = -arai2 - brbi2;
        *rmat++ = -arbr2 + aibi2;
        *rmat++ =  arai2 - brbi2;
        *rmat++ = arar - aiai + brbr - bibi;
        *rmat++ = -aibr2 - arbi2;
        *rmat++ =  arbr2 + aibi2;
        *rmat++ =  arbi2 - aibr2;
        *rmat++ = arar + aiai - brbr - bibi;
    }
}


/*	Set a 3x1 vector to all zeros */
void zerovec(double *vec)
{
    *vec++ = 0;
    *vec++ = 0;
    *vec++ = 0;
}

/* ------------------------------------------------------------
	Function takes the given endtimes of intervals, and
	returns the interval lengths in an array, assuming that
	the first interval starts at 0.

	If the intervals are all greater than 0, then this
	returns 1, otherwise it returns 0.
   ------------------------------------------------------------ */
int times2intervals( double *endtimes, double *intervals, long n)
{
    int allpos;
    int count;
    double lasttime;

    allpos = 1;
    lasttime = 0.0;

    for (count = 0; count < n; count++)
	{
        intervals[count] = endtimes[count]-lasttime;
        lasttime = endtimes[count];
        if (intervals[count] <= 0)
        {
            allpos = 0;
        }
	}

    return (allpos);
}

/* Go through time for one df and one dx, dy, dz */
void blochsim(double *b1real, double *b1imag,
              double *xgrad, double *ygrad, double *zgrad, double *tsteps,
              int ntime, double *e1, double *e2, double df,
              double dx, double dy, double dz,
              double *mx, double *my, double *mz, int mode, double B0, double alpha, double g) /* NGL: Maxwell */
{
    int tcount;
    double gammadx;
    double gammady;
    double gammadz;
    double dt;
    double rotmat[9];
    double amat[9], bvec[3];  /* A and B propagation matrix and vector */
    double arot[9], brot[3];  /* A and B after rotation step */
    double decmat[9];         /* Decay matrix for each time step */
    double decvec[3];         /* Recovery vector for each time step */
    double rotx,roty,rotz;    /* Rotation axis coordinates */
    double imat[9], mvec[3];
    double mcurr0[3];        /* Current magnetization before rotation */
    double mcurr1[3];		 /* Current magnetization before decay */

    /* NGL: Maxwell -- start 05/03/2020 */
    double alpha2;
    double g2;
    double Gx;
    double Gy;
    double Gz;
    double Gx2;
    double Gy2;
    double Gz2;
    double B02;
    double Gz3;
    double x2;
    double y2;
    double z2;
    double xy;
    double xz;
    double yz;
    double x2_term;
    double y2_term;
    double z2_term;
    double xy_term;
    double yz_term;
    double xz_term;
    double x3;
    double y3;
    double z3;
    double x2y;
    double x2z;
    double xy2;
    double y2z;
    double xz2;
    double yz2;
    double xyz;
    double x3_term;
    double y3_term;
    double z3_term;
    double x2y_term;
    double x2z_term;
    double xy2_term;
    double y2z_term;
    double xz2_term;
    double yz2_term;
    double xyz_term;
    double sum_concomitant_terms;
    /* NGL: Maxwell -- end 05/03/2020 */

    eyemat(amat);            /* A is the identity matrix */
    eyemat(imat);            /* I is the identity matrix */

    zerovec(bvec);
    zerovec(decvec);
    zeromat(decmat);

    gammadx = dx * GAMMA;  /* Convert to Hz/cm(?), NGL => [rad/sec/G] * [cm] */
    gammady = dy * GAMMA;  /* Convert to Hz/cm */
    gammadz = dz * GAMMA;  /* Convert to Hz/cm */

    mcurr0[0] = *mx;  /* Set starting x magnetization */
    mcurr0[1] = *my;  /* Set starting y magnetization */
    mcurr0[2] = *mz;  /* Set starting z magnetization */

    /* NGL: Maxwell -- start 05/03/2020 */
    g = g * 1e-4; /* [G/cm] * [T/1e4G] => [T/cm] */
    alpha2 = alpha * alpha;
    g2 = g * g; /* [T/cm]^2 */
    /* NGL: Maxwell -- end 05/03/2020 */
    
    for (tcount = 0; tcount < ntime; tcount++)
	{
        /* NGL: Maxwell -- start 05/03/2020 */
        Gx = xgrad[tcount] * 1e-4; /* [G/cm] * [T/1e4G] => [T/cm] */
        Gy = ygrad[tcount] * 1e-4; /* [G/cm] * [T/1e4G] => [T/cm]  */
        Gz = zgrad[tcount] * 1e-4; /* [G/cm] * [T/1e4G] => [T/cm]  */
        
        /*---------------------------------------------------------------*/
        /* Calculate 1/B0 order concomitant gradient terms               */
        /* (x2, y2, z2, xy, yz, xz)                                      */
        /* 1 / [T] * (([G/cm] * [T/1e4G])^2 * [cm]^2) => [T]             */
        /* 1 / [T] * ([G/cm]^2 * [T/1e4G]^2 * [cm]^2) => [T]             */
        /*---------------------------------------------------------------*/
        Gx2 = Gx * Gx; /* [T/cm]^2 */
        Gy2 = Gy * Gy; /* [T/cm]^2 */
        Gz2 = Gz * Gz; /* [T/cm]^2 */

        x2 = dx * dx; /* [cm]^2 */
        y2 = dy * dy; /* [cm]^2 */
        z2 = dz * dz; /* [cm]^2 */
        xy = dx * dy; /* [cm]^2 */
        xz = dx * dz; /* [cm]^2 */
        yz = dy * dz; /* [cm]^2 */

        x2_term = (1 / (2 * B0) * (alpha2 * Gz2 + g2)) * x2;
        y2_term = (1 / (2 * B0) * ((alpha - 1) * (alpha - 1) * Gz2 + g2)) * y2;
        z2_term = (1 / (2 * B0) * (Gx2 + Gy2)) * z2;
        xy_term = (1 / B0 * (-g * Gz)) * xy;
        yz_term = (1 / B0 * (g * Gx + (alpha - 1) * Gy * Gz)) * yz;
        xz_term = (1 / B0 * (g * Gy - alpha * Gx * Gz)) * xz;

        /*---------------------------------------------------------------*/
        /* Calcualte 1/B0^2 order concomitant gradient terms (10 terms)  */
        /* (x3, y3, z3, x2y, x2z, xy2, y2z, xz2, yz2, xyz)               */
        /* 1 / [T]^2 * (([G/cm] * [T/1e4G])^3 * [cm]^3) => [T]           */
        /* 1 / [T]^2 * ([G/cm]^3 * [T/1e4G]^3 * [cm]^2) => [T]           */
        /*---------------------------------------------------------------*/
        B02 = B0 * B0; /* [T]^2 */
        Gz3 = Gz * Gz * Gz; /* [T/cm]^3 */

        x3  = dx * dx * dx; /* [cm]^3 */
        y3  = dy * dy * dy; /* [cm]^3 */
        z3  = dz * dz * dz; /* [cm]^3 */
        x2y = dx * dx * dy; /* [cm]^3 */
        x2z = dx * dx * dz; /* [cm]^3 */
        xy2 = dx * dy * dy; /* [cm]^3 */
        y2z = dy * dy * dz; /* [cm]^3 */
        xz2 = dx * dz * dz; /* [cm]^3 */
        yz2 = dy * dz * dz; /* [cm]^3 */
        xyz = dx * dy * dz; /* [cm]^3 */

        x3_term  = (-1 / (2 * B02) * (Gx * (alpha2 * Gz2 + g2))) * x3;
        y3_term  = (-1 / (2 * B02) * (Gy * ((alpha - 1) * (alpha - 1) * Gz2 + g2))) * y3;
        z3_term  = (-1 / (2 * B02) * (Gz * (Gx2 + Gy2))) * z3;
        x2y_term = (-1 / (2 * B02) * (Gy * (alpha2 * Gz2 + g2) - 2 * g * Gx * Gy)) * x2y;
        x2z_term = (-1 / (2 * B02) * (Gz * (alpha2 * Gz2 + g2) + 2 * Gx * (g * Gy - alpha * Gx * Gz))) * x2z;
        xy2_term = (-1 / (2 * B02) * (Gx * (alpha - 1) * (alpha - 1) * Gz2 + g2 * Gx - 2 * g * Gy * Gz)) * xy2;
        y2z_term = (-1 / (2 * B02) * ((alpha - 1) * (alpha - 1) * Gz3 + g2 * Gz + 2 * Gy * (g * Gx + (alpha - 1) * Gy * Gz))) * y2z;
        xz2_term = (-1 / (2 * B02) * (Gx * (Gx2 + Gy2) + 2 * Gz * (g * Gy - alpha * Gx * Gz))) * xz2;
        yz2_term = (-1 / (2 * B02) * (Gy * (Gx2 + Gy2) + 2 * Gz * (g * Gx + (alpha - 1) * Gy * Gz))) * yz2;
        xyz_term = (-1 / B02 * (g * (Gx2 + Gy2 - Gz2) - (Gx * Gy * Gz))) * xyz;

        /*---------------------------------------------------------------*/
        /* Calculate the sum of all concomitant terms                    */ 
        /*---------------------------------------------------------------*/
        sum_concomitant_terms = x2_term + y2_term + z2_term + xy_term + yz_term + xz_term + 
                                x3_term + y3_term + z3_term + x2y_term + x2z_term + xy2_term + y2z_term + xz2_term + yz2_term + xyz_term; /* [T] */
        /* NGL: Maxwell -- end 05/03/2020 */

        /* Rotation */
        dt = tsteps[tcount];
        /* [G/cm] * [rad/sec/G] * [cm] => [rad/sec] */ 
        /* [rad/sec/G] * [1e4G/T] * [T] => [rad/sec] */
        rotz = -(xgrad[tcount] * gammadx + ygrad[tcount] * gammady + zgrad[tcount] * gammadz + df * TWOPI + (GAMMA * sum_concomitant_terms * 1e4)) * dt;
        rotx = (- b1real[tcount] * GAMMA * dt);
        //roty = (+ *b1imag++ * GAMMA * *tsteps++);
        roty = (- b1imag[tcount] * GAMMA * dt); /* NGL: changed to minus */
        calcrotmat(rotx, roty, rotz, rotmat);

        if (mode == 1)
        {
            multmats(rotmat, amat, arot);
            multmatvec(rotmat, bvec, brot);
		} else
        {
            multmatvec(rotmat, mcurr0, mcurr1);
        }

        /* Decay */
        decvec[2] = 1 - *e1;
        decmat[0] = *e2;
        decmat[4] = *e2++;
        decmat[8] = *e1++;

        if (mode == 1)
        {
            multmats(decmat, arot, amat);
            multmatvec(decmat, brot, bvec);
            addvecs(bvec, decvec, bvec);
		} else
        {
            multmatvec(decmat, mcurr1, mcurr0);
            addvecs(mcurr0, decvec, mcurr0);
        }

        /* printf("rotmat = [%6.3f  %6.3f  %6.3f ] \n", rotmat[0], rotmat[3], rotmat[6]);
         * printf("         [%6.3f  %6.3f  %6.3f ] \n", rotmat[1], rotmat[4], rotmat[7]);
         * printf("         [%6.3f  %6.3f  %6.3f ] \n", rotmat[2], rotmat[5], rotmat[8]);
         * printf("A = [%6.3f  %6.3f  %6.3f ] \n", amat[0], amat[3], amat[6]);
         * printf("    [%6.3f  %6.3f  %6.3f ] \n", amat[1], amat[4], amat[7]);
         * printf("    [%6.3f  %6.3f  %6.3f ] \n", amat[2], amat[5], amat[8]);
         * printf(" B = <%6.3f,%6.3f,%6.3f> \n", bvec[0], bvec[1], bvec[2]);
         * printf("<mx,my,mz> = <%6.3f,%6.3f,%6.3f> \n", amat[6] + bvec[0], amat[7] + bvec[1], amat[8] + bvec[2]);
         * printf("\n");
         */

        /* Sample output at times.  */
        /* Only do this if transient! */
        if (mode == 2)
        {
            *mx = mcurr0[0];
            *my = mcurr0[1];
            *mz = mcurr0[2];
            mx++;
            my++;
            mz++;
        }
    }

    /* If only recording the endpoint, either store the last point, */
    /* or calculate the steady-state endpoint.                      */
    if (mode == 0) /* Indicates start at given m, or m0. */
    {
        *mx = mcurr0[0];
        *my = mcurr0[1];
        *mz = mcurr0[2];
    } else if (mode == 1) /* Indicates to find steady-state magnetization */
    {
        scalemat(amat,-1.0);          /* Negate A matrix */
        addmats(amat, imat, amat);    /* Now amat = (I-A) */
        invmat(amat,imat);            /* Reuse imat as inv(I-A) */
        multmatvec(imat, bvec, mvec); /* Now M = inv(I-A)*B */
        *mx = mvec[0];
        *my = mvec[1];
        *mz = mvec[2];
    }
}



void blochsimfz(double *b1real, double *b1imag, double *xgrad, double *ygrad, double *zgrad,
                double *tsteps, int ntime, double t1, double t2, double *dfreq, int nfreq,
                double *dxpos, double *dypos, double *dzpos, int npos,
                double *mx, double *my, double *mz, int mode, double B0, double alpha, double g)
{
    int count;
    int poscount;
    int fcount;
    int totpoints;
    int totcount = 0;

    int ntout;

    double *e1;
    double *e2;
    double *e1ptr;
    double *e2ptr;
    double *tstepsptr;
    double *dxptr, *dyptr, *dzptr;

    if (mode & 2)
    {
        ntout = ntime;
    } else
    {
        ntout = 1;
    }

    /* First calculate the E1 and E2 values at each time step */
    e1 = (double *) malloc(ntime * sizeof(double));
    e2 = (double *) malloc(ntime * sizeof(double));

    e1ptr = e1;
    e2ptr = e2;
    tstepsptr = tsteps;

    for (count = 0; count < ntime; count++)
    {
        *e1ptr++ = exp(- *tstepsptr / t1);
        *e2ptr++ = exp(- *tstepsptr++ / t2);
    }

    totpoints = npos * nfreq;

    for (fcount = 0; fcount < nfreq; fcount++)
    {
        dxptr = dxpos;
        dyptr = dypos;
        dzptr = dzpos;
        for (poscount = 0; poscount < npos; poscount++)
        {
            if (mode == 3)	/* Steady state AND record all time points */
            {	/* First go through and find steady state, then
				repeat as if transient starting at steady st.*/
                blochsim(b1real, b1imag, xgrad, ygrad, zgrad, tsteps, ntime,
                         e1, e2, *dfreq, *dxptr, *dyptr, *dzptr, mx, my, mz, 1, B0, alpha, g);

                blochsim(b1real, b1imag, xgrad, ygrad, zgrad, tsteps, ntime,
                         e1, e2, *dfreq, *dxptr++, *dyptr++, *dzptr++, mx, my, mz, 2, B0, alpha, g);
            } else
            {
                blochsim(b1real, b1imag, xgrad, ygrad, zgrad, tsteps, ntime,
                         e1, e2, *dfreq, *dxptr++, *dyptr++, *dzptr++, mx, my, mz, mode, B0, alpha, g);
            }

            mx += ntout;
            my += ntout;
            mz += ntout;

            totcount++;
            if ((totpoints > 40000) && ( ((10*totcount)/totpoints)> (10*(totcount-1)/totpoints) ))
            {
                printf("%d%% Complete.\n",(100*totcount/totpoints));
            }
        }
        dfreq++;
    }

    free(e1);
    free(e2);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* bloch(b1, gradxyz, dt, t1, t2, df, dx, dy, dz, mode, B0, alpha, g, mx, my, mz) */
    double *b1r;   /* Real-part of B1 field [G]         */
    double *b1i;   /* Imag-part of B1 field [G]         */
    double *gx;    /* X-axis gradient [G/cm]            */
    double *gy;    /* Y-axis gradient [G/cm]            */
    double *gz;    /* Z-axis gradient [G/cm]            */
    double *tp;    /* Time steps [sec]                  */
    double *ti;    /* Time intervals [sec]              */
    double t1;     /* T1 time constant [sec]            */
    double t2;     /* T2 time constant [sec]            */
    double *df;    /* Off-resonance frequencies [Hz]    */
    double *dx;    /* X Positions [cm]                  */
    double *dy;    /* Y Positions [cm]                  */
    double *dz;    /* Z Positions [cm]                  */
    int md;        /* Mode - 0=from M0, 1=steady-state	*/
    double tstep;  /* Time step, if single parameter    */
    double B0;     /* main magnetic field strength      */
    double alpha;  /* dimensionless symmetry parameter  */
    double g;      /* transverse gradient               */

    double *mxin;  /* Input points */
    double *myin;
    double *mzin;

    double *mxout;  /* Input points */
    double *myout;
    double *mzout;

    double *mx;    /* Output Arrays */
    double *my;
    double *mz;

    int gyaflag = 0;  /* 1 if gy was allocated */
    int gzaflag = 0;  /* 1 if gy was allocated */
    int dyaflag = 0;  /* 1 if dy was allocated */
    int dzaflag = 0;  /* 1 if dy was allocated */

    int ntime;      /* Number of time points */
    int ntout;      /* Number of time poitns at output */
    mwSize outsize[3]; /* Output matrix sizes */

    int ngrad;    /* Number of gradient dimensions */
    int nf;
    int npos;     /* Number of positions.  Calculated from nposN and nposM, depends on them */
    int nposM;    /* Height of passed position matrix */
    int nposN;    /* Width of passed position matrix  */
    int nfnpos;   /* Number of frequencies * number of positions */
    int ntnfnpos; /* Number of output times *frequencies*number of positions */
    int count;


    #ifdef DEBUG
    printf("---------------------------------------\n");
    printf("3D-position + frequency Bloch Simulator\n");
    printf("---------------------------------------\n\n");
    #endif

    /* Number of Time, RF, and Grad points */
    ntime = (int) mxGetM(prhs[0]) * (int) mxGetN(prhs[0]);

    /* ====================== RF (B1) =========================
    /* :  If complex, split up.  If real, allocate an imaginary part ==== */
    if (mxIsComplex(prhs[0]))
    {
        b1r = mxGetPr(prhs[0]);
        b1i = mxGetPi(prhs[0]);
    } else
    {
        b1r = mxGetPr(prhs[0]);
        b1i = (double *) malloc(ntime * sizeof(double));
        for (count = 0; count < ntime; count++)
        {
            b1i[count] = 0.0;
        }
    }
    #ifdef DEBUG
    printf("%d B1 points.\n", ntime);
    #endif

    /* ======================= Gradients ========================= */
    /* Number of Time, RF, and Grad points */
    ngrad = (int) mxGetM(prhs[1]) * (int) mxGetN(prhs[1]);
    gx = mxGetPr(prhs[1]); /* X-gradient is first N points */

    if (ngrad < 2*ntime) /* Need to allocate Y-gradient */
	{
        #ifdef DEBUG
        printf("Assuming 1-Dimensional Gradient\n");
        #endif
        gy = (double *) malloc(ntime * sizeof(double));
        gyaflag = 1;
        for (count = 0; count < ntime; count++)
        {
            gy[count] = 0.0;
        }
    } else
    {
        #ifdef DEBUG
        printf("Assuming (at least) 2-Dimensional Gradient\n");
	    #endif
        gy = gx + ntime; /* Assign from Nx3 input array */
	}

    if (ngrad < 3*ntime) /* Need to allocate Z-gradient */
    {
        gz = (double *) malloc(ntime * sizeof(double));
        gzaflag = 1;
        for (count = 0; count < ntime; count++)
        {
            gz[count] = 0.0;
        }
    } else
    {
        #ifdef DEBUG
        printf("Assuming 3-Dimensional Gradient\n");
        #endif
        gz = gx + 2*ntime; /* Assign from Nx3 input array */
    }

    /* Warning if Gradient length is not 1x, 2x, or 3x RF length */
    #ifdef DEBUG
    printf("%d Gradient Points (total) \n", ngrad);
    #endif
    if ((ngrad != ntime) && (ngrad != 2*ntime) && (ngrad != 3*ntime))
    {
        printf("Gradient length differs from B1 length\n");
    }

    if (gx == NULL) 
    {
        printf("ERROR:  gx is not allocated. \n");
    }
    if (gy == NULL)
    {
        printf("ERROR:  gy is not allocated. \n");
    }
    if (gz == NULL)
    {
        printf("ERROR:  gz is not allocated. \n");
    }

    /* === Time points ===== */
    /*	THREE Cases:
		1) Single value given -> this is the interval length for all.
		2) List of intervals given.
		3) Monotonically INCREASING list of end times given.
     *
     * For all cases, the goal is for tp to have the intervals.
     */

    ti = NULL;
    tp = mxGetPr(prhs[2]);
    if (mxGetM(prhs[2]) * mxGetN(prhs[2]) == 1)	/* === Case 1 === */
    {
        tp = (double *) malloc(ntime * sizeof(double));
        tstep = *(mxGetPr(prhs[2]));
        for (count = 0; count < ntime; count++)
        {
            tp[count] = tstep;
        }
    } else if (mxGetM(prhs[2]) * mxGetN(prhs[2]) != ntime)
    {
        printf("Time-point length differs from B1 length\n");
    } else
    {
        tp = mxGetPr(prhs[2]);
        ti = (double *) malloc(ntime * sizeof(double));
        if (( times2intervals(tp, ti, ntime)))
        {
            printf("Times are monotonically increasing. \n");
            tp = ti;
        }
    }

    /* === Relaxation Times ===== */
    t1 = *mxGetPr(prhs[3]);
    t2 = *mxGetPr(prhs[4]);

    /* === Frequency Points ===== */
    df = mxGetPr(prhs[5]);
    nf = (int) mxGetM(prhs[5]) * (int) mxGetN(prhs[5]);

    #ifdef DEBUG
    printf("%d Frequency points.\n", nf);
    #endif

    /* === Position Points ===== */
    nposM = (int) mxGetM(prhs[6]);
    nposN = (int) mxGetN(prhs[6]);

    #ifdef DEBUG
    printf("Position vector is %d x %d.\n", nposM, nposN);
    #endif

    if (nposN == 3) /* Assume 3 position dimensions given */
    {
        npos = nposM;
        #ifdef DEBUG
        printf("Assuming %d 3-Dimensional Positions\n", npos);
        #endif
        dx = mxGetPr(prhs[6]);
        dy = dx + npos;
        dz = dy + npos;
    } else if (nposN == 2) /* Assume only 2 position dimensions given */
    {
        npos = nposM;
        #ifdef DEBUG
        printf("Assuming %d 2-Dimensional Positions\n", npos);
        #endif
        dx = mxGetPr(prhs[6]);
        dy = dx + npos;
        dz = (double *) malloc(npos * sizeof(double));
        dzaflag = 1;
        for (count = 0; count < npos; count++)
        {
            dz[count] = 0.0;
        }
    } else /* Either 1xN, Nx1 or something random.  In all these
           /* cases we assume that 1 position is given, because it
           /* is too much work to try to figure out anything else! */
    {
        npos = nposM * nposN;
        #ifdef DEBUG
        printf("Assuming %d 1-Dimensional Positions\n", npos);
        #endif
        dx = mxGetPr(prhs[6]);
        dy = (double *) malloc(npos * sizeof(double));
        dz = (double *) malloc(npos * sizeof(double));
        dyaflag = 1;
        dzaflag = 1;
        for (count = 0; count < npos; count++)
        {
            dy[count] = 0.0;
            dz[count] = 0.0;
        }
        #ifdef DEBUG
        if ((nposM != 1) && (nposN != 1))
        {
            printf("Position vector should be 1xN, Nx1, Nx2 or Nx3. \n");
            printf(" -> Assuming 1 position dimension is given. \n");
        }
        #endif
    }

    if (dx == NULL)
    {
        printf("ERROR:  dx is not allocated. \n");
    }
    if (dy == NULL)
    {
        printf("ERROR:  dy is not allocated. \n");
    }
    if (dz == NULL)
    {
        printf("ERROR:  dz is not allocated. \n");
    }
    nfnpos = nf * npos; /* Just used to speed things up below */ 

    /* ===== Mode, defaults to 0 (simulate single endpoint, transient). ==== */
    if (nrhs > 7)
    {
        md = (int)(*mxGetPr(prhs[7]));
    } else
    {
        md = 0;
    }

    if (md & 2)
    {
        ntout = ntime; /* Include time points.	*/
    } else
    {
        ntout = 1;
    }

    #ifdef DEBUG
    printf("Mode = %d, %d Output Time Points \n", md, ntout);
    #endif

    ntnfnpos = ntout * nfnpos;

    #ifdef DEBUG
    if ((md & 1) == 0)
    {
        printf("Simulation from Initial Condition.\n");
    } else
    {
        printf("Simulation of Steady-State.\n");
    }    
    if ((md & 2) == 0)
    {
        printf("Simulation to Endpoint. \n");
    } else
    {
        printf("Simulation over Time.\n");
    }
    #endif

    /* NGL: Maxwell -- start 05/03/2020 */
    B0    = mxGetScalar(prhs[8]);
    alpha = mxGetScalar(prhs[9]);
    g     = mxGetScalar(prhs[10]);
    //printf("B0 = %g [T], alpha = %g, g = %d\n", B0, alpha, g);
    /* NGL: Maxwell -- end 05/03/2020 */

    /* ===== Allocate Output Magnetization vectors arrays */
    plhs[0] = mxCreateDoubleMatrix(ntnfnpos, 1, mxREAL); /* Mx, output */
    plhs[1] = mxCreateDoubleMatrix(ntnfnpos, 1, mxREAL); /* My, output */
    plhs[2] = mxCreateDoubleMatrix(ntnfnpos, 1, mxREAL); /* Mz, output */

    mx = mxGetPr(plhs[0]);
    my = mxGetPr(plhs[1]);
    mz = mxGetPr(plhs[2]);

    mxout = mx;
    myout = my;
    mzout = mz;

    /* ===== If Initial Magnetization is given... */
    if ((nrhs > 13) &&
        (mxGetM(prhs[11]) * mxGetN(prhs[11]) == nfnpos) &&
        (mxGetM(prhs[12]) * mxGetN(prhs[12]) == nfnpos) &&
        (mxGetM(prhs[13]) * mxGetN(prhs[13]) == nfnpos))

        /* Set output magnetization to that passed.
         * If multiple time points, then just the
         * first is set. */
    {
        #ifdef DEBUG
        printf("Using Specified Initial Magnetization.\n");
        #endif

        mxin = mxGetPr(prhs[11]);
        myin = mxGetPr(prhs[12]);
        mzin = mxGetPr(prhs[13]);
        for (count = 0; count < nfnpos; count++)
        {
            *mxout = *mxin++;
            *myout = *myin++;
            *mzout = *mzin++;
            mxout += ntout;
            myout += ntout;
            mzout += ntout;
        }
    } else
    {
        #ifdef DEBUG
        if (nrhs > 13) /* Magnetization given, but wrong size! */
        {
            printf("Initial magnetization passed, but not Npositions x Nfreq.\n");
        }
        printf(" --> Using [0; 0; 1] for initial magnetization. \n");
        #endif
        for (count = 0; count < nfnpos; count++)
        {
            *mxout = 0;	/* Set magnetization to Equilibrium */
            *myout = 0;
            *mzout = 1;
            mxout += ntout;
            myout += ntout;
            mzout += ntout;
        }
    }

    /* ======= Do The Simulation! ====== */
    #ifdef DEBUG
    printf("Calling blochsimfz() function in Mex file.\n");
    #endif
    blochsimfz(b1r, b1i, gx, gy, gz, tp, ntime, t1, t2, df, nf, dx, dy, dz, npos, mx, my, mz, md, B0, alpha, g);

    /* ======= Reshape Output Matrices ====== */
    if ((ntout > 1) && (nf > 1) && (npos > 1))
    {
        //outsize = (mwSize *) mxMalloc(3 * sizeof(*outsize)); /* NGL 2019-08-29 */
        outsize[0] = ntout;
        outsize[1] = npos;
        outsize[2] = nf;
        mxSetDimensions(plhs[0], outsize, 3); /* Set to 3D array */
        mxSetDimensions(plhs[1], outsize, 3); /* Set to 3D array */
        mxSetDimensions(plhs[2], outsize, 3); /* Set to 3D array */
    } else /* Basically "squeeze" the matrix */
    {
        if (ntout > 1)
        {
            //outsize = (mwSize *) mxMalloc(2 * sizeof(*outsize)); /* NGL 2019-08-29 */
            outsize[0] = ntout;
            outsize[1] = npos * nf;
        } else
        {
            //outsize = (mwSize *) mxMalloc(2 * sizeof(*outsize)); /* NGL 2019-08-29 */
            outsize[0] = npos;
            outsize[1] = nf;
        }
        mxSetDimensions(plhs[0], outsize, 2); /* Set to 2D array */
        mxSetDimensions(plhs[1], outsize, 2); /* Set to 2D array */
        mxSetDimensions(plhs[2], outsize, 2); /* Set to 2D array */
    }

    /* ====== Free up allocated memory, if necessary. ===== */
    //mxFree(outsize); /* NGL 08/29/2019 */
    
    if (!mxIsComplex(prhs[0]))
    {
        free(b1i);	/* We had to allocate this before */
    }

    if (mxGetM(prhs[2]) * mxGetN(prhs[2]) == 1)
    {
        free(tp);	/* We had to allocate this */
    }

    if (ti != NULL)
    {
        free(ti);
    }

    if (dyaflag == 1)
    {
        free(dy);
    }
    if (dzaflag == 1)
    {
        free(dz);
    }
    if (gyaflag == 1)
    {
        free(gy);
    }
    if (gzaflag == 1)
    {
        free(gz);
    }
}