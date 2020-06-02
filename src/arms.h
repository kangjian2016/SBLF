/* header file for arms function */

int arms_simple (int ninit, float *xl, float *xr,
	         float (*myfunc)(float x, void *mydata), void *mydata,
                 int dometrop, float *xprev, float *xsamp);

int arms (float *xinit, int ninit, float *xl, float *xr,
	 float (*myfunc)(float x, void *mydata), void *mydata,
         float *convex, int npoint, int dometrop, float *xprev, float *xsamp,
         int nsamp, float *qcent, float *xcent, int ncent,
         int *neval);

float expshift(float y, float y0);

#define YCEIL 50.                /* maximum y avoiding overflow in exp(y) */

