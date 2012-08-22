/* 
 * Compute global cost matrix - companion
 * to the dtw R package
 * (c) Toni Giorgino  2007
 * Distributed under GPL-2 with NO WARRANTY.
 *  
 */

int argmin(double *list, int n);

void computeCM(			/* IN */
	       int *s,		/* mtrx dimensions */
	       int *wm,		/* windowing matrix */
	       double *lm,	/* local cost mtrx */
	       int *nstepsp,	/* no of steps in stepPattern */
	       double *dir,	/* stepPattern description */
				/* OUT */
	       double *cm,	/* cost matrix */
	       int *sm		/* direction mtrx */

				);


#ifdef TEST_UNIT
void tm_print(int *s, double *mm, double *r);
#endif
