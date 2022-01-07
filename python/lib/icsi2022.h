// define for basic functions
void ellips_func(double *, double *, int , int ,double *, int ,double *); /* Ellipsoidal */
void bent_cigar_func(double *, double *, int , int ,double *, int ,double *); /* Discus */
void discus_func(double *, double *, int , int ,double *, int ,double *);  /* Bent_Cigar */
void rosenbrock_func(double *, double *, int , int ,double *, int ,double *); /* Rosenbrock's */
void ackley_func(double *, double *, int , int ,double *, int ,double *); /* Ackley's */
void rastrigin_func(double *, double *, int , int ,double *, int ,double *); /* Rastrigin's  */
void griewank_func(double *, double *, int , int ,double *, int ,double *); /* Griewank's  */
void schwefel_func(double *, double *, int , int ,double *, int ,double *); /* Schwefel's */
void grie_rosen_func(double *, double *, int , int ,double *, int ,double *); /* Griewank-Rosenbrock  */
void escaffer6_func(double *, double *, int , int ,double *, int ,double *); /* Expanded Scaffers F6  */
void happycat_func(double *, double *, int , int ,double *, int ,double *); /* HappyCat */
void weierstrass_func(double *, double *, int , int ,double *, int ,double *); /* Weierstrass's  */
void katsuura_func(double *, double *, int , int ,double *, int ,double *); /* Katsuura */


void alple_func(double *, double *, int , int ,double *, int ,double *); /* alple 01 */
void dixon_func(double *, double *, int , int ,double *, int ,double *); /* Dixon-Prince's Function */

// not use 


void ex3_func(double *, double *, int , int ,double *, int ,double *); 
void logexp_func(double *, double *, int , int ,double *, int ,double *); 
void invert_cos_wave_func(double *, double *, int , int ,double *, int ,double *); 
void patho_func(double *, double *, int , int ,double *, int ,double *); 
void salomon_func(double *, double *, int , int ,double *, int ,double *); 
void sargan_func(double *, double *, int , int ,double *, int ,double *); 
void wavy_func(double *, double *, int , int ,double *, int ,double *); 



// for composition functions

void cf01 (double *, double *, int ,  int ,double *, int ,double *); /* Composition Function 1 */
void cf02 (double *, double *, int ,  int ,double *, int ,double *); /* Composition Function 2 */
void cf03 (double *, double *, int ,  int ,double *, int ,double *); /* Composition Function 3 */
// void cf04 (double *, double *, int ,  int ,double *, int ,double *); /* Composition Function 4 */


// define for basic operator
void shiftfunc (double*,double*,int,double*);
void rotatefunc (double*,double*,int, double*);
void sr_func (double *, double *, int, double*, double*, double, int, int); /* shift and rotate */
// composition operator
void cf_cal2(double *, double *,double *,double *, int);

// api
void ceval(double *x, int nx, int mx, double *f, int func_num0, char * path);