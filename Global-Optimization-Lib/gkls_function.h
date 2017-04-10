
#ifndef GKLS_FUNCTION_H
#define GKLS_FUNCTION_H

#include <memory>

#include "gkls/gkls.h"
#include "gkls/rnd_gen.h"

#include "map/map.h"
#include "math_function.h"

class GklsFunction : public MathFunction
{
public:
	typedef std::shared_ptr<GklsFunction> GklsFunctionPtr;

	enum FunctionType
	{
		ND, D, D2
	};

	/*---------------- Variables accessible by the user -------------------- */
	static unsigned int GKLS_dim;    /* dimension of the problem,        */
	/* 2<=GKLS_dim<NUM_RND (see random) */
	static unsigned int GKLS_num_minima; /* number of local minima, >=2  */

	static double GKLS_global_dist;  /* distance from the paraboloid minimizer  */
	/* to the global minimizer                 */
	static double GKLS_global_radius;/* radius of the global minimizer          */
	/* attraction region                       */
	static double GKLS_global_value; /* global minimum value,                   */

	static int GKLS_parameters_check(void);/* test the validity of the input parameters*/

private:
	FunctionType func_type = ND;
	unsigned int sequential_number;
	static const int KEY = 1, PRECISION = 10;
	static void linear_transform(double point[]);

	/*---------------- Variables accessible by the user -------------------- */
	static double *GKLS_domain_left; /* left boundary vector of D  */
	/* D=[GKLS_domain_left; GKLS_domain_right] */
	static double *GKLS_domain_right;/* right boundary vector of D */

	/* GKLS_global_value < GKLS_PARABOLOID_MIN */
	T_GKLS_Minima GKLS_minima;
	/* see the structures type description     */
	T_GKLS_GlobalMinima GKLS_glob;


	/*--------------------------- Global variables ----------------------*/
	int isArgSet = 0; /* isArgSet == 1 if all necessary parameters are set */

	double delta; /* parameter using in D2-type function generation;     */
	/* it is chosen randomly from the                      */
	/* open interval (0,GKLS_DELTA_MAX_VALUE)              */
	unsigned long rnd_counter; /* index of random array elements */


	void GKLS_free(void);        /* deallocate memory needed for the generator  */



	int GKLS_arg_generate(unsigned int); /* test function generator */

	double GKLS_ND_func(double *);  /* evaluation of an ND-typed test function  */

	double GKLS_D_func(double *);   /* evaluation of a D-typed test function    */

	double GKLS_D2_func(double *);  /* evaluation of a D2-type test function    */

	/*------------------ Auxiliary functions prototypes -----------------*/
	int GKLS_alloc(void);
	int GKLS_coincidence_check(void);
	int GKLS_set_basins(void);
	int GKLS_initialize_rnd(unsigned int, unsigned int, int);

public:
	GklsFunction(unsigned int _sequential_number = 1);
	GklsFunction(FunctionType, unsigned int _sequential_number = 1);
	~GklsFunction();
	GklsFunction(const GklsFunction&);

	void setFunctionType(FunctionType);
	static void mapScalar_To_nDSpace(double x, double point[]);
	void getFunctionMinimum(double min_pt[]);

	double getValue(double arg) override;

	static double GKLS_norm(double *, double *);

	/*------------------------User function prototypes -------------------------*/

	static int GKLS_domain_alloc(void); /* allocate boundary vectors   */

	static int  GKLS_set_default(void); /* set default values of the input parameters  */
	/* and allocate the boundary vectors if necessary */

	static void GKLS_domain_free(void); /* deallocate boundary vectors */
};

#endif
