/*
 * Incorporation of filter into a Simulink S-function
 * using Simulink API.
 * This file is to be compiled to a mex-file
 */
#include"kffilter.h"
#include"kfudf.h"
#define S_FUNCTION_NAME KalmanFilter
#define S_FUNCTION_LEVEL 2
#include "simstruc.h"
using namespace mvkf;

#define MDL_INITIAL_SIZES
static void mdlInitializeSizes(SimStruct *S)
{
	ssSetNumSFcnParams(S, 0);					// no parameters
	if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S))
		return;	// Parameter mismatch reported by the Simulink engine

	if (!ssSetNumInputPorts(S, 2)) return;		// 2 input ports
	ssSetInputPortWidth(S, 0, 1);				// input port 0 width = 1
	ssSetInputPortWidth(S, 1, 1);				// input port 1 width = 1

	ssSetInputPortDirectFeedThrough(S, 0, 1);	// input 0 is used
	ssSetInputPortDirectFeedThrough(S, 1, 1);	// input 1 is used

	if (!ssSetNumOutputPorts(S, 2)) return;		// 2 output ports
	ssSetOutputPortWidth(S, 0, 1);				// output 0 width = 1
	ssSetOutputPortWidth(S, 1, 1);				// output 1 width = 1

	ssSetNumSampleTimes(S, 1);					// 1 block-based sample time

	ssSetNumPWork(S, 2);						// reserve 2 pointers for our objects
	
	ssSetSimStateCompliance(S, USE_CUSTOM_SIM_STATE);

	ssSetOptions(S, SS_OPTION_CALL_TERMINATE_ON_EXIT | SS_OPTION_WORKS_WITH_CODE_REUSE);
}
#define MDL_INITIALIZE_SAMPLE_TIMES
static void mdlInitializeSampleTimes(SimStruct *S)
{
	ssSetSampleTime(S, 0, 0.1);					// frequency = 10Hz
	ssSetOffsetTime(S, 0, 0.0);					// no offset
}
#define MDL_START
static void mdlStart(SimStruct *S)
{
	ISysMat *psm = new SysMat;
	Filter *pf = new Filter(psm);
	pf->InitStates();

	ssGetPWork(S)[0] = (void*)pf;
	ssGetPWork(S)[1] = (void*)psm;
}
#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid)
{
	Filter *pf = (Filter*)ssGetPWork(S)[0];
	double time = (double)ssGetT(S);
	pf->ProjectAhead();	// må hente input!
	pf->set_Time(time);
}
#define MDL_OUTPUTS
static void mdlOutputs(SimStruct *S, int_T tid)
{
	Filter *pf = (Filter*)ssGetPWork(S)[0];
	pf->ComputeGain();

	// Getting inputs
	InputRealPtrsType u = ssGetInputPortRealSignalPtrs(S, 0);

	double inputs[2] = {
		(double)(*u[0]) , (double)(*u[1]) };	// let's hope they are convertible

	IMatrix *pinp = new Matrix<2, 1>; // wrong dimensions! should be <1,1>
	pinp->at(0, 0) = inputs[0];
	pinp->at(1, 0) = inputs[1];

	pf->UpdateEstimate(pinp);
	delete pinp;

	pf->ComputeCovariance();

	// Setting outputs
	real_T *y = ssGetOutputPortRealSignal(S, 0);
	y[0] = (real_T)(pf->get_State(2));			// let's hope they are convertible
	y[1] = (real_T)(pf->get_State(4));			// 2 = psi, 4 = b
}
static void mdlTerminate(SimStruct *S)
{
	Filter *pf = (Filter*)ssGetPWork(S)[0];
	ISysMat *psm = (ISysMat*)ssGetPWork(S)[1];
	delete pf; delete psm;
}
#ifdef MATLAB_MEX_FILE
 #include"simulink.c"
#else
 #include"cg_fun.h"
#endif