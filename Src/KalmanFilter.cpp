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
	ssSetNumSFcnParams(S, 0);                   // no parameters
	if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S))
		return;	// Parameter mismatch reported by the Simulink engine

	if (!ssSetNumInputPorts(S, 2)) return;      // 2 input ports
	ssSetInputPortWidth(S, 0, 1);               // input port 0 width = 1
	ssSetInputPortWidth(S, 1, 1);               // input port 1 width = 1

	ssSetInputPortDirectFeedThrough(S, 0, 1);   // input 0 is used
	ssSetInputPortDirectFeedThrough(S, 1, 1);   // input 1 is used

	if (!ssSetNumOutputPorts(S, 2)) return;     // 2 output ports
	ssSetOutputPortWidth(S, 0, 1);              // output 0 width = 1
	ssSetOutputPortWidth(S, 1, 1);              // output 1 width = 1

	ssSetNumSampleTimes(S, 1);                  // 1 block-based sample time

	ssSetNumPWork(S, 2);                        // reserve 2 pointers for our objects
	
	ssSetSimStateCompliance(S, USE_CUSTOM_SIM_STATE);

	ssSetOptions(S, SS_OPTION_CALL_TERMINATE_ON_EXIT | SS_OPTION_WORKS_WITH_CODE_REUSE);
}
#define MDL_INITIALIZE_SAMPLE_TIMES
static void mdlInitializeSampleTimes(SimStruct *S)
{
	ssSetSampleTime(S, 0, 0.1);                 // frequency = 10Hz
	ssSetOffsetTime(S, 0, 0.0);                 // no offset
}
#define MDL_START
static void mdlStart(SimStruct *S)
{
	ISysMat *psm = new SysMat;
	Filter *pf = new Filter(psm);
	ssGetPWork(S)[0] = (void*)pf;
	ssGetPWork(S)[1] = (void*)psm;
}
//#define MDL_UPDATE
//static void mdlUpdate(SimStruct *S, int_T tid)    // mdlUpdate() shouldn't be here as 
//{                                                 // we have direct feedthrough and 
//	Filter *pf = (Filter*)ssGetPWork(S)[0];         // no discrete states
//	double time = (double)ssGetT(S);                // We don't use state vectors either
//	pf->ProjectAhead();
//	pf->set_Time(time);
//}
#define MDL_OUTPUTS
static void mdlOutputs(SimStruct *S, int_T task_id)
{
	Filter *pf = (Filter*)ssGetPWork(S)[0];

	// Step 1:
	pf->ComputeGain();

	// Retrieving inputs:
	InputRealPtrsType input = ssGetInputPortRealSignalPtrs(S, 0);
	double inputs[2] = 
	{ (double)(*input[0]) , (double)(*input[1]) };  // let's hope they are convertible
	Matrix y(1,1);
	y.at(0, 0) = inputs[0];                         // measurement is first input
	Matrix u(1, 1);
	u.at(0, 0) = inputs[1];                         // reference signal is second input

	// Step 2:
	pf->UpdateEstimate(y);

	// Step 3:
	pf->ComputeCovariance();

	// Step 4:
	pf->ProjectAhead(u);

	// Updating time:
	double time = (double)ssGetT(S);
	pf->set_Time(time);

	// Setting outputs:
	real_T *output = ssGetOutputPortRealSignal(S, 0);
	output[0] = (real_T)(pf->get_State(2));         // first output is psi
	output[1] = (real_T)(pf->get_State(4));         // second output is b
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