//
// Incorporation of filter for a particular system into
// a Simulink S-function using Simulink/Matlab API.
// This file is to be compiled as a mex-file
// 
#include"filter.h"
#include"udf.h"
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

	ssSetOptions(S, SS_OPTION_WORKS_WITH_CODE_REUSE |
                 SS_OPTION_RUNTIME_EXCEPTION_FREE_CODE);
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

    Sysmat *psm = new Sysmat;
    Filter<Dim>* pf = new Filter<Dim>(psm);
	ssGetPWork(S)[0] = (void*)pf;
	ssGetPWork(S)[1] = (void*)psm;
}
//#define MDL_UPDATE
//static void mdlUpdate(SimStruct *S, int_T tid)    // mdlUpdate() shouldn't be here as 
//{                                                 // we have direct feedthrough and 
//    Filter *pf = (Filter*)ssGetPWork(S)[0];       // no discrete states
//    double time = (double)ssGetT(S);              // We don't use state vectors either
//    pf->ProjectAhead();
//    pf->set_Time(time);
//}
#define MDL_OUTPUTS
static void mdlOutputs(SimStruct *S, int_T task_id)
{
	Filter<Dim> *pf = (Filter<Dim>*)ssGetPWork(S)[0];

	// Step 1:
	pf->ComputeGain();

	// Retrieving inputs:
	InputRealPtrsType input = ssGetInputPortRealSignalPtrs(S, 0);
	double inputs[2] = 
	{ (double)(*input[0]) , (double)(*input[1]) };  // let's hope they are convertible
	Filter<Dim>::VecY y;
	y(0, 0) = inputs[0];                            // measurement is first input
    Filter<Dim>::VecU u;
	u(0, 0) = inputs[1];                            // reference signal is second input

	// Step 2:
	pf->UpdateEstimate(y, u);

	// Step 3:
	pf->ComputeCovariance();

	// Step 4:
	pf->ProjectAhead(u);

	// Updating time:
	double time = (double)(ssGetT(S));
	pf->set_Time(time);

	// Setting outputs:
	real_T *output = ssGetOutputPortRealSignal(S, 0);
	output[0] = (real_T)(pf->get_State(2));         // first output is psi
	output[1] = (real_T)(pf->get_State(4));         // second output is b
}
static void mdlTerminate(SimStruct *S)
{
	Filter<Dim> *pf = (Filter<Dim>*)ssGetPWork(S)[0];
	Sysmat *psm = (Sysmat*)ssGetPWork(S)[1];
	delete pf; delete psm;
}
#ifdef MATLAB_MEX_FILE
 #include"simulink.c"
#else
 #include"cg_fun.h"
#endif