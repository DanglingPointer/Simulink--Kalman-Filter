# Kalman filter
Generic implementation of [Kalman filter](https://en.wikipedia.org/wiki/Kalman_filter) written in C++ and integrated with Simulink Api in order to be used as a Simulink S-function for a particular linearized system. The compiled mex-file as well as Matlab and Simulink files are to be found in the folder `Simulink&Matlab/`. 

The C++ code consist of the following files:
* `matrix.h`: implementation of matrices which may be reused for other purposes
* `filter.h`: implementation of the Kalman filter, may be used for any linear system. User has to parametrize the `Mattypes` type wrapper and implement a realization of the `ISysmat` interface
* `udf.h`: user-defined functionalities described above, implemented for the particuar system
* `KalmanFilter.cpp`: integration with Simulink/Matlab API, to be compiled as a mex-file in order to use in a Simulink S-function
* `ConsoleDebug.cpp`: for testing the filter in Visual Studio, should be ignored otherwise. When used, the content of the `KalmanFilter.cpp` should be commented out

The Simulink model with the filter was used in a project as a part of the course in TTK4115 Linear System Theory at NTNU.
