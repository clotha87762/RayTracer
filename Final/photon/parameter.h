#ifndef __parameter_h__
#define __parameter_h__


#define MAX_PHOTON 1000000
#define PI 3.141592

#define epsilon 0.000001
#define TOLERANCE_DEPTH 0.001

#define KERNEL_SIZE 6

const bool culling = true;

const int PHOTON_COUNT = 20000;
const int CAUSTIC_COUNT = 20000;

const float QUERY_MAX_DIS = 1.0f;
const int PHOTON_QUERY_COUNT = 100;
const int GLOBAL_ILLUMINATION_SAMPLE = 50;

const int MIN_PHOTON_TO_ESTIMATE = 8;

const float totalLightArea = 2.0f;

const int PHOTON_TRACE_MAX_DEPTH = 6;
const int IRRADIANCE_MAX_DEPTH = 3;

enum QueryMode {NAIVE, BALANCE_TREE};

const QueryMode QUERY_MODE = BALANCE_TREE;

#endif
