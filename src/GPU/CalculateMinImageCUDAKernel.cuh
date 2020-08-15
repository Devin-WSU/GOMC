/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#pragma once
#ifdef GOMC_CUDA

#include <cuda.h>
#include <cuda_runtime.h>
#include "ConstantDefinitionsCUDAKernel.cuh"

__device__ inline double3 Difference(double * x, double * y, double * z, uint i, uint j){
  return make_double3(x[i] - x[j], y[i] - y[j], z[i] - z[j]);
}

__device__ inline void TransformSlantGPU(double &tx, double &ty, double &tz,
    double x, double y, double z,
    double *gpu_cell_x,
    double *gpu_cell_y,
    double *gpu_cell_z)
{
  tx = x * gpu_cell_x[0] + y * gpu_cell_x[1] + z * gpu_cell_x[2];
  ty = x * gpu_cell_y[0] + y * gpu_cell_y[1] + z * gpu_cell_y[2];
  tz = x * gpu_cell_z[0] + y * gpu_cell_z[1] + z * gpu_cell_z[2];
}

__device__ inline void TransformUnSlantGPU(double &tx, double &ty, double &tz,
    double x, double y, double z,
    double *gpu_Invcell_x,
    double *gpu_Invcell_y,
    double *gpu_Invcell_z)
{
  tx = x * gpu_Invcell_x[0] + y * gpu_Invcell_x[1] + z * gpu_Invcell_x[2];
  ty = x * gpu_Invcell_y[0] + y * gpu_Invcell_y[1] + z * gpu_Invcell_y[2];
  tz = x * gpu_Invcell_z[0] + y * gpu_Invcell_z[1] + z * gpu_Invcell_z[2];
}

__device__ inline double MinImageSignedGPU(double raw, double ax, double halfAx)
{
  if (raw > halfAx)
    raw -= ax;
  else if (raw < -halfAx)
    raw += ax;
  return raw;
}

__device__ inline double3 MinImageGPU(double3 rawVec, double3 axis, double3 halfAx){
  rawVec.x = MinImageSignedGPU(rawVec.x, axis.x, halfAx.x);
  rawVec.y = MinImageSignedGPU(rawVec.y, axis.y, halfAx.y);
  rawVec.z = MinImageSignedGPU(rawVec.z, axis.z, halfAx.z);
  return rawVec;
}

// Call by calculate energy whether it is in rCut
  __device__ inline bool InRcutGPU(double &distSq, 
                                   double * x, double * y, double * z,
                                   uint i, uint j,
                                   double3 axis, double3 halfAx, 
                                   double gpu_rCut, int gpu_nonOrth,
                                   double *gpu_cell_x, double *gpu_cell_y,
                                   double *gpu_cell_z, double *gpu_Invcell_x,
                                   double *gpu_Invcell_y, double *gpu_Invcell_z)
{
  distSq = 0;
  double3 t, dist;
  dist = Difference(x, y, z, i, j);
  // Do a binary print here of dist


  printf("Before MinImage dist (%d, %d) = \n\tx - %f, in binary : %lu", i,  j, dist.x, (union { double d; uint64_t u; }) {dist.x} .u);
  printf("\n\ty - %f, in binary : %lu", dist.y, (union { double d; uint64_t u; }) {dist.y} .u);
  printf("\n\tz - %f, in binary : %lu", dist.z, (union { double d; uint64_t u; }) {dist.z} .u);
  printf("\n");



  if(gpu_nonOrth) {
    TransformUnSlantGPU(t.x, t.y, t.z, dist.x, dist.y, dist.z, gpu_Invcell_x,
                        gpu_Invcell_y, gpu_Invcell_z);
    t = MinImageGPU(t, axis, halfAx);
    TransformSlantGPU(dist.x, dist.y, dist.z, t.x, t.y, t.z, gpu_cell_x, gpu_cell_y,
                      gpu_cell_z);
  } else {
    dist = MinImageGPU(dist, axis, halfAx);
    distSq = dist.x * dist.x + dist.y * dist.y + dist.z * dist.z;
  }

  printf("After MinImage dist (%d, %d) = \n\tx - %f, in binary : %lu", i,  j, dist.x, (union { double d; uint64_t u; }) {dist.x} .u);
  printf("\n\ty - %f, in binary : %lu", dist.y, (union { double d; uint64_t u; }) {dist.y} .u);
  printf("\n\tz - %f, in binary : %lu", dist.z, (union { double d; uint64_t u; }) {dist.z} .u);
  printf("\n");

  printf("distSq : %f, in binary : %lu\n", distSq, (union { double d; uint64_t u; }) {distSq} .u);
  printf("rCutSq : %f, in binary : %lu\n", gpu_rCut * gpu_rCut, (union { double d; uint64_t u; }) {gpu_rCut * gpu_rCut} .u);


  return ((gpu_rCut * gpu_rCut) > distSq);
}


// Call by force calculate to return the distance and virial component
__device__ inline bool InRcutGPU(double &distSq, double3 & dist, 
                                 double * x, double * y, double * z,
                                 uint i, uint j,
                                 double3 axis, double3 halfAx, 
                                 double gpu_rCut, int gpu_nonOrth,
                                 double *gpu_cell_x, double *gpu_cell_y,
                                 double *gpu_cell_z, double *gpu_Invcell_x,
                                 double *gpu_Invcell_y, double *gpu_Invcell_z)
{
  distSq = 0;
  double3 t;
  dist = Difference(x, y, z, i, j);
  // Do a binary print here of dist


  printf("Before MinImage dist (%d, %d) = \n\tx - %f, in binary : %lu", i,  j, dist.x, (union { double d; uint64_t u; }) {dist.x} .u);
  printf("\n\ty - %f, in binary : %lu", dist.y, (union { double d; uint64_t u; }) {dist.y} .u);
  printf("\n\tz - %f, in binary : %lu", dist.z, (union { double d; uint64_t u; }) {dist.z} .u);
  printf("\n");



  if(gpu_nonOrth) {
    printf("nonorth\n");
    TransformUnSlantGPU(t.x, t.y, t.z, dist.x, dist.y, dist.z, gpu_Invcell_x,
                        gpu_Invcell_y, gpu_Invcell_z);
    t = MinImageGPU(t, axis, halfAx);
    TransformSlantGPU(dist.x, dist.y, dist.z, t.x, t.y, t.z, gpu_cell_x, gpu_cell_y,
                      gpu_cell_z);
  } else {
    printf("orth\n");
    dist = MinImageGPU(dist, axis, halfAx);
    //distSq = dist.x * dist.x + dist.y * dist.y + dist.z * dist.z;
    distSq = fma(dist.x, dist.x, fma(dist.y, dist.y, dist.z * dist.z));
  }

  printf("After MinImage dist (%d, %d) = \n\tx - %f, in binary : %lu", i,  j, dist.x, (union { double d; uint64_t u; }) {dist.x} .u);
  printf("\n\ty - %f, in binary : %lu", dist.y, (union { double d; uint64_t u; }) {dist.y} .u);
  printf("\n\tz - %f, in binary : %lu", dist.z, (union { double d; uint64_t u; }) {dist.z} .u);
  printf("\n");

  printf("distSq : %f, in binary : %lu\n", distSq, (union { double d; uint64_t u; }) {distSq} .u);
  printf("rCutSq : %f, in binary : %lu\n", gpu_rCut * gpu_rCut, (union { double d; uint64_t u; }) {gpu_rCut * gpu_rCut} .u);

  if (gpu_rCut * gpu_rCut > distSq){
    printf("InRCut : true\n");
  } else {
    printf("InRCut : false\n");
  }

  return ((gpu_rCut * gpu_rCut) > distSq);
}

__device__ inline int FlatIndexGPU(int i, int j, int gpu_count)
{
  return i + j * gpu_count;
}

__device__ inline double DotProductGPU(double kx, double ky, double kz,
                                       double x, double y, double z)
{
  return (kx * x + ky * y + kz * z);
}

__device__ inline double DeviceGetLambdaVDW(int molA, int kindA, int molB,
    int kindB, int box,
    bool *gpu_isFraction,
    int *gpu_molIndex,
    int *gpu_kindIndex,
    double *gpu_lambdaVDW)
{
  double lambda = 1.0;
  if(gpu_isFraction[box]) {
    if((gpu_molIndex[box] == molA) && (gpu_kindIndex[box] == kindA)) {
      lambda *= gpu_lambdaVDW[box];
    }
    if((gpu_molIndex[box] == molB) && (gpu_kindIndex[box] == kindB)) {
      lambda *= gpu_lambdaVDW[box];
    }
  }
  return lambda;
}

__device__ inline double DeviceGetLambdaCoulomb(int molA, int kindA, int molB,
    int kindB, int box,
    bool *gpu_isFraction,
    int *gpu_molIndex,
    int *gpu_kindIndex,
    double *gpu_lambdaCoulomb)
{
  double lambda = 1.0;
  if(gpu_isFraction[box]) {
    if((gpu_molIndex[box] == molA) && (gpu_kindIndex[box] == kindA)) {
      lambda *= gpu_lambdaCoulomb[box];
    }
    if((gpu_molIndex[box] == molB) && (gpu_kindIndex[box] == kindB)) {
      lambda *= gpu_lambdaCoulomb[box];
    }
  }
  return lambda;
}

// Add atomic operations for GPUs that do not support it
// atomicAdd and atomicSub only support double for Compute Capability >= 6.0
#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
static __inline__ __device__ double atomicAdd(double *address, double val)
{
  unsigned long long int* address_as_ull = (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  if (val == 0.0)
    return __longlong_as_double(old);
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif

#endif /*GOMC_CUDA*/
