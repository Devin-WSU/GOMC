/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifdef GOMC_CUDA
#include <cuda.h>
#include "cub/cub.cuh"
#include <stdio.h>
#include "ConstantDefinitionsCUDAKernel.cuh"
#include "CalculateMinImageCUDAKernel.cuh"
#include "CalculateForceCUDAKernel.cuh"
#include "CalculateEnergyCUDAKernel.cuh"
#include "CUDAMemoryManager.cuh"
#include <vector>
#define NUMBER_OF_NEIGHBOR_CELL 27

using namespace cub;

struct Lock;

void CallBoxInterGPU(VariablesCUDA *vars,
                     std::vector<int> cellVector,
                     std::vector<int> cellStartIndex,
                     std::vector<std::vector<int> > neighborList,
                     XYZArray const &coords,
                     BoxDimensions const &boxAxes,
                     bool electrostatic,
                     std::vector<double> particleCharge,
                     std::vector<int> particleKind,
                     std::vector<int> particleMol,
                     double &REn,
                     double &LJEn,
                     bool sc_coul,
                     double sc_sigma_6,
                     double sc_alpha,
                     uint sc_power,
                     uint const box,
                     int * hostEnergyVectorLJKeys,
                     int * hostEnergyVectorREnKeys,
                     double * hostEnergyVectorLJValues,
                     double * hostEnergyVectorREnValues,
                     uint * numberOfInters
                    )
{
  int atomNumber = coords.Count();
  int cellVectorCount = cellVector.size();
  int neighborListCount = neighborList.size() * NUMBER_OF_NEIGHBOR_CELL;
  int numberOfCells = neighborList.size();
  int *gpu_particleKind, *gpu_particleMol;
  int *gpu_neighborList, *gpu_cellStartIndex;
  int blocksPerGrid, threadsPerBlock;
  int energyVectorLen = 0;
  double *gpu_particleCharge;
  double *gpu_REn, *gpu_LJEn;
  double *gpu_final_REn, *gpu_final_LJEn;

  int *gpu_REnForSortAndCPUReductionKeys, *gpu_LJEnForSortAndCPUReductionKeys;
  double *gpu_REnForSortAndCPUReductionValues, *gpu_LJEnForSortAndCPUReductionValues, *gpu_REnSortedValues, *gpu_LJEnSortedValues;



  // Run the kernel
  threadsPerBlock = 256;
  blocksPerGrid = (int)(numberOfCells * NUMBER_OF_NEIGHBOR_CELL);
  energyVectorLen = numberOfCells * NUMBER_OF_NEIGHBOR_CELL * threadsPerBlock;

  // Convert neighbor list to 1D array
  std::vector<int> neighborlist1D(neighborListCount);
  for(int i = 0; i < neighborList.size(); i++) {
    for(int j = 0; j < NUMBER_OF_NEIGHBOR_CELL; j++) {
      neighborlist1D[i * NUMBER_OF_NEIGHBOR_CELL + j] = neighborList[i][j];
    }
  }


  CUMALLOC((void**) &gpu_REnForSortAndCPUReductionKeys, energyVectorLen * sizeof(int));
  CUMALLOC((void**) &gpu_LJEnForSortAndCPUReductionKeys, energyVectorLen * sizeof(int));
  CUMALLOC((void**) &gpu_REnForSortAndCPUReductionValues, energyVectorLen * sizeof(double));
  CUMALLOC((void**) &gpu_LJEnForSortAndCPUReductionValues, energyVectorLen * sizeof(double));

  CUMALLOC((void**) &gpu_REnSortedValues, energyVectorLen * sizeof(double));
  CUMALLOC((void**) &gpu_LJEnSortedValues, energyVectorLen * sizeof(double));

  CUMALLOC((void**) &gpu_neighborList, neighborListCount * sizeof(int));
  CUMALLOC((void**) &gpu_cellStartIndex, cellStartIndex.size() * sizeof(int));
  CUMALLOC((void**) &gpu_particleCharge, particleCharge.size() * sizeof(double));
  CUMALLOC((void**) &gpu_particleKind, particleKind.size() * sizeof(int));
  CUMALLOC((void**) &gpu_particleMol, particleMol.size() * sizeof(int));
  CUMALLOC((void**) &gpu_REn, energyVectorLen * sizeof(double));
  CUMALLOC((void**) &gpu_LJEn, energyVectorLen * sizeof(double));
  CUMALLOC((void**) &gpu_final_REn, sizeof(double));
  CUMALLOC((void**) &gpu_final_LJEn, sizeof(double));

  // Copy necessary data to GPU
  cudaMemcpy(gpu_neighborList, &neighborlist1D[0], neighborListCount * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_cellStartIndex, &cellStartIndex[0], cellStartIndex.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_cellVector, &cellVector[0], atomNumber * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleCharge, &particleCharge[0], particleCharge.size() * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleKind, &particleKind[0], particleKind.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleMol, &particleMol[0], particleMol.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(double), cudaMemcpyHostToDevice);

  uint * numIntersPerCell = new uint[blocksPerGrid];
  uint * gpu_numIntersPerCell;
  CUMALLOC((void**) &gpu_numIntersPerCell, blocksPerGrid * sizeof(uint));
  
  CalcNumberOfInteractions <<< blocksPerGrid, threadsPerBlock>>>(gpu_cellStartIndex,
    vars->gpu_cellVector,
    gpu_neighborList,
    numberOfCells,
    cellVectorCount,
    vars->gpu_x,
    vars->gpu_y,
    vars->gpu_z,
    boxAxes.GetAxis(box).x,
    boxAxes.GetAxis(box).y,
    boxAxes.GetAxis(box).z,
    electrostatic,
    gpu_particleCharge,
    gpu_particleKind,
    gpu_particleMol,
    gpu_REn,
    gpu_LJEn,
    vars->gpu_sigmaSq,
    vars->gpu_epsilon_Cn,
    vars->gpu_n,
    vars->gpu_VDW_Kind,
    vars->gpu_isMartini,
    vars->gpu_count,
    vars->gpu_rCut,
    vars->gpu_rCutCoulomb,
    vars->gpu_rCutLow,
    vars->gpu_rOn,
    vars->gpu_alpha,
    vars->gpu_ewald,
    vars->gpu_diElectric_1,
    vars->gpu_nonOrth,
    vars->gpu_cell_x[box],
    vars->gpu_cell_y[box],
    vars->gpu_cell_z[box],
    vars->gpu_Invcell_x[box],
    vars->gpu_Invcell_y[box],
    vars->gpu_Invcell_z[box],
    sc_coul,
    sc_sigma_6,
    sc_alpha,
    sc_power,
    vars->gpu_rMin,
    vars->gpu_rMaxSq,
    vars->gpu_expConst,
    vars->gpu_molIndex,
    vars->gpu_kindIndex,
    vars->gpu_lambdaVDW,
    vars->gpu_lambdaCoulomb,
    vars->gpu_isFraction,
    box,
    gpu_numIntersPerCell);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  CubDebugExit(cudaMemcpy(numIntersPerCell, gpu_numIntersPerCell, blocksPerGrid * sizeof(uint),
  cudaMemcpyDeviceToHost));


  uint totalNumInters = 0;
  for (int i = 0; i < blocksPerGrid; i++){
    totalNumInters+=numIntersPerCell[i];
  }  
  CUMALLOC((void**) &vars->numberOfInters, sizeof(uint));
  // Set this for the flattened force calc
  CubDebugExit(cudaMemcpy(vars->numberOfInters, &totalNumInters, sizeof(uint), cudaMemcpyHostToDevice));

/*
  for (int i = 0; i < blocksPerGrid; i++){
    std::cout << "numIntersPerCell[" << i << "] : " << numIntersPerCell[i] << std::endl;
  }
*/
  uint tmp = 0;
  for (int i = 1; i <= blocksPerGrid; i++){
    for (int j = i-1; j < i; j++){
      numIntersPerCell[j] = tmp;
      tmp += numIntersPerCell[i];
    }
  }
/*
  std::cout << "number of interactions : " << totalNumInters << std::endl;
  for (int i = 0; i < blocksPerGrid; i++){
    std::cout << "numIntersPerCell[" << i << "] : " << numIntersPerCell[i] << std::endl;
  }
*/
  double * flatgpu_REn;
  double * flatgpu_LJEn;
  CUMALLOC((void**) &flatgpu_REn, totalNumInters * sizeof(double));
  CUMALLOC((void**) &flatgpu_LJEn, totalNumInters * sizeof(double));

  hostEnergyVectorLJKeys = (int*) malloc( sizeof(int) * totalNumInters);  
  hostEnergyVectorREnKeys = (int*) malloc( sizeof(int) * totalNumInters);
  
  int flatIndexREn = -1;
  int flatIndexLJEn = -1;

  int * gpu_flatIndexREn;
  int * gpu_flatIndexLJEn;
  CUMALLOC((void**) &gpu_flatIndexREn, sizeof(int));
  CUMALLOC((void**) &gpu_flatIndexLJEn, sizeof(int));

  CubDebugExit(cudaMemcpy(gpu_flatIndexREn, &flatIndexREn, sizeof(int), cudaMemcpyHostToDevice));
  CubDebugExit(cudaMemcpy(gpu_flatIndexLJEn, &flatIndexLJEn, sizeof(int), cudaMemcpyHostToDevice));

  BoxInterGPUFlattened <<< blocksPerGrid, threadsPerBlock>>>(gpu_cellStartIndex,
    vars->gpu_cellVector,
    gpu_neighborList,
    numberOfCells,
    cellVectorCount,
    vars->gpu_x,
    vars->gpu_y,
    vars->gpu_z,
    boxAxes.GetAxis(box).x,
    boxAxes.GetAxis(box).y,
    boxAxes.GetAxis(box).z,
    electrostatic,
    gpu_particleCharge,
    gpu_particleKind,
    gpu_particleMol,
    gpu_REn,
    gpu_LJEn,
    vars->gpu_sigmaSq,
    vars->gpu_epsilon_Cn,
    vars->gpu_n,
    vars->gpu_VDW_Kind,
    vars->gpu_isMartini,
    vars->gpu_count,
    vars->gpu_rCut,
    vars->gpu_rCutCoulomb,
    vars->gpu_rCutLow,
    vars->gpu_rOn,
    vars->gpu_alpha,
    vars->gpu_ewald,
    vars->gpu_diElectric_1,
    vars->gpu_nonOrth,
    vars->gpu_cell_x[box],
    vars->gpu_cell_y[box],
    vars->gpu_cell_z[box],
    vars->gpu_Invcell_x[box],
    vars->gpu_Invcell_y[box],
    vars->gpu_Invcell_z[box],
    sc_coul,
    sc_sigma_6,
    sc_alpha,
    sc_power,
    vars->gpu_rMin,
    vars->gpu_rMaxSq,
    vars->gpu_expConst,
    vars->gpu_molIndex,
    vars->gpu_kindIndex,
    vars->gpu_lambdaVDW,
    vars->gpu_lambdaCoulomb,
    vars->gpu_isFraction,
    box,
    gpu_numIntersPerCell,
    flatgpu_REn,
    flatgpu_LJEn,
    gpu_flatIndexREn,
    gpu_flatIndexLJEn);
cudaDeviceSynchronize();
checkLastErrorCUDA(__FILE__, __LINE__);



/*
  BoxInterGPU <<< blocksPerGrid, threadsPerBlock>>>(gpu_cellStartIndex,
      vars->gpu_cellVector,
      gpu_neighborList,
      numberOfCells,
      cellVectorCount,
      vars->gpu_x,
      vars->gpu_y,
      vars->gpu_z,
      boxAxes.GetAxis(box).x,
      boxAxes.GetAxis(box).y,
      boxAxes.GetAxis(box).z,
      electrostatic,
      gpu_particleCharge,
      gpu_particleKind,
      gpu_particleMol,
      gpu_REn,
      gpu_LJEn,
      vars->gpu_sigmaSq,
      vars->gpu_epsilon_Cn,
      vars->gpu_n,
      vars->gpu_VDW_Kind,
      vars->gpu_isMartini,
      vars->gpu_count,
      vars->gpu_rCut,
      vars->gpu_rCutCoulomb,
      vars->gpu_rCutLow,
      vars->gpu_rOn,
      vars->gpu_alpha,
      vars->gpu_ewald,
      vars->gpu_diElectric_1,
      vars->gpu_nonOrth,
      vars->gpu_cell_x[box],
      vars->gpu_cell_y[box],
      vars->gpu_cell_z[box],
      vars->gpu_Invcell_x[box],
      vars->gpu_Invcell_y[box],
      vars->gpu_Invcell_z[box],
      sc_coul,
      sc_sigma_6,
      sc_alpha,
      sc_power,
      vars->gpu_rMin,
      vars->gpu_rMaxSq,
      vars->gpu_expConst,
      vars->gpu_molIndex,
      vars->gpu_kindIndex,
      vars->gpu_lambdaVDW,
      vars->gpu_lambdaCoulomb,
      vars->gpu_isFraction,
      box);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
*/
  // Copy Keys to device
  /*
  CubDebugExit(cudaMemcpy(gpu_REnForSortAndCPUReductionKeys, hostEnergyVectorREnKeys, energyVectorLen * sizeof(int),
                          cudaMemcpyHostToDevice));
  CubDebugExit(cudaMemcpy(gpu_LJEnForSortAndCPUReductionKeys, hostEnergyVectorLJKeys, energyVectorLen * sizeof(int),
                          cudaMemcpyHostToDevice));
*/
/*
  // Clone energy vector for sorting
  CubDebugExit(cudaMemcpy(gpu_REnForSortAndCPUReductionValues, flatgpu_REn, totalNumInters * sizeof(double),
                          cudaMemcpyDeviceToDevice));
  CubDebugExit(cudaMemcpy(gpu_LJEnForSortAndCPUReductionValues, flatgpu_LJEn, totalNumInters * sizeof(double),
                          cudaMemcpyDeviceToDevice));

  // Determine temporary device storage requirements
void     *d_temp_storage1 = NULL;
size_t   temp_storage_bytes1 = 0;
cub::DeviceRadixSort::SortPairs(d_temp_storage1, temp_storage_bytes1,
  gpu_REnForSortAndCPUReductionKeys, gpu_REnForSortAndCPUReductionKeys, gpu_REnForSortAndCPUReductionValues, gpu_REnSortedValues, energyVectorLen);
// Allocate temporary storage
cudaMalloc(&d_temp_storage1, temp_storage_bytes1);
// Run sorting operation
cub::DeviceRadixSort::SortPairs(d_temp_storage1, temp_storage_bytes1,
  gpu_REnForSortAndCPUReductionKeys, gpu_REnForSortAndCPUReductionKeys, gpu_REnForSortAndCPUReductionValues, gpu_REnSortedValues, energyVectorLen);

// Copy sorted Keys to back to host
CubDebugExit(cudaMemcpy(hostEnergyVectorREnKeys, gpu_REnForSortAndCPUReductionKeys, totalNumInters * sizeof(int),
    cudaMemcpyDeviceToHost));
CubDebugExit(cudaMemcpy(hostEnergyVectorLJKeys, gpu_LJEnForSortAndCPUReductionKeys, totalNumInters * sizeof(int),
    cudaMemcpyDeviceToHost));

*/
// Copy sorted Values back to host
CubDebugExit(cudaMemcpy(hostEnergyVectorREnValues, flatgpu_REn, totalNumInters * sizeof(double),
    cudaMemcpyDeviceToHost));
CubDebugExit(cudaMemcpy(hostEnergyVectorLJValues, flatgpu_LJEn, totalNumInters * sizeof(double),
    cudaMemcpyDeviceToHost));

  *numberOfInters = totalNumInters;

  // ReduceSum
  void * d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_REn,
                    gpu_final_REn, energyVectorLen);
  CubDebugExit(CUMALLOC(&d_temp_storage, temp_storage_bytes));
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_REn,
                    gpu_final_REn, energyVectorLen);
  CUFREE(d_temp_storage);

  // LJ ReduceSum
  d_temp_storage = NULL;
  temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_LJEn,
                    gpu_final_LJEn, energyVectorLen);
  CubDebugExit(CUMALLOC(&d_temp_storage, temp_storage_bytes));
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_LJEn,
                    gpu_final_LJEn, energyVectorLen);
  CUFREE(d_temp_storage);
  // Copy back the result to CPU ! :)
  CubDebugExit(cudaMemcpy(&REn, gpu_final_REn, sizeof(double),
                          cudaMemcpyDeviceToHost));
  CubDebugExit(cudaMemcpy(&LJEn, gpu_final_LJEn, sizeof(double),
                          cudaMemcpyDeviceToHost));                        
/*
  std::cout <<   "gpu's box : " << box << std::endl;
       
  double cutoff = fmax(vars->gpu_rCut[0], vars->gpu_rCutCoulomb[0]);
  std::cout <<   "gpu's cutoff : " << cutoff << std::endl;
  std::cout <<   "gpu's gpu_rCut[0] : " << vars->gpu_rCut << std::endl;
  std::cout <<   "gpu's gpu_rCutCoulomb[box] : " << vars->gpu_rCutCoulomb << std::endl;*/

  CUFREE(gpu_particleCharge);
  CUFREE(gpu_particleKind);
  CUFREE(gpu_particleMol);
  CUFREE(gpu_REn);
  CUFREE(gpu_LJEn);
  CUFREE(gpu_final_REn);
  CUFREE(gpu_final_LJEn);
  CUFREE(gpu_neighborList);
  CUFREE(gpu_cellStartIndex);


  CUFREE(gpu_REnForSortAndCPUReductionKeys);
  CUFREE(gpu_LJEnForSortAndCPUReductionKeys);
  CUFREE(gpu_REnForSortAndCPUReductionValues);
  CUFREE(gpu_LJEnForSortAndCPUReductionValues);
  CUFREE(gpu_REnSortedValues);
  CUFREE(gpu_LJEnSortedValues);
}

void GetNumberOfInters(VariablesCUDA *vars,
  std::vector<int> cellVector,
  std::vector<int> cellStartIndex,
  std::vector<std::vector<int> > neighborList,
  XYZArray const &coords,
  BoxDimensions const &boxAxes,
  bool electrostatic,
  std::vector<double> particleCharge,
  std::vector<int> particleKind,
  std::vector<int> particleMol,
  double &REn,
  double &LJEn,
  bool sc_coul,
  double sc_sigma_6,
  double sc_alpha,
  uint sc_power,
  uint const box,
  uint * numberOfInters
 )
{
int atomNumber = coords.Count();
int cellVectorCount = cellVector.size();
int neighborListCount = neighborList.size() * NUMBER_OF_NEIGHBOR_CELL;
int numberOfCells = neighborList.size();
int *gpu_particleKind, *gpu_particleMol;
int *gpu_neighborList, *gpu_cellStartIndex;
int blocksPerGrid, threadsPerBlock;
int energyVectorLen = 0;
double *gpu_particleCharge;
double *gpu_REn, *gpu_LJEn;
double *gpu_final_REn, *gpu_final_LJEn;

int *gpu_REnForSortAndCPUReductionKeys, *gpu_LJEnForSortAndCPUReductionKeys;
double *gpu_REnForSortAndCPUReductionValues, *gpu_LJEnForSortAndCPUReductionValues, *gpu_REnSortedValues, *gpu_LJEnSortedValues;



// Run the kernel
threadsPerBlock = 256;
blocksPerGrid = (int)(numberOfCells * NUMBER_OF_NEIGHBOR_CELL);
energyVectorLen = numberOfCells * NUMBER_OF_NEIGHBOR_CELL * threadsPerBlock;

// Convert neighbor list to 1D array
std::vector<int> neighborlist1D(neighborListCount);
for(int i = 0; i < neighborList.size(); i++) {
for(int j = 0; j < NUMBER_OF_NEIGHBOR_CELL; j++) {
neighborlist1D[i * NUMBER_OF_NEIGHBOR_CELL + j] = neighborList[i][j];
}
}


CUMALLOC((void**) &gpu_REnForSortAndCPUReductionKeys, energyVectorLen * sizeof(int));
CUMALLOC((void**) &gpu_LJEnForSortAndCPUReductionKeys, energyVectorLen * sizeof(int));
CUMALLOC((void**) &gpu_REnForSortAndCPUReductionValues, energyVectorLen * sizeof(double));
CUMALLOC((void**) &gpu_LJEnForSortAndCPUReductionValues, energyVectorLen * sizeof(double));

CUMALLOC((void**) &gpu_REnSortedValues, energyVectorLen * sizeof(double));
CUMALLOC((void**) &gpu_LJEnSortedValues, energyVectorLen * sizeof(double));

CUMALLOC((void**) &gpu_neighborList, neighborListCount * sizeof(int));
CUMALLOC((void**) &gpu_cellStartIndex, cellStartIndex.size() * sizeof(int));
CUMALLOC((void**) &gpu_particleCharge, particleCharge.size() * sizeof(double));
CUMALLOC((void**) &gpu_particleKind, particleKind.size() * sizeof(int));
CUMALLOC((void**) &gpu_particleMol, particleMol.size() * sizeof(int));
CUMALLOC((void**) &gpu_REn, energyVectorLen * sizeof(double));
CUMALLOC((void**) &gpu_LJEn, energyVectorLen * sizeof(double));
CUMALLOC((void**) &gpu_final_REn, sizeof(double));
CUMALLOC((void**) &gpu_final_LJEn, sizeof(double));

// Copy necessary data to GPU
cudaMemcpy(gpu_neighborList, &neighborlist1D[0], neighborListCount * sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(gpu_cellStartIndex, &cellStartIndex[0], cellStartIndex.size() * sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(vars->gpu_cellVector, &cellVector[0], atomNumber * sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(gpu_particleCharge, &particleCharge[0], particleCharge.size() * sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy(gpu_particleKind, &particleKind[0], particleKind.size() * sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(gpu_particleMol, &particleMol[0], particleMol.size() * sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(double), cudaMemcpyHostToDevice);

uint * numIntersPerCell = new uint[blocksPerGrid];
uint * gpu_numIntersPerCell;
CUMALLOC((void**) &gpu_numIntersPerCell, blocksPerGrid * sizeof(uint));

CalcNumberOfInteractions <<< blocksPerGrid, threadsPerBlock>>>(gpu_cellStartIndex,
vars->gpu_cellVector,
gpu_neighborList,
numberOfCells,
cellVectorCount,
vars->gpu_x,
vars->gpu_y,
vars->gpu_z,
boxAxes.GetAxis(box).x,
boxAxes.GetAxis(box).y,
boxAxes.GetAxis(box).z,
electrostatic,
gpu_particleCharge,
gpu_particleKind,
gpu_particleMol,
gpu_REn,
gpu_LJEn,
vars->gpu_sigmaSq,
vars->gpu_epsilon_Cn,
vars->gpu_n,
vars->gpu_VDW_Kind,
vars->gpu_isMartini,
vars->gpu_count,
vars->gpu_rCut,
vars->gpu_rCutCoulomb,
vars->gpu_rCutLow,
vars->gpu_rOn,
vars->gpu_alpha,
vars->gpu_ewald,
vars->gpu_diElectric_1,
vars->gpu_nonOrth,
vars->gpu_cell_x[box],
vars->gpu_cell_y[box],
vars->gpu_cell_z[box],
vars->gpu_Invcell_x[box],
vars->gpu_Invcell_y[box],
vars->gpu_Invcell_z[box],
sc_coul,
sc_sigma_6,
sc_alpha,
sc_power,
vars->gpu_rMin,
vars->gpu_rMaxSq,
vars->gpu_expConst,
vars->gpu_molIndex,
vars->gpu_kindIndex,
vars->gpu_lambdaVDW,
vars->gpu_lambdaCoulomb,
vars->gpu_isFraction,
box,
gpu_numIntersPerCell);
cudaDeviceSynchronize();
checkLastErrorCUDA(__FILE__, __LINE__);

CubDebugExit(cudaMemcpy(numIntersPerCell, gpu_numIntersPerCell, blocksPerGrid * sizeof(uint),
cudaMemcpyDeviceToHost));


uint totalNumInters = 0;
for (int i = 0; i < blocksPerGrid; i++){
totalNumInters+=numIntersPerCell[i];
}  

*numberOfInters = totalNumInters;

CUFREE(gpu_particleCharge);
CUFREE(gpu_particleKind);
CUFREE(gpu_particleMol);
CUFREE(gpu_REn);
CUFREE(gpu_LJEn);
CUFREE(gpu_final_REn);
CUFREE(gpu_final_LJEn);
CUFREE(gpu_neighborList);
CUFREE(gpu_cellStartIndex);


CUFREE(gpu_REnForSortAndCPUReductionKeys);
CUFREE(gpu_LJEnForSortAndCPUReductionKeys);
CUFREE(gpu_REnForSortAndCPUReductionValues);
CUFREE(gpu_LJEnForSortAndCPUReductionValues);
CUFREE(gpu_REnSortedValues);
CUFREE(gpu_LJEnSortedValues);
}


__global__ void CalcNumberOfInteractions( int *gpu_cellStartIndex,
                                          int *gpu_cellVector,
                                          int *gpu_neighborList,
                                          int numberOfCells,
                                          int cellVectorCount,
                                          double *gpu_x,
                                          double *gpu_y,
                                          double *gpu_z,
                                          double xAxes,
                                          double yAxes,
                                          double zAxes,
                                          bool electrostatic,
                                          double *gpu_particleCharge,
                                          int *gpu_particleKind,
                                          int *gpu_particleMol,
                                          double *gpu_REn,
                                          double *gpu_LJEn,
                                          double *gpu_sigmaSq,
                                          double *gpu_epsilon_Cn,
                                          double *gpu_n,
                                          int *gpu_VDW_Kind,
                                          int *gpu_isMartini,
                                          int *gpu_count,
                                          double *gpu_rCut,
                                          double *gpu_rCutCoulomb,
                                          double *gpu_rCutLow,
                                          double *gpu_rOn,
                                          double *gpu_alpha,
                                          int *gpu_ewald,
                                          double *gpu_diElectric_1,
                                          int *gpu_nonOrth,
                                          double *gpu_cell_x,
                                          double *gpu_cell_y,
                                          double *gpu_cell_z,
                                          double *gpu_Invcell_x,
                                          double *gpu_Invcell_y,
                                          double *gpu_Invcell_z,
                                          bool sc_coul,
                                          double sc_sigma_6,
                                          double sc_alpha,
                                          uint sc_power,
                                          double *gpu_rMin,
                                          double *gpu_rMaxSq,
                                          double *gpu_expConst,
                                          int *gpu_molIndex,
                                          int *gpu_kindIndex,
                                          double *gpu_lambdaVDW,
                                          double *gpu_lambdaCoulomb,
                                          bool *gpu_isFraction,
                                          int box,
                                          uint * numIntersPerCell)
{

  // GJS For flattening //
  const int threadsPerBlock = 256;
  __shared__ uint cache[threadsPerBlock];
  uint myCounter = 0;

  double distSq;
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  gpu_REn[threadID] = 0.0;
  gpu_LJEn[threadID] = 0.0;
  double cutoff = fmax(gpu_rCut[0], gpu_rCutCoulomb[box]);

  int currentCell = blockIdx.x / 27;
  int nCellIndex = blockIdx.x;
  int neighborCell = gpu_neighborList[nCellIndex];

  // calculate number of particles inside neighbor Cell
  int particlesInsideCurrentCell, particlesInsideNeighboringCells;
  int endIndex = neighborCell != numberOfCells - 1 ?
                 gpu_cellStartIndex[neighborCell + 1] : cellVectorCount;
  particlesInsideNeighboringCells = endIndex - gpu_cellStartIndex[neighborCell];

  // Calculate number of particles inside current Cell
  endIndex = currentCell != numberOfCells - 1 ?
             gpu_cellStartIndex[currentCell + 1] : cellVectorCount;
  particlesInsideCurrentCell = endIndex - gpu_cellStartIndex[currentCell];

  // total number of pairs
  int numberOfPairs = particlesInsideCurrentCell * particlesInsideNeighboringCells;

  for(int pairIndex = threadIdx.x; pairIndex < numberOfPairs; pairIndex += blockDim.x) {
    int neighborParticleIndex = pairIndex / particlesInsideCurrentCell;
    int currentParticleIndex = pairIndex % particlesInsideCurrentCell;

    int currentParticle = gpu_cellVector[gpu_cellStartIndex[currentCell] + currentParticleIndex];
    int neighborParticle = gpu_cellVector[gpu_cellStartIndex[neighborCell] + neighborParticleIndex];

    if(currentParticle < neighborParticle && gpu_particleMol[currentParticle] != gpu_particleMol[neighborParticle]) {
      // Check if they are within rcut
      distSq = 0;
      double dx = gpu_x[currentParticle] - gpu_x[neighborParticle];
      double dy = gpu_y[currentParticle] - gpu_y[neighborParticle];
      double dz = gpu_z[currentParticle] - gpu_z[neighborParticle];

      dx = min(fabs(dx), xAxes - fabs(dx));
      dy = min(fabs(dy), yAxes - fabs(dy));
      dz = min(fabs(dz), zAxes - fabs(dz));

      distSq = dx * dx + dy * dy + dz * dz;

      if((cutoff * cutoff) > distSq) {
        myCounter++;
      }
    }
  }

  int cacheIndex = threadIdx.x;
  cache[cacheIndex] = myCounter;
  __syncthreads();
  int i = blockDim.x/2;
  while(i != 0){
    if(cacheIndex < i)
      cache[cacheIndex] += cache[cacheIndex+i];
    __syncthreads();
    i /= 2;
  }
  if (cacheIndex == 0)
    numIntersPerCell[blockIdx.x] = cache[0];
}

__global__ void BoxInterGPU(int *gpu_cellStartIndex,
                            int *gpu_cellVector,
                            int *gpu_neighborList,
                            int numberOfCells,
                            int cellVectorCount,
                            double *gpu_x,
                            double *gpu_y,
                            double *gpu_z,
                            double xAxes,
                            double yAxes,
                            double zAxes,
                            bool electrostatic,
                            double *gpu_particleCharge,
                            int *gpu_particleKind,
                            int *gpu_particleMol,
                            double *gpu_REn,
                            double *gpu_LJEn,
                            double *gpu_sigmaSq,
                            double *gpu_epsilon_Cn,
                            double *gpu_n,
                            int *gpu_VDW_Kind,
                            int *gpu_isMartini,
                            int *gpu_count,
                            double *gpu_rCut,
                            double *gpu_rCutCoulomb,
                            double *gpu_rCutLow,
                            double *gpu_rOn,
                            double *gpu_alpha,
                            int *gpu_ewald,
                            double *gpu_diElectric_1,
                            int *gpu_nonOrth,
                            double *gpu_cell_x,
                            double *gpu_cell_y,
                            double *gpu_cell_z,
                            double *gpu_Invcell_x,
                            double *gpu_Invcell_y,
                            double *gpu_Invcell_z,
                            bool sc_coul,
                            double sc_sigma_6,
                            double sc_alpha,
                            uint sc_power,
                            double *gpu_rMin,
                            double *gpu_rMaxSq,
                            double *gpu_expConst,
                            int *gpu_molIndex,
                            int *gpu_kindIndex,
                            double *gpu_lambdaVDW,
                            double *gpu_lambdaCoulomb,
                            bool *gpu_isFraction,
                            int box)
{
  double distSq;
  double qi_qj_fact;
  double qqFact = 167000.0;
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  gpu_REn[threadID] = 0.0;
  gpu_LJEn[threadID] = 0.0;
  double lambdaVDW = 0.0, lambdaCoulomb = 0.0;
  double cutoff = fmax(gpu_rCut[0], gpu_rCutCoulomb[box]);

  int currentCell = blockIdx.x / 27;
  int nCellIndex = blockIdx.x;
  int neighborCell = gpu_neighborList[nCellIndex];

  // calculate number of particles inside neighbor Cell
  int particlesInsideCurrentCell, particlesInsideNeighboringCells;
  int endIndex = neighborCell != numberOfCells - 1 ?
                 gpu_cellStartIndex[neighborCell + 1] : cellVectorCount;
  particlesInsideNeighboringCells = endIndex - gpu_cellStartIndex[neighborCell];

  // Calculate number of particles inside current Cell
  endIndex = currentCell != numberOfCells - 1 ?
             gpu_cellStartIndex[currentCell + 1] : cellVectorCount;
  particlesInsideCurrentCell = endIndex - gpu_cellStartIndex[currentCell];

  // total number of pairs
  int numberOfPairs = particlesInsideCurrentCell * particlesInsideNeighboringCells;

  for(int pairIndex = threadIdx.x; pairIndex < numberOfPairs; pairIndex += blockDim.x) {
    int neighborParticleIndex = pairIndex / particlesInsideCurrentCell;
    int currentParticleIndex = pairIndex % particlesInsideCurrentCell;

    int currentParticle = gpu_cellVector[gpu_cellStartIndex[currentCell] + currentParticleIndex];
    int neighborParticle = gpu_cellVector[gpu_cellStartIndex[neighborCell] + neighborParticleIndex];

    if(currentParticle < neighborParticle && gpu_particleMol[currentParticle] != gpu_particleMol[neighborParticle]) {
      // Check if they are within rcut
      distSq = 0;
      double dx = gpu_x[currentParticle] - gpu_x[neighborParticle];
      double dy = gpu_y[currentParticle] - gpu_y[neighborParticle];
      double dz = gpu_z[currentParticle] - gpu_z[neighborParticle];

      dx = min(fabs(dx), xAxes - fabs(dx));
      dy = min(fabs(dy), yAxes - fabs(dy));
      dz = min(fabs(dz), zAxes - fabs(dz));

      distSq = dx * dx + dy * dy + dz * dz;

      if((cutoff * cutoff) > distSq) {
        double cA = gpu_particleCharge[currentParticle];
        double cB = gpu_particleCharge[neighborParticle];
        int kA = gpu_particleKind[currentParticle];
        int kB = gpu_particleKind[neighborParticle];
        int mA = gpu_particleMol[currentParticle];
        int mB = gpu_particleMol[neighborParticle];

        lambdaVDW = DeviceGetLambdaVDW(mA, kA, mB, kB, box, gpu_isFraction,
                                       gpu_molIndex, gpu_kindIndex, gpu_lambdaVDW);

        if(electrostatic) {
          qi_qj_fact = cA * cB * qqFact;
          lambdaCoulomb = DeviceGetLambdaCoulomb(mA, kA, mB, kB, box,
                                                 gpu_isFraction, gpu_molIndex,
                                                 gpu_kindIndex, gpu_lambdaCoulomb);
          gpu_REn[threadID] += CalcCoulombGPU(distSq, kA, kB,
                                              qi_qj_fact, gpu_rCutLow[0],
                                              gpu_ewald[0], gpu_VDW_Kind[0],
                                              gpu_alpha[box],
                                              gpu_rCutCoulomb[box],
                                              gpu_isMartini[0],
                                              gpu_diElectric_1[0],
                                              lambdaCoulomb,
                                              sc_coul,
                                              sc_sigma_6,
                                              sc_alpha,
                                              sc_power,
                                              gpu_sigmaSq[box],
                                              gpu_count[0]);
        }
        gpu_LJEn[threadID] += CalcEnGPU(distSq, kA, kB,
                                        gpu_sigmaSq, gpu_n, gpu_epsilon_Cn,
                                        gpu_VDW_Kind[0], gpu_isMartini[0],
                                        gpu_rCut[0], gpu_rOn[0], gpu_count[0], lambdaVDW,
                                        sc_sigma_6, sc_alpha, sc_power, gpu_rMin,
                                        gpu_rMaxSq, gpu_expConst);
      }
    }
  }
}

__global__ void BoxInterGPUFlattened(int *gpu_cellStartIndex,
                            int *gpu_cellVector,
                            int *gpu_neighborList,
                            int numberOfCells,
                            int cellVectorCount,
                            double *gpu_x,
                            double *gpu_y,
                            double *gpu_z,
                            double xAxes,
                            double yAxes,
                            double zAxes,
                            bool electrostatic,
                            double *gpu_particleCharge,
                            int *gpu_particleKind,
                            int *gpu_particleMol,
                            double *gpu_REn,
                            double *gpu_LJEn,
                            double *gpu_sigmaSq,
                            double *gpu_epsilon_Cn,
                            double *gpu_n,
                            int *gpu_VDW_Kind,
                            int *gpu_isMartini,
                            int *gpu_count,
                            double *gpu_rCut,
                            double *gpu_rCutCoulomb,
                            double *gpu_rCutLow,
                            double *gpu_rOn,
                            double *gpu_alpha,
                            int *gpu_ewald,
                            double *gpu_diElectric_1,
                            int *gpu_nonOrth,
                            double *gpu_cell_x,
                            double *gpu_cell_y,
                            double *gpu_cell_z,
                            double *gpu_Invcell_x,
                            double *gpu_Invcell_y,
                            double *gpu_Invcell_z,
                            bool sc_coul,
                            double sc_sigma_6,
                            double sc_alpha,
                            uint sc_power,
                            double *gpu_rMin,
                            double *gpu_rMaxSq,
                            double *gpu_expConst,
                            int *gpu_molIndex,
                            int *gpu_kindIndex,
                            double *gpu_lambdaVDW,
                            double *gpu_lambdaCoulomb,
                            bool *gpu_isFraction,
                            int box,
                            uint * numIntersPerCell, 
                            double *flatgpu_REn,
                            double *flatgpu_LJEn,
                            int * flatIndexREn,
                            int * flatIndexLJEn
                          )
{
  double distSq;
  double qi_qj_fact;
  double qqFact = 167000.0;
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  gpu_REn[threadID] = 0.0;
  gpu_LJEn[threadID] = 0.0;
  double lambdaVDW = 0.0, lambdaCoulomb = 0.0;
  double cutoff = fmax(gpu_rCut[0], gpu_rCutCoulomb[box]);

  int currentCell = blockIdx.x / 27;
  int nCellIndex = blockIdx.x;
  int neighborCell = gpu_neighborList[nCellIndex];

  // calculate number of particles inside neighbor Cell
  int particlesInsideCurrentCell, particlesInsideNeighboringCells;
  int endIndex = neighborCell != numberOfCells - 1 ?
  gpu_cellStartIndex[neighborCell + 1] : cellVectorCount;
  particlesInsideNeighboringCells = endIndex - gpu_cellStartIndex[neighborCell];

  // Calculate number of particles inside current Cell
  endIndex = currentCell != numberOfCells - 1 ?
  gpu_cellStartIndex[currentCell + 1] : cellVectorCount;
  particlesInsideCurrentCell = endIndex - gpu_cellStartIndex[currentCell];

  // total number of pairs aka entries in our NxM matrix or r_ij
  int numberOfPairs = particlesInsideCurrentCell * particlesInsideNeighboringCells;

  for(int pairIndex = threadIdx.x; pairIndex < numberOfPairs; pairIndex += blockDim.x) {
  int neighborParticleIndex = pairIndex / particlesInsideCurrentCell;
  int currentParticleIndex = pairIndex % particlesInsideCurrentCell;

  int currentParticle = gpu_cellVector[gpu_cellStartIndex[currentCell] + currentParticleIndex];
  int neighborParticle = gpu_cellVector[gpu_cellStartIndex[neighborCell] + neighborParticleIndex];

    if(currentParticle < neighborParticle && gpu_particleMol[currentParticle] != gpu_particleMol[neighborParticle]) {
      // Check if they are within rcut
      distSq = 0;
      double dx = gpu_x[currentParticle] - gpu_x[neighborParticle];
      double dy = gpu_y[currentParticle] - gpu_y[neighborParticle];
      double dz = gpu_z[currentParticle] - gpu_z[neighborParticle];

      dx = min(fabs(dx), xAxes - fabs(dx));
      dy = min(fabs(dy), yAxes - fabs(dy));
      dz = min(fabs(dz), zAxes - fabs(dz));

      distSq = dx * dx + dy * dy + dz * dz;

      if((cutoff * cutoff) > distSq) {
        double cA = gpu_particleCharge[currentParticle];
        double cB = gpu_particleCharge[neighborParticle];
        int kA = gpu_particleKind[currentParticle];
        int kB = gpu_particleKind[neighborParticle];
        int mA = gpu_particleMol[currentParticle];
        int mB = gpu_particleMol[neighborParticle];

        lambdaVDW = DeviceGetLambdaVDW(mA, kA, mB, kB, box, gpu_isFraction,
                    gpu_molIndex, gpu_kindIndex, gpu_lambdaVDW);

        if(electrostatic) {
          qi_qj_fact = cA * cB * qqFact;
          lambdaCoulomb = DeviceGetLambdaCoulomb(mA, kA, mB, kB, box,
                                gpu_isFraction, gpu_molIndex,
                                gpu_kindIndex, gpu_lambdaCoulomb);
          gpu_REn[threadID] += CalcCoulombGPU(distSq, kA, kB,
                              qi_qj_fact, gpu_rCutLow[0],
                              gpu_ewald[0], gpu_VDW_Kind[0],
                              gpu_alpha[box],
                              gpu_rCutCoulomb[box],
                              gpu_isMartini[0],
                              gpu_diElectric_1[0],
                              lambdaCoulomb,
                              sc_coul,
                              sc_sigma_6,
                              sc_alpha,
                              sc_power,
                              gpu_sigmaSq[box],
                              gpu_count[0]);
              flatgpu_REn[atomicAdd(flatIndexREn, 1)] = CalcCoulombGPU( distSq, kA, kB,
                                                        qi_qj_fact, gpu_rCutLow[0],
                                                        gpu_ewald[0], gpu_VDW_Kind[0],
                                                        gpu_alpha[box],
                                                        gpu_rCutCoulomb[box],
                                                        gpu_isMartini[0],
                                                        gpu_diElectric_1[0],
                                                        lambdaCoulomb,
                                                        sc_coul,
                                                        sc_sigma_6,
                                                        sc_alpha,
                                                        sc_power,
                                                        gpu_sigmaSq[box],
                                                        gpu_count[0]);
                                                       
            
        }
        gpu_LJEn[threadID] += CalcEnGPU(distSq, kA, kB,
                      gpu_sigmaSq, gpu_n, gpu_epsilon_Cn,
                      gpu_VDW_Kind[0], gpu_isMartini[0],
                      gpu_rCut[0], gpu_rOn[0], gpu_count[0], lambdaVDW,
                      sc_sigma_6, sc_alpha, sc_power, gpu_rMin,
                      gpu_rMaxSq, gpu_expConst);
          flatgpu_LJEn[atomicAdd(flatIndexLJEn, 1)] = CalcEnGPU(distSq, kA, kB,
                      gpu_sigmaSq, gpu_n, gpu_epsilon_Cn,
                      gpu_VDW_Kind[0], gpu_isMartini[0],
                      gpu_rCut[0], gpu_rOn[0], gpu_count[0], lambdaVDW,
                      sc_sigma_6, sc_alpha, sc_power, gpu_rMin,
                      gpu_rMaxSq, gpu_expConst);
      }
    }
    //flatThreadID += blockDim.x;
  }
}

__device__ double CalcCoulombGPU(double distSq,
                                 int kind1,
                                 int kind2,
                                 double qi_qj_fact,
                                 double gpu_rCutLow,
                                 int gpu_ewald,
                                 int gpu_VDW_Kind,
                                 double gpu_alpha,
                                 double gpu_rCutCoulomb,
                                 int gpu_isMartini,
                                 double gpu_diElectric_1,
                                 double gpu_lambdaCoulomb,
                                 bool sc_coul,
                                 double sc_sigma_6,
                                 double sc_alpha,
                                 uint sc_power,
                                 double gpu_sigmaSq,
                                 int gpu_count)
{
  if((gpu_rCutCoulomb * gpu_rCutCoulomb) < distSq) {
    return 0.0;
  }

  int index = FlatIndexGPU(kind1, kind2, gpu_count);
  if(gpu_VDW_Kind == GPU_VDW_STD_KIND) {
    return CalcCoulombParticleGPU(distSq, qi_qj_fact, gpu_ewald, gpu_alpha,
                                  gpu_lambdaCoulomb, sc_coul, sc_sigma_6,
                                  sc_alpha, sc_power, gpu_sigmaSq);
  } else if(gpu_VDW_Kind == GPU_VDW_SHIFT_KIND) {
    return CalcCoulombShiftGPU(distSq, qi_qj_fact, gpu_ewald, gpu_alpha,
                               gpu_rCutCoulomb, gpu_lambdaCoulomb, sc_coul,
                               sc_sigma_6, sc_alpha, sc_power, gpu_sigmaSq);
  } else if(gpu_VDW_Kind == GPU_VDW_EXP6_KIND) {
    return CalcCoulombExp6GPU(distSq, qi_qj_fact, gpu_ewald, gpu_alpha,
                              gpu_lambdaCoulomb, sc_coul, sc_sigma_6, sc_alpha,
                              sc_power, gpu_sigmaSq);
  } else if(gpu_VDW_Kind == GPU_VDW_SWITCH_KIND && gpu_isMartini) {
    return CalcCoulombSwitchMartiniGPU(distSq, qi_qj_fact, gpu_ewald, gpu_alpha,
                                       gpu_rCutCoulomb, gpu_diElectric_1,
                                       gpu_lambdaCoulomb, sc_coul, sc_sigma_6,
                                       sc_alpha, sc_power, gpu_sigmaSq);
  } else
    return CalcCoulombSwitchGPU(distSq, qi_qj_fact, gpu_alpha, gpu_ewald,
                                gpu_rCutCoulomb, gpu_lambdaCoulomb,
                                sc_coul, sc_sigma_6, sc_alpha, sc_power,
                                gpu_sigmaSq);
}

__device__ double CalcEnGPU(double distSq, int kind1, int kind2,
                            double *gpu_sigmaSq, double *gpu_n,
                            double *gpu_epsilon_Cn, int gpu_VDW_Kind,
                            int gpu_isMartini, double gpu_rCut, double gpu_rOn,
                            int gpu_count, double gpu_lambdaVDW,
                            double sc_sigma_6, double sc_alpha, uint sc_power,
                            double *gpu_rMin, double *gpu_rMaxSq,
                            double *gpu_expConst)
{
  if((gpu_rCut * gpu_rCut) < distSq) {
    return 0.0;
  }

  int index = FlatIndexGPU(kind1, kind2, gpu_count);
  if(gpu_VDW_Kind == GPU_VDW_STD_KIND) {
    return CalcEnParticleGPU(distSq, index, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn,
                             gpu_lambdaVDW, sc_sigma_6, sc_alpha, sc_power);
  } else if(gpu_VDW_Kind == GPU_VDW_SHIFT_KIND) {
    return CalcEnShiftGPU(distSq, index, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn,
                          gpu_rCut, gpu_lambdaVDW, sc_sigma_6, sc_alpha,
                          sc_power);
  } else if(gpu_VDW_Kind == GPU_VDW_EXP6_KIND) {
    return CalcEnExp6GPU(distSq, index, gpu_sigmaSq[index], gpu_n[index],
                         gpu_lambdaVDW, sc_sigma_6,
                         sc_alpha, sc_power, gpu_rMin[index],
                         gpu_rMaxSq[index], gpu_expConst[index]);
  } else if(gpu_VDW_Kind == GPU_VDW_SWITCH_KIND && gpu_isMartini) {
    return CalcEnSwitchMartiniGPU(distSq, index, gpu_sigmaSq, gpu_n,
                                  gpu_epsilon_Cn, gpu_rCut, gpu_rOn,
                                  gpu_lambdaVDW, sc_sigma_6, sc_alpha,
                                  sc_power);
  } else
    return CalcEnSwitchGPU(distSq, index, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn,
                           gpu_rCut, gpu_rOn, gpu_lambdaVDW, sc_sigma_6,
                           sc_alpha, sc_power);
}

//ElectroStatic Calculation
//**************************************************************//
__device__ double CalcCoulombParticleGPU(double distSq, double qi_qj_fact,
    double gpu_ewald, double gpu_alpha,
    double gpu_lambdaCoulomb, bool sc_coul,
    double sc_sigma_6, double sc_alpha,
    uint sc_power, double gpu_sigmaSq)
{
  if(gpu_lambdaCoulomb >= 0.999999) {
    return CalcCoulombParticleGPUNoLambda(distSq, qi_qj_fact, gpu_ewald, gpu_alpha);
  }
  if(sc_coul) {
    double sigma6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
    sigma6 = max(sigma6, sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaCoulomb), (double)sc_power);
    double softDist6 = lambdaCoef * sigma6 * dist6;
    double softRsq = pow(softDist6, 1.0 / 3.0);
    return gpu_lambdaCoulomb * CalcCoulombParticleGPUNoLambda(softRsq, qi_qj_fact, gpu_ewald, gpu_alpha);
  } else {
    return gpu_lambdaCoulomb * CalcCoulombParticleGPUNoLambda(distSq, qi_qj_fact, gpu_ewald, gpu_alpha);
  }
}

__device__ double CalcCoulombParticleGPUNoLambda(double distSq,
    double qi_qj_fact,
    double gpu_ewald,
    double gpu_alpha)
{
  if(gpu_ewald) {
    double dist = sqrt(distSq);
    double value = gpu_alpha * dist;
    return qi_qj_fact * erfc(value) / dist;
  } else {
    double dist = sqrt(distSq);
    return qi_qj_fact / dist;
  }
}

__device__ double CalcCoulombShiftGPU(double distSq, double qi_qj_fact,
                                      int gpu_ewald, double gpu_alpha,
                                      double gpu_rCut, double gpu_lambdaCoulomb,
                                      bool sc_coul, double sc_sigma_6,
                                      double sc_alpha, uint sc_power,
                                      double gpu_sigmaSq)
{

  if(gpu_lambdaCoulomb >= 0.999999) {
    return CalcCoulombShiftGPUNoLambda(distSq, qi_qj_fact, gpu_ewald, gpu_alpha,
                                       gpu_rCut);
  }

  if(sc_coul) {
    double sigma6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
    sigma6 = max(sigma6, sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaCoulomb), (double)sc_power);
    double softDist6 = lambdaCoef * sigma6 * dist6;
    double softRsq = pow(softDist6, (double)1.0 / 3.0);
    return gpu_lambdaCoulomb * CalcCoulombShiftGPUNoLambda(softRsq, qi_qj_fact,
           gpu_ewald, gpu_alpha,
           gpu_rCut);
  } else {
    return gpu_lambdaCoulomb * CalcCoulombShiftGPUNoLambda(distSq, qi_qj_fact,
           gpu_ewald, gpu_alpha,
           gpu_rCut);
  }
}

__device__ double CalcCoulombShiftGPUNoLambda(double distSq, double qi_qj_fact,
    int gpu_ewald, double gpu_alpha,
    double gpu_rCut)
{
  if(gpu_ewald) {
    double dist = sqrt(distSq);
    double value = gpu_alpha * dist;
    return qi_qj_fact * (1 - erf(value)) / dist;
  } else {
    double dist = sqrt(distSq);
    return qi_qj_fact * (1.0 / dist - 1.0 / gpu_rCut);
  }
}

__device__ double CalcCoulombExp6GPU(double distSq, double qi_qj_fact,
                                     int gpu_ewald, double gpu_alpha,
                                     double gpu_lambdaCoulomb, bool sc_coul,
                                     double sc_sigma_6, double sc_alpha,
                                     uint sc_power, double gpu_sigmaSq)
{
  if(gpu_lambdaCoulomb >= 0.999999) {
    return CalcCoulombExp6GPUNoLambda(distSq, qi_qj_fact, gpu_ewald, gpu_alpha);
  }

  if(sc_coul) {
    double sigma6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
    sigma6 = max(sigma6, sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaCoulomb),
                                       (double)sc_power);
    double softDist6 = lambdaCoef * sigma6 * dist6;
    double softRsq = pow(softDist6, (double)1.0 / 3.0);
    return gpu_lambdaCoulomb * CalcCoulombExp6GPUNoLambda(softRsq, qi_qj_fact,
           gpu_ewald, gpu_alpha);
  } else {
    return gpu_lambdaCoulomb * CalcCoulombExp6GPUNoLambda(distSq, qi_qj_fact,
           gpu_ewald, gpu_alpha);
  }
}

__device__ double CalcCoulombExp6GPUNoLambda(double distSq, double qi_qj_fact,
    int gpu_ewald, double gpu_alpha)
{
  if(gpu_ewald) {
    double dist = sqrt(distSq);
    double val = gpu_alpha * dist;
    return qi_qj_fact * erfc(val) / dist;
  } else {
    double dist = sqrt(distSq);
    return qi_qj_fact / dist;
  }
}

__device__ double CalcCoulombSwitchMartiniGPU(double distSq, double qi_qj_fact,
    int gpu_ewald, double gpu_alpha,
    double gpu_rCut,
    double gpu_diElectric_1,
    double gpu_lambdaCoulomb,
    bool sc_coul, double sc_sigma_6,
    double sc_alpha, uint sc_power,
    double gpu_sigmaSq)
{
  if(gpu_lambdaCoulomb >= 0.999999) {
    return CalcCoulombSwitchMartiniGPUNoLambda(distSq, qi_qj_fact, gpu_ewald, gpu_alpha, gpu_rCut, gpu_diElectric_1);
  }

  if(sc_coul) {
    double sigma6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
    sigma6 = max(sigma6, sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaCoulomb), (double)sc_power);
    double softDist6 = lambdaCoef * sigma6 * dist6;
    double softRsq = pow(softDist6, 1.0 / 3.0);
    return gpu_lambdaCoulomb * CalcCoulombSwitchMartiniGPUNoLambda(softRsq, qi_qj_fact, gpu_ewald, gpu_alpha, gpu_rCut, gpu_diElectric_1);
  } else {
    return gpu_lambdaCoulomb * CalcCoulombSwitchMartiniGPUNoLambda(distSq, qi_qj_fact, gpu_ewald, gpu_alpha, gpu_rCut, gpu_diElectric_1);
  }
}

__device__ double CalcCoulombSwitchMartiniGPUNoLambda(double distSq,
    double qi_qj_fact,
    int gpu_ewald,
    double gpu_alpha,
    double gpu_rCut,
    double gpu_diElectric_1)
{
  if(gpu_ewald) {
    double dist = sqrt(distSq);
    double value = gpu_alpha * dist;
    return qi_qj_fact * (1 - erf(value)) / dist;
  } else {
    // in Martini, the Coulomb switching distance is zero, so we will have
    // sqrt(distSq) - rOnCoul =  sqrt(distSq)
    double dist = sqrt(distSq);
    double rij_ronCoul_3 = dist * distSq;
    double rij_ronCoul_4 = distSq * distSq;

    double A1 = 1.0 * (-(1.0 + 4) * gpu_rCut) / (pow(gpu_rCut, 1.0 + 2) *
                pow(gpu_rCut, 2));
    double B1 = -1.0 * (-(1.0 + 3) * gpu_rCut) / (pow(gpu_rCut, 1.0 + 2) *
                pow(gpu_rCut, 3));
    double C1 = 1.0 / pow(gpu_rCut, 1.0) - A1 / 3.0 * pow(gpu_rCut, 3) -
                B1 / 4.0 * pow(gpu_rCut, 4);

    double coul = -(A1 / 3.0) * rij_ronCoul_3 - (B1 / 4.0) * rij_ronCoul_4 - C1;
    return qi_qj_fact * gpu_diElectric_1 * (1.0 / dist + coul);
  }
}

__device__ double CalcCoulombSwitchGPU(double distSq, double qi_qj_fact,
                                       double gpu_alpha, int gpu_ewald,
                                       double gpu_rCut,
                                       double gpu_lambdaCoulomb, bool sc_coul,
                                       double sc_sigma_6, double sc_alpha,
                                       uint sc_power, double gpu_sigmaSq)
{
  if(gpu_lambdaCoulomb >= 0.999999) {
    return CalcCoulombSwitchGPUNoLambda(distSq, qi_qj_fact, gpu_ewald, gpu_alpha, gpu_rCut);
  }

  if(sc_coul) {
    double sigma6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
    sigma6 = max(sigma6, sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaCoulomb), (double)sc_power);
    double softDist6 = lambdaCoef * sigma6 * dist6;
    double softRsq = pow(softDist6, 1.0 / 3.0);
    return gpu_lambdaCoulomb * CalcCoulombSwitchGPUNoLambda(softRsq, qi_qj_fact, gpu_ewald, gpu_alpha, gpu_rCut);
  } else {
    return gpu_lambdaCoulomb * CalcCoulombSwitchGPUNoLambda(distSq, qi_qj_fact, gpu_ewald, gpu_alpha, gpu_rCut);
  }
}

__device__ double CalcCoulombSwitchGPUNoLambda(double distSq, double qi_qj_fact,
    double gpu_alpha, int gpu_ewald,
    double gpu_rCut)
{
  if(gpu_ewald) {
    double dist = sqrt(distSq);
    double value = gpu_alpha * dist;
    return qi_qj_fact * (1 - erf(value)) / dist;
  } else {
    double rCutSq = gpu_rCut * gpu_rCut;
    double dist = sqrt(distSq);
    double switchVal = distSq / rCutSq - 1.0;
    switchVal *= switchVal;
    return qi_qj_fact * switchVal / dist;
  }
}

//VDW Calculation
//**************************************************************//
__device__ double CalcEnParticleGPU(double distSq, int index,
                                    double *gpu_sigmaSq, double *gpu_n,
                                    double *gpu_epsilon_Cn,
                                    double gpu_lambdaVDW,
                                    double sc_sigma_6,
                                    double sc_alpha,
                                    uint sc_power)
{
  if(gpu_lambdaVDW >= 0.999999) {
    return CalcEnParticleGPUNoLambda(distSq, index, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn);
  }

  double sigma6 = gpu_sigmaSq[index] * gpu_sigmaSq[index] * gpu_sigmaSq[index];
  sigma6 = max(sigma6, sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaVDW), (double)sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = pow(softDist6, 1.0 / 3.0);

  return gpu_lambdaVDW * CalcEnParticleGPUNoLambda(softRsq, index, gpu_sigmaSq,
         gpu_n, gpu_epsilon_Cn);
}

__device__ double CalcEnParticleGPUNoLambda(double distSq, int index,
    double *gpu_sigmaSq, double *gpu_n,
    double *gpu_epsilon_Cn)
{
  double rRat2 = gpu_sigmaSq[index] / distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double repulse = pow(rRat2, gpu_n[index] / 2.0);
  return gpu_epsilon_Cn[index] * (repulse - attract);
}

__device__ double CalcEnShiftGPU(double distSq, int index, double *gpu_sigmaSq,
                                 double *gpu_n, double *gpu_epsilon_Cn,
                                 double gpu_rCut,
                                 double gpu_lambdaVDW,
                                 double sc_sigma_6,
                                 double sc_alpha,
                                 uint sc_power)
{
  if(gpu_lambdaVDW >= 0.999999) {
    return CalcEnShiftGPUNoLambda(distSq, index, gpu_sigmaSq, gpu_n,
                                  gpu_epsilon_Cn, gpu_rCut);
  }

  double sigma6 = gpu_sigmaSq[index] * gpu_sigmaSq[index] * gpu_sigmaSq[index];
  sigma6 = max(sigma6, sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaVDW), (double)sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = pow(softDist6, (double)1.0 / 3.0);

  return gpu_lambdaVDW * CalcEnShiftGPUNoLambda(softRsq, index, gpu_sigmaSq,
         gpu_n, gpu_epsilon_Cn,
         gpu_rCut);
}

__device__ double CalcEnShiftGPUNoLambda(double distSq, int index,
    double *gpu_sigmaSq,
    double *gpu_n, double *gpu_epsilon_Cn,
    double gpu_rCut)
{
  double rRat2 = gpu_sigmaSq[index] / distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double repulse = pow(rRat2, gpu_n[index] / 2.0);

  double shiftRRat2 = gpu_sigmaSq[index] / (gpu_rCut * gpu_rCut);
  double shiftRRat4 = shiftRRat2 * shiftRRat2;
  double shiftAttract = shiftRRat4 * shiftRRat2;
  double shiftRepulse = pow(shiftRRat2, gpu_n[index] / 2.0);
  double shiftConst = gpu_epsilon_Cn[index] * (shiftRepulse - shiftAttract);

  return (gpu_epsilon_Cn[index] * (repulse - attract) - shiftConst);
}

__device__ double CalcEnExp6GPU(double distSq, int index, double gpu_sigmaSq,
                                double gpu_n, double gpu_lambdaVDW,
                                double sc_sigma_6, double sc_alpha,
                                uint sc_power, double gpu_rMin,
                                double gpu_rMaxSq, double gpu_expConst)
{
  if(distSq < gpu_rMaxSq) {
    return num::BIGNUM;
  }
  if(gpu_lambdaVDW >= 0.999999) {
    return CalcEnExp6GPUNoLambda(distSq, gpu_n, gpu_rMin, gpu_expConst);
  }
  double sigma6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
  sigma6 = max(sigma6, sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaVDW), (double)sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = pow(softDist6, (double)1.0 / 3.0);

  return gpu_lambdaVDW * CalcEnExp6GPUNoLambda(softRsq,  gpu_n, gpu_rMin,
         gpu_expConst);
}

__device__ double CalcEnExp6GPUNoLambda(double distSq, double gpu_n,
                                        double gpu_rMin, double gpu_expConst)
{
  double dist = sqrt(distSq);
  double rRat = gpu_rMin / dist;
  double rRat2 = rRat * rRat;
  double attract = rRat2 * rRat2 * rRat2;

  uint alph_ij = gpu_n;
  double repulse = (6.0 / alph_ij) * exp(alph_ij * (1.0 - dist / gpu_rMin));
  return gpu_expConst * (repulse - attract);
}

__device__ double CalcEnSwitchMartiniGPU(double distSq, int index,
    double *gpu_sigmaSq, double *gpu_n,
    double *gpu_epsilon_Cn,
    double gpu_rCut, double gpu_rOn,
    double gpu_lambdaVDW,
    double sc_sigma_6,
    double sc_alpha,
    uint sc_power)
{
  if(gpu_lambdaVDW >= 0.999999) {
    return CalcEnSwitchMartiniGPUNoLambda(distSq, index, gpu_sigmaSq, gpu_n,
                                          gpu_epsilon_Cn, gpu_rCut, gpu_rOn);
  }

  double sigma6 = gpu_sigmaSq[index] * gpu_sigmaSq[index] * gpu_sigmaSq[index];
  sigma6 = max(sigma6, sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaVDW), (double)sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = pow(softDist6, (double)1.0 / 3.0);

  return gpu_lambdaVDW * CalcEnSwitchMartiniGPUNoLambda(softRsq, index,
         gpu_sigmaSq, gpu_n,
         gpu_epsilon_Cn,
         gpu_rCut, gpu_rOn);
}

__device__ double CalcEnSwitchMartiniGPUNoLambda(double distSq, int index,
    double *gpu_sigmaSq,
    double *gpu_n,
    double *gpu_epsilon_Cn,
    double gpu_rCut,
    double gpu_rOn)
{
  double r_2 = 1.0 / distSq;
  double r_4 = r_2 * r_2;
  double r_6 = r_4 * r_2;
  double r_n = pow(r_2, gpu_n[index] / 2.0);

  double rij_ron = sqrt(distSq) - gpu_rOn;
  double rij_ron_2 = rij_ron * rij_ron;
  double rij_ron_3 = rij_ron_2 * rij_ron;
  double rij_ron_4 = rij_ron_2 * rij_ron_2;

  double pn = gpu_n[index];
  double An = pn * ((pn + 1) * gpu_rOn - (pn + 4) * gpu_rCut) /
              (pow(gpu_rCut, pn + 2) * pow(gpu_rCut - gpu_rOn, 2));
  double Bn = -pn * ((pn + 1) * gpu_rOn - (pn + 3) * gpu_rCut) /
              (pow(gpu_rCut, pn + 2) * pow(gpu_rCut - gpu_rOn, 3));
  double Cn = 1.0 / pow(gpu_rCut, pn) - An / 3.0 * pow(gpu_rCut - gpu_rOn, 3) -
              Bn / 4.0 * pow(gpu_rCut - gpu_rOn, 4);

  double A6 = 6.0 * ((6.0 + 1) * gpu_rOn - (6.0 + 4) * gpu_rCut) /
              (pow(gpu_rCut, 6.0 + 2) * pow(gpu_rCut - gpu_rOn, 2));
  double B6 = -6.0 * ((6.0 + 1) * gpu_rOn - (6.0 + 3) * gpu_rCut) /
              (pow(gpu_rCut, 6.0 + 2) * pow(gpu_rCut - gpu_rOn, 3));
  double C6 = 1.0 / pow(gpu_rCut, 6.0) - A6 / 3.0 * pow(gpu_rCut - gpu_rOn, 3) -
              B6 / 4.0 * pow(gpu_rCut - gpu_rOn, 4);

  double shifttempRep = -(An / 3.0) * rij_ron_3 -
                        (Bn / 4.0) * rij_ron_4 - Cn;
  double shifttempAtt = -(A6 / 3.0) * rij_ron_3 - (B6 / 4.0) * rij_ron_4 - C6;

  const double shiftRep = ( distSq > gpu_rOn * gpu_rOn ? shifttempRep : -Cn);
  const double shiftAtt = ( distSq > gpu_rOn * gpu_rOn ? shifttempAtt : -C6);

  double sig6 = pow(gpu_sigmaSq[index], 3);
  double sign = pow(gpu_sigmaSq[index], pn / 2);
  double Eij = gpu_epsilon_Cn[index] * (sign * (r_n + shiftRep) -
                                        sig6 * (r_6 + shiftAtt));
  return Eij;
}

__device__ double CalcEnSwitchGPU(double distSq, int index, double *gpu_sigmaSq,
                                  double *gpu_n, double *gpu_epsilon_Cn,
                                  double gpu_rCut, double gpu_rOn,
                                  double gpu_lambdaVDW, double sc_sigma_6,
                                  double sc_alpha, uint sc_power)
{
  if(gpu_lambdaVDW >= 0.999999) {
    return CalcEnSwitchGPUNoLambda(distSq, index, gpu_sigmaSq, gpu_n,
                                   gpu_epsilon_Cn, gpu_rCut, gpu_rOn);
  }
  double sigma6 = gpu_sigmaSq[index] * gpu_sigmaSq[index] * gpu_sigmaSq[index];
  sigma6 = max(sigma6, sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaVDW), (double)sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = pow(softDist6, (double)1.0 / 3.0);

  return gpu_lambdaVDW * CalcEnSwitchGPUNoLambda(softRsq, index, gpu_sigmaSq,
         gpu_n, gpu_epsilon_Cn,
         gpu_rCut, gpu_rOn);
}

__device__ double CalcEnSwitchGPUNoLambda(double distSq, int index,
    double *gpu_sigmaSq, double *gpu_n,
    double *gpu_epsilon_Cn,
    double gpu_rCut, double gpu_rOn)
{
  double rCutSq = gpu_rCut * gpu_rCut;
  double rOnSq = gpu_rOn * gpu_rOn;

  double rCutSq_rijSq = rCutSq  - distSq;
  double rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;

  double rRat2 = gpu_sigmaSq[index] / distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;

  double repulse = pow(rRat2, gpu_n[index] / 2.0);

  double factor1 = rCutSq - 3 * rOnSq;
  double factor2 = pow((rCutSq - rOnSq), -3);
  double fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2 * distSq);

  const double factE = ( distSq > rOnSq ? fE : 1.0);

  return (gpu_epsilon_Cn[index] * (repulse - attract)) * factE;
}

struct Lock {
  int * mutex;
  Lock (void) {
    int state = 0;
    CUMALLOC((void**)&mutex, sizeof(int));
    cudaMemcpy(mutex, &state, sizeof(int), cudaMemcpyHostToDevice);
  }

  ~Lock (void){
    cudaFree(mutex);
  }


__device__ void lock( void ) {
  while( atomicCAS( mutex, 0, 1 ) != 0 );
  __threadfence();}
  
  __device__ void unlock( void ) {
    __threadfence();
    atomicExch( mutex, 0 );
  }
};

#endif
