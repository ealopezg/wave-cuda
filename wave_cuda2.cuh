#ifndef WAVE_H
#define WAVE_H

float *initializeGrid(int N);
__device__ void schroedinger_(float *grid_0, float *grid_1, float *grid_2, int N, int i, int j, int t);
void printGrid(float *grid, int N);
__global__ void schroedinger(float *grid_0, float *grid_1, float *grid_2, int N, int T, int t);
void freeGrid(float *grid, int N);
void save(float *grid, int N, char *filename);
#endif