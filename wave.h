#ifndef WAVE_H
#define WAVE_H

float ***initializeGrid(int N);
void schroedinger_(float ***grid, int N, int i, int j, int t);
void printGrid(float ***grid, int N);
void schroedinger(float ***grid, int N, int T);
void freeGrid(float ***grid, int N);
void save(float ***grid, int N, char *filename);
#endif