#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <time.h>

#include "wave_cuda2.cuh"

/**
 * @brief Estructura donde se guardan los parámetros de entrada
 *
 */
typedef struct
{
    int grid_size;  // Tamaño de la matriz
    int steps;      // Cantidad de iteraciones
    int x;
    int y;
    char *filename; // Nombre archivo de salida
} options_t;

/**
 * @brief Crea una matriz de 3 dimensiones , N*N*3, donde la tercera dimensión corresponde a los pasos anteriores (0: t,1: t-1,2: t-2)
 *
 * @param N Tamaño de la matriz
 * @return float*** Matriz
 */
float* initializeGrid(int N){
    float* grid = (float*)malloc(N*N*sizeof(float));
    for(int i = 0;i <N*N;i++){
        grid[i] = 0;
    }
    return grid;
}


// Calcula la iteración t en cada caso
/**
 * @brief Calcula la iteración en cada caso
 *
 * @param grid Matriz de flotantes de NxN
 * @param N Tamaño de la matriz
 * @param i Posición fila
 * @param j Posición Columna
 * @param t Número Iteración
 */
__device__ void schroedinger_(float *grid_0, float *grid_1, float *grid_2, int N, int i, int j, int t)
{

    switch (t)
    {
    case 0:
        // Tiene que ser N-1 porque o sino el impulso inicial no estaría centrado
        if ((0.4 * (N - 1) < i) && (i < 0.6 * (N - 1)) && (0.4 * (N - 1) < j) && (j < 0.6 * (N - 1)))
        {

            grid_0[N * i + j] = 20;
        }
        else
        {
            grid_0[N * i + j] = 0;
        }
        break;
    case 1:
        grid_0[N * i + j] = grid_1[N * i + j] + (0.00125 * (grid_1[N * (i+1) + j] + grid_1[N * (i-1) + j] + grid_1[N * i + j-1] + grid_1[N * i + j+1] - 4*grid_1[N * i + j]));
        break;
    default:
        grid_0[N * i + j] = (2 * grid_1[N * i + j]) - grid_2[N * i + j] + (0.0025 * (grid_1[N * (i + 1) + j] + grid_1[N * (i - 1) + j] + grid_1[N * i + j - 1] + grid_1[N * i + j + 1] - 4 * grid_1[N * i + j]));
        break;
    }
}

/**
 * @brief Calcula la difusión de una onda según la equación de Schroendinger
 *
 * @param grid Matriz de flotantes de NxN
 * @param N Tamaño de la matriz
 * @param T Número de iteración de salida
 */
__global__ void schroedinger(float *grid_0, float *grid_1, float *grid_2, int N, int T, int t)
{
    /*for(int t = 0; t <= T; t++){
        for(int i = 1;i < N-1;i++){
            for(int j = 1;j <N-1;j++){
                schroedinger_(grid,N,i,j,t);
            }
        }
        if(t != T){
            for (int i = 1; i < N - 1; i++)
            {
                for (int j = 1; j < N - 1; j++)
                {
                    grid[3*(N*i +j)+ 2] = grid[3*(N*i +j)+ 1];
                    grid[3*(N*i +j)+ 1] = grid[3*(N*i +j)+ 0];
                }
            }
        }
        
    }*/
    int blocksize = blockDim.y * blockDim.x; // number of threads in a TB
    //printf("tamano bloque: %d\n", blocksize);
    int blockId = gridDim.x * blockIdx.y + blockIdx.x; // unique block Id
    
    //int tid = blockId * blocksize + blockDim.x*threadIdx.y + threadIdx.x; // global tid
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    int j = blockDim.y*blockIdx.y + threadIdx.y;
    if (i < N && j < N)
    {

        schroedinger_(grid_0, grid_1, grid_2, N, i, j, t);
        //printf("t: %d tIdX: %d tIdY %d\n",t,tIdX,tIdY);
        //printf("blockId: %d\n", blockId);
        // if (t != T)
        // {
        //     grid_2[N * i + j] = grid_1[N * i + j];
        //     grid_1[N * i + j] = grid_0[N * i + j];
        // }
    }
        
}

__global__ void copyValues(float *grid_0, float *grid_1, float *grid_2, int N, int T, int t){
    int blocksize = blockDim.y * blockDim.x; // number of threads in a TB
    // printf("tamano bloque: %d\n", blocksize);
    int blockId = gridDim.x * blockIdx.y + blockIdx.x; // unique block Id

    // int tid = blockId * blocksize + blockDim.x*threadIdx.y + threadIdx.x; // global tid
    int tIdX = blockDim.x * blockIdx.x + threadIdx.x;
    int tIdY = blockDim.y * blockIdx.y + threadIdx.y;
    if (tIdX < N && tIdY < N)
    {
        if (t != T)
        {
            grid_2[N * tIdX + tIdY] = grid_1[N * tIdX + tIdY];
            grid_1[N * tIdX + tIdY] = grid_0[N * tIdX + tIdY];
        }
    }
}


/**
 * @brief Imprime la matriz, junto con los valores de las iteraciones
 * anteriores
 *
 * @param grid Matriz de flotantes de NxN
 * @param N Tamaño de la matriz
 */
void printGrid(float* grid,int N){
    for(int i = 0;i <N;i++){
        for(int j = 0;j <N;j++){
            printf("%f ",grid[N*i +j]);
        }
        printf("\n");
    }
}

/**
 * @brief Guarda el valor de la matriz en un archivo en formato binario
 *
 * @param grid Matriz de flotantes de NxN
 * @param N Tamaño de la matriz
 * @param filename Nombre del archivo de salida
 */
void save(float* grid,int N,char*filename){
    FILE *f = fopen(filename,"w");
    fwrite(grid,sizeof(float),N*N,f);
    fclose(f);
}

__host__ int main(int argc, char *const *argv)
{

    options_t options;
    options.filename = (char *)malloc(15 * sizeof(char));

    if (argc != 11)
    {
        printf("Argumentos Incorrectos\nUSO: ./wave -N tamano_grilla -x tamano_bloque_en_X -y tamano_bloque_en_Y -T numero_de_pasos -f archivo_de_salida\n");
        exit(1);
    }
    int c;

    while (((c = getopt(argc, argv, "N:x:y:T:f:")) != -1))
    {
        switch (c)
        {
        case 'N':
            options.grid_size = atof(optarg);
            break;
        case 'x':
            options.x = atof(optarg);
            break;
        case 'y':
            options.y = atof(optarg);
            break;
        case 'T':
            options.steps = atof(optarg);
            break;
        case 'f':
            options.filename = optarg;
            break;
        default:
            break;
        }
    }

    dim3 blocksize;
    dim3 gridsize;

    gridsize.x = options.grid_size / blocksize.x;
    gridsize.y = options.grid_size / blocksize.y;
    blocksize.x = options.x;
    blocksize.y = options.y;
    
    //printf("N:%d T:%d H:%d f:%s\n", options.grid_size, options.steps, options.threads, options.filename, options.steps);
    float *grid_0 = initializeGrid(options.grid_size);
    float *grid_1 = initializeGrid(options.grid_size);
    float *grid_2 = initializeGrid(options.grid_size);

    float *d_grid_0;
    float *d_grid_1;
    float *d_grid_2;
    cudaMalloc((void **) &d_grid_0, options.grid_size*options.grid_size*sizeof(float));
    cudaMemcpy(d_grid_0, grid_0, options.grid_size * options.grid_size * sizeof(float), cudaMemcpyHostToDevice);

    cudaMalloc((void **) &d_grid_1, options.grid_size * options.grid_size * sizeof(float));
    cudaMemcpy(d_grid_1, grid_1,options.grid_size * options.grid_size * sizeof(float), cudaMemcpyHostToDevice);

    cudaMalloc((void **) &d_grid_2, options.grid_size * options.grid_size * sizeof(float));
    cudaMemcpy(d_grid_2, grid_2,options.grid_size * options.grid_size * sizeof(float), cudaMemcpyHostToDevice);

    for(int t = 0; t <= options.steps; t++){
        schroedinger<<<gridsize,blocksize>>>(d_grid_0, d_grid_1, d_grid_2, options.grid_size, options.steps, t);
        //cudaDeviceSynchronize();
        copyValues<<<gridsize, blocksize>>>(d_grid_0, d_grid_1, d_grid_2, options.grid_size, options.steps, t);
        cudaDeviceSynchronize();
    }
    cudaMemcpy(grid_0, d_grid_0, options.grid_size * options.grid_size * sizeof(float), cudaMemcpyDeviceToHost);
    //printGrid(grid, options.grid_size);
    
    cudaFree(d_grid_0);
    cudaFree(d_grid_1);
    cudaFree(d_grid_2);
    save(grid_0, options.grid_size, options.filename);
    free(grid_0);
    free(grid_1);
    free(grid_2);
}