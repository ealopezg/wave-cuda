/**
 * @file wave.c
 * @author Esteban Cruces (esteban.cruces@usach.cl)
 * @author Esteban López (esteban.lopez.g@usach.cl)
 * @brief Laboratorio 1 HPC 2022-2
 * Cálculo de la difusión de una onda según la equación de Schroendinger
 * de forma paralela usando OpenMP
 * @version 1
 * @date 2022-10-09
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <omp.h>

#include "wave.h"




/**
 * @brief Estructura donde se guardan los parámetros de entrada
 * 
 */
typedef struct{
    int grid_size; // Tamaño de la matriz
    int steps; // Cantidad de iteraciones
    int threads; // Número de hebras a utilizar
    char *filename; // Nombre archivo de salida
} options_t;



/**
 * @brief Crea una matriz de 3 dimensiones , N*N*3, donde la tercera dimensión corresponde a los pasos anteriores (0: t,1: t-1,2: t-2)
 *
 * @param N Tamaño de la matriz
 * @return float*** Matriz
 */
float*** initializeGrid(int N){
    float*** grid = (float***)malloc(N*sizeof(float**));
    int i, j, k;
    #pragma omp parallel shared(grid, N) private(i, j, k)
    {
        #pragma omp for
        for (i = 0; i < N; i++)
        {
            grid[i] = (float **)malloc(N * sizeof(float *));
            for (j = 0; j < N; j++)
            {
                grid[i][j] = (float *)malloc(3 * sizeof(float));
                for (k = 0; k < 3; k++)
                {
                    grid[i][j][k] = 0;
                }
            }
        }
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
void schroedinger_(float ***grid, int N, int i, int j, int t)
{

    switch (t)
    {
    case 0:
        // Tiene que ser N-1 porque o sino el impulso inicial no estaría centrado
        if ((0.4 * (N-1) < i) && (i < 0.6 * (N-1)) && (0.4 * (N-1) < j) && (j < 0.6 * (N-1)))
        {
            grid[i][j][0] = 20;
        }
        else
        {
            grid[i][j][0] = 0;
        }
        break;
    case 1:
        grid[i][j][0] = grid[i][j][1] + (0.00125 * (
                grid[i + 1][j][1]
                + grid[i - 1][j][1]
                + grid[i][j - 1][1]
                + grid[i][j + 1][1]
                - (4 * grid[i][j][1])
                )
            );
        break;
    default:
        grid[i][j][0] = (2 * grid[i][j][1]) - grid[i][j][2] + (0.0025 * (
                grid[i + 1][j][1]
                + grid[i - 1][j][1]
                + grid[i][j - 1][1]
                + grid[i][j + 1][1]
                - (4 * grid[i][j][1])
                )
            );
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
void schroedinger(float*** grid, int N, int T){
    int i,j,t;
    for (t = 0; t <= T; t++)
    {
        #pragma omp parallel shared(grid, N, T, t) private(i, j)
        {
            #pragma omp for collapse(2)
            // De 1 a N-1 porque cada valor borde siempre va a ser 0
            for (i = 1; i < N-1; i++)
            {
                for (j = 1; j < N-1; j++)
                {
                    schroedinger_(grid, N, i, j, t);
                }
            }
            #pragma omp barrier
            // Si la iteración t es la final, no es necesario copiar los valores
            if(t != T){
                // Copia los valores de la iteración actual al bloque anterior
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        grid[i][j][2] = grid[i][j][1];
                        grid[i][j][1] = grid[i][j][0];
                    }
                }
            }
            
           
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
void printGrid(float ***grid, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("(%f,%f,%f) ", grid[i][j][0], grid[i][j][1], grid[i][j][2]);
        }
        printf("\n");
    }
}

/**
 * @brief Libera la memoria de la matriz
 *
 * @param grid Matriz de flotantes de NxN
 * @param N Tamaño de la matriz
 */
void freeGrid(float*** grid, int N){
    int i, j;
    #pragma omp parallel shared(grid, N) private(i, j)
    {
        #pragma omp for
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                free(grid[i][j]);
            }
            free(grid[i]);
        }
    }
    
    free(grid);
}

/**
 * @brief Guarda el valor de la matriz en un archivo en formato binario
 *
 * @param grid Matriz de flotantes de NxN
 * @param N Tamaño de la matriz
 * @param filename Nombre del archivo de salida
 */
void save(float*** grid,int N,char*filename){
    float *H = (float*)malloc(N*N*sizeof(float));
    int i, j;
    #pragma omp parallel shared(grid, N) private(i, j)
    {
        #pragma omp for collapse(2)
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                H[i * N + j] = grid[i][j][0];
            }
        }
    }
    
    FILE *f = fopen(filename,"w");
    fwrite(H,sizeof(float),N*N,f);
    fclose(f);
    free(H);
}

int main(int argc, char * const*argv)
{

    options_t options;
    options.filename = (char *)malloc(15*sizeof(char));

    if(argc != 9){
        printf("Argumentos Incorrectos\nUSO: ./wave -N tamano_grilla -T numero_de_pasos -H numero_de_hebras -f archivo_de_salida\n");
        exit(1);
    }
    int c;

    while (( (c = getopt(argc, argv, "N:T:H:f:")) != -1)){
        switch (c)
        {
            case 'N':
                options.grid_size = atof(optarg);
                break;
            case 'T':
                options.steps = atof(optarg);
                break;
            case 'H':
                options.threads = atof(optarg);
                break;
            case 'f':
                options.filename = optarg;
                break;
            default:
                break;
        }
    }
    //printf("N:%d T:%d H:%d f:%s\n",options.grid_size,options.steps,options.threads,options.filename,options.steps);
    omp_set_dynamic(0); // Explicitly disable dynamic teams
    omp_set_num_threads(options.threads); // Sets threads count
    float*** grid = initializeGrid(options.grid_size);
    schroedinger(grid,options.grid_size,options.steps);
    save(grid,options.grid_size,options.filename);
    freeGrid(grid,options.grid_size);
}