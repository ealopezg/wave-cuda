# Labororatorio 1 HPC 2022-2 (OpenMP)
Laboratorio 1 de la asignatura High Performance Computing 2022-2

El objetivo de este laboratorio es implementar un simulador paralelo de la difusión de una onda según la equación de Schroendinger, usando OpenMP como tecnología de paralelización.

## Autores
* Esteban Cruces (esteban.cruces@usach.cl)
* Esteban López (esteban.lopez.g@usach.cl)

## Compilación
```
make
```
## Ejecución
```
./wave -N tamano_grilla -T numero_de_pasos -H numero_de_hebras -f archivo_de_salida
```
### Ejemplo
```
./wave -N 256 -T 10000  -H 12 -f out.raw
```