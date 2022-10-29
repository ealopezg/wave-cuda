CC = gcc
CFLAGS=-I. -fopenmp
wave:
	$(CC) -o wave wave.c $(CFLAGS)

clean:
	rm -f wave *.raw

test:
	./wave -N 256 -T 10000 -H 12 -f out.raw