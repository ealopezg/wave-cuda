CC = nvcc
CFLAGS=
wave:
	$(CC) -o wave wave.cu $(CFLAGS)

clean:
	rm -f wave *.raw

test:
	./wave -N 256 -T 10000 -x 16 -y 16 -f out.raw