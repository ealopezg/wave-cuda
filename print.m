


N = 256;
colormap("jet");
f = fopen('salida.raw','r');
I = fread(f,N*N,'float32');
I = reshape(I,N,N);
I = I';
%surf(I);
%hold on;
imagesc(I); axis('off');axis('square');