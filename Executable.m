clear all
clc
close all
%%
tic
xo=[
    0.70, 0.20, 0.05, 0.05;
    0.60, 0.10, 0.10, 0.20;
    0.40, 0.40, 0.10, 0.10;
    0.40, 0.05, 0.40, 0.05;    
    0.30, 0.30, 0.10, 0.30;
    0.20, 0.20, 0.20, 0.40;
    0.10, 0.80, 0.05, 0.05;
    0.10, 0.60, 0.10, 0.20;
    0.10, 0.40, 0.20, 0.30;
    0.10, 0.10, 0.60, 0.20;
    0.10, 0.10, 0.20, 0.60;
    0.10, 0.05, 0.80, 0.05;
    0.10, 0.05, 0.05, 0.80;
    0.05, 0.05, 0.80, 0.10;
    0.01,0.01,0.01,0.97;
    ];

L=length(xo);

for i=1:L
x10=xo(i,1);
x20=xo(i,2);
x30=xo(i,3);
l=100;
cr_NRTL(x10,x20,x30,l)
end
toc

