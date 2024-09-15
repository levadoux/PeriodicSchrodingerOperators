% main2D.m - eigenvalues of periodic Schrödinger operator in 1D
% The period of the potential V is equal to 2pi along both axes. 
% The used parameters to solve problems on [0,2pi]*[0,2pi]
% with generalized periodic BCs are :
% 1) a random draw of n_rand values in the Brillouin zone [0,1]*[0,1]
% --> the founded eigenvalues are plotted in BLUE
% 2) a discretization (n_lines values for each edge)
% of the boundaries of the square [0,0.5]*[0,0.5]
% (a quarter of the Brillouin zone) 
% --> the founded eigenvalues are plotted in BLACK


clear all;

N=10;
Neig=12;
h = 2*pi/N; x = h*(1:N); y = x;

n_rand = 1000; n_lines = 400;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONSTRUCTION OF THE DIFFERENTIATION MATRICES 
% --> see Trefethen's book: "Spectral Method in MATLAB"

% first order
column1 = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';
D1 = toeplitz(column1, -column1); % antisymetric matrix

% second order
column2 = [-pi^2/(3*h^2)-1/6 ...
    -.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
D2 = toeplitz(column2);  % 2nd-order differentiation
% symetric matrix

% !!! 2D CASE: MATRICES OF SIZE N^2 (KRONECKER PRODUCT) !!!

D1X = kron(eye(N),D1);
D1Y = kron(D1,eye(N));
D2X = kron(eye(N),D2);
D2Y = kron(D2,eye(N));
I = eye(N^2);

[xx,yy] = meshgrid(x,y);
xd=xx(:);
yd=yy(:);
potential = diag(pot2D(xd,yd));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOM LOOP

spect_rand=[]; spect_square = [];

disp("Random loop...")
disp("====================")

kvalues = rand(2,n_rand);

tic
it = 1;

for it=1:n_rand

    k1=kvalues(1,it);
    k2=kvalues(2,it);

    H = -(D2X+2*1i*k1*D1X-(k1^2)*I)-(D2Y+2*1i*k2*D1Y-(k2^2)*I);
    H = H + potential;

    v = eigs(H,Neig,'smallestreal');
    spect_rand = [spect_rand,v];

    if mod(it, n_rand/20)==0
        fprintf("=")
    end
end

fprintf(newline)
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SQUARE LOOP


disp("Square loop...")
disp("====================")

% Vertices of the square
A = 0; B = 0.5; C = 0.5+0.5*1i; D = 0.5*1i;
edge1 = linspace(A,B,n_lines);
edge2 = linspace(B,C,n_lines);
edge3 = linspace(C,D,n_lines);
edge4 = linspace(D,A,n_lines);

kvalues = [edge1,edge2,edge3,edge4];

tic
it = 1;

for k=kvalues

    k1=real(k);
    k2=imag(k);

    H = -(D2X+2*1i*k1*D1X-(k1^2)*I)-(D2Y+2*1i*k2*D1Y-(k2^2)*I);
    H = H + potential;

    v = eigs(H,Neig,'smallestreal');
    spect_square = [spect_square,v];

    if mod(it, 4*n_lines/20)==0
        fprintf("=")
    end
    it = it+1;
end

fprintf(newline)
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()

scatter(real(spect_rand), imag(spect_rand),2, 'filled','b'); hold on;
scatter(real(spect_square), imag(spect_square),3, 'filled','k'); hold on;

spect = [spect_rand,spect_square];

xm=min(real(spect),[],'all');
xM=max(real(spect),[],'all');
xl=xM-xm;
xlim([xm-0.1*xl xM+0.1*xl]);
ym=min(imag(spect),[],'all');
yM=max(imag(spect),[],'all');
yl=yM-ym;


if yl==0 % Used if the potential is real
    yl=1;
    ax = gca; % Get current axis handle
    ax.YTick = -1:1:1;
end

ylim([ym-0.1*yl yM+0.1*yl]);
title("Spectrum, "+Neig+" eigenvalues")


