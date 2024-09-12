% main1D.m - eigenvalues of periodic Schrödinger operator in 1D
% The period of the potential V is equal to 2pi. 
% The used parameters to solve problems on [0,2pi]
% with generalized periodic BCs is a discretization of [0,0.5] 
% (half of the Brillouin zone)
% To sort the eigenvalues with different colors, use sort_eig = true !

clear all;

N=50; % Size of the differantiation matrices (number of points for the trigonometric intepolant)
Neig = 6; % Number of eigenvalues 
h = 2*pi/N; x = h*(1:N); y = x;

n = 1000; % Number of used parameters in [0,0.5]

spect = [];
kvalues = linspace(0,0.5,n); % Used parameters (discretization of [0,0.5] with n points regularly spaced)

sort_eig = false; % If false, all the eigenvalues will be plotted in black.
                  % If true, a different color for each eigenvalue.


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

% Construction of the potential matrix
potential = diag(pot1D(x));
I = eye(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Main Loop...")
disp("====================")

tic
it = 1;

for k1=kvalues

    H = -(D2+2*1i*k1*D1-(k1^2)*I);
    H = H + potential;
    % H is our periodic Schrödinger operator
    v = eigs(H, Neig,'smallestreal'); % To get the Neig first eigenvalues ranked by smallestreal part
    spect = [spect,v];

    if mod(it, n/20)==0
        fprintf("=")
    end
    it = it+1;

end

fprintf(newline)
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',[100 100 1100 600])

subplot(2,1,1)
xx = linspace(0,2*pi,1000);
plot(xx, real(pot1D(xx)),'LineWidth',1,'Color','r');
hold on
plot(xx,imag(pot1D(xx)),'LineWidth',1,'Color','b');
hold on
scatter(x, real(pot1D(x)),10, 'filled', 'r');
hold on
scatter(x, imag(pot1D(x)),10, 'filled', 'b');
legend("real","imag","Npts = "+N)
title("Potential over a period [0,2pi]")



subplot(2,1,2)
if sort_eig
    for l=1:Neig
        scatter(real(spect(l,:)), imag(spect(l,:)),2, 'filled');
        hold on;
    end
else
scatter(real(spect), imag(spect),2, 'filled', 'k');
end

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
title("Spectrum, dk = 1/"+n+", "+Neig+" eigenvalues")
    
  


