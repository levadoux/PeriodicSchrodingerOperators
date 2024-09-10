function Z = pot1D(x)


% Real potential
% Z = sqrt(x);
% Z = cos(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Complex potential
% Z = 1i*sqrt(x);
% Z = 1i*cos(x);

% Very interesting: Mathieu potentials with mu on the unit circle !
mu = exp(1i*5);
Z = mu*cos(x);





end