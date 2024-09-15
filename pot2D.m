function Z = pot2D(x,y)
% Only interesting to see the lines of the square for complex potentials...
% Use Ctr+R / Ctr+T to add/remove '%'



% Complex potentials - Special case V(x,y) = V(x)+V(y)
% Z = 1i*sqrt(x) + 1i*sqrt(y);
% Z = 1i*sqrt(x) + cos(y);

% Complex potentials - Other cases


Z =  1i*cos(x).*sin(y); % Spectrum only real around Re = 2,5 !!!

% Z = sqrt(x) + 1i*cos(x).*sin(y);  % Spectrum only real around Re = 2,5 !!!
                                    % The lines still give some boundaries !!!

end