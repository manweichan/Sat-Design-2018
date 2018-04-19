function [r, v] = mee2eci(mu, mee)

% convert modified equinoctial orbital
% elements to eci position and velocity vectors

% input

%  mu     = gravitational constant (km**3/sec**2)
%  mee(1) = semilatus rectum of orbit (kilometers)
%  mee(2) = f equinoctial element
%  mee(3) = g equinoctial element
%  mee(4) = h equinoctial element
%  mee(5) = k equinoctial element
%  mee(6) = true longitude (radians)

% output

%  r = eci position vector (kilometers)
%  v = eci velocity vector (kilometers/second)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% unload equinoctial orbital elements

pmee = mee(1);
fmee = mee(2);
gmee = mee(3);
hmee = mee(4);
kmee = mee(5);
lmee = mee(6);

smovrp = sqrt(mu / pmee);

tani2s = hmee^2 + kmee^2;

cosl = cos(lmee);

sinl = sin(lmee);

wmee = 1 + fmee * cosl + gmee * sinl;

radius = pmee / wmee;

hsmks = hmee^2 - kmee^2;

ssqrd = 1 + tani2s;

% compute eci position vector

r(1) = radius * (cosl + hsmks * cosl + 2 * hmee * kmee * sinl) / ssqrd;

r(2) = radius * (sinl - hsmks * sinl + 2 * hmee * kmee * cosl) / ssqrd;

r(3) = 2 * radius * (hmee * sinl - kmee * cosl) / ssqrd;

% compute eci velocity vector

v(1) = - smovrp * (sinl + hsmks * sinl - 2 * hmee * kmee * cosl + gmee ...
       - 2 * fmee * hmee * kmee + hsmks * gmee) / ssqrd;

v(2) = - smovrp * (-cosl + hsmks * cosl + 2 * hmee * kmee * sinl - fmee ...
       + 2 * gmee * hmee * kmee + hsmks * fmee) / ssqrd;

v(3) = 2 * smovrp * (hmee * cosl + kmee * sinl + fmee * hmee ...
       + gmee * kmee) / ssqrd;


