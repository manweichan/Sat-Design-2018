function om_constants

% astrodynamic and utility constants

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global dtr rtd mu req j2 egrav

% angular conversion factors

dtr = pi / 180.0;

rtd = 180.0 / pi;

% earth gravitational constant (km**3/sec**2)

mu = 398600.4415;

% earth equatorial radius (kilometers)

req = 6378.1363;

% earth oblateness gravity coefficient (non-dimensional)

j2 = 0.00108263;

% earth surface gravity (meters/second**2)

egrav = 9.80665;


