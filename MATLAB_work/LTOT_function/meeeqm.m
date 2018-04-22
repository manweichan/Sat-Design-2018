function ydot = meeeqm(t, y)

% first-order modified equinoctial equations of motion

% input

%  t    = current simulation time (seconds)
%  y(1) = semilatus rectum of orbit (kilometers)
%  y(2) = f equinoctial element
%  y(3) = g equinoctial element
%  y(4) = h equinoctial element
%  y(5) = k equinoctial element
%  y(6) = true longitude (radians)

% output

%  ydot(1:6) = first time derivatives of modified equinoctial elements

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mu v1 beta0 thracc

% unload current modified equinoctial orbital elements

pmee = y(1);
fmee = y(2);
gmee = y(3);
hmee = y(4);
xkmee = y(5);
xlmee = y(6);

% compute current absolute value of yaw angle (radians)

beta_tmp = atan3(v1 * sin(beta0), (v1 * cos(beta0) - thracc * t));

% compute current eci state vector

[reci, veci] = mee2eci(mu, y);

% compute current classical orbital elements

oev_wrk = eci2orb1 (mu, reci, veci);

% current argument of latitude (radians)

arglat = mod(oev_wrk(4) + oev_wrk(6), 2.0 * pi);

% compute yaw angle based on current orbital position (radians)

if (arglat >= 0.0 && arglat < 0.5 * pi)
    
    beta_wrk = -beta_tmp;
    
end

if (arglat > 0.5 * pi && arglat < pi)
    
    beta_wrk = beta_tmp;
    
end

if (arglat > pi && arglat < 1.5 * pi)
    
   beta_wrk = beta_tmp;
   
end

if (arglat > 1.5 * pi && arglat < 2.0 * pi)
    
   beta_wrk = -beta_tmp;
   
end

% compute modified equinoctial elements equations of motion

sinl = sin(xlmee);

cosl = cos(xlmee);

wmee = 1.0 + fmee * cosl + gmee * sinl;

sesqr = 1.0 + hmee * hmee + xkmee * xkmee;

ydot(1) = (2.0 * pmee / wmee) * sqrt(pmee / mu) * thracc * cos(beta_wrk);

ydot(2) = sqrt(pmee / mu) * (((wmee + 1.0) * cosl + fmee) * (thracc * cos(beta_wrk) / wmee) ...
    -(hmee * sinl - xkmee * cosl) * (gmee * thracc * sin(beta_wrk) / wmee));

ydot(3) = sqrt(pmee / mu) * (((wmee + 1.0) * sinl + gmee) * (thracc * cos(beta_wrk) / wmee) ...
    -(hmee * sinl - xkmee * cosl) * (fmee * thracc * sin(beta_wrk) / wmee));

ydot(4) = sqrt(pmee / mu) * (sesqr * thracc * sin(beta_wrk) / (2.0 * wmee)) * cosl;

ydot(5) = sqrt(pmee / mu) * (sesqr * thracc * sin(beta_wrk) / (2.0 * wmee)) * sinl;

ydot(6) = sqrt(mu * pmee) * (wmee / pmee)^2 + (1.0 / wmee) * sqrt(pmee / mu) ...
    * (hmee * sinl - xkmee * cosl) * thracc * sin(beta_wrk);




