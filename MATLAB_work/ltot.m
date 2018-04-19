% ltot.m            October 22, 2013

% low-thrust orbit transfer between non-coplanar circular orbits

% integrated solution with Edelbaum yaw steering

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

global req mu dtr rtd rkcoef beta0

global v1 tdflag thracc

% initialize rkf78 method

rkcoef = 1;

% read astrodynamic constants and conversion factors

om_constants;

% request inputs

clc; home;

fprintf('\nLow-thrust Orbit Transfer Between Non-coplanar Circular Orbits\n');

while (1)

    fprintf('\n\nplease input the initial altitude (kilometers)\n');

    alt1 = input('? ');

    if (alt1 > 0.0)
        break;
    end

end

while (1)

    fprintf('\nplease input the final altitude (kilometers)\n');

    alt2 = input('? ');

    if (alt2 > 0.0)
        break;
    end

end

while (1)

    fprintf('\nplease input the initial orbital inclination (degrees)');
    fprintf('\n(0 <= inclination <= 180)\n');

    inc1 = input('? ');

    if (inc1 >= 0.0 && inc1 <= 180.0)
        break;
    end

end

while (1)

    fprintf('\nplease input the final orbital inclination (degrees)');
    fprintf('\n(0 <= inclination <= 180)\n');

    inc2 = input('? ');

    if (inc2 >= 0.0 && inc2 <= 180.0)
        break;
    end

end

while (1)

    fprintf('\nplease input the initial RAAN (degrees)');
    fprintf('\n(0 <= RAAN <= 360)\n');

    raan0 = input('? ');

    if (raan0 >= 0.0 && raan0 <= 360.0)
        break;
    end

end

while (1)

    fprintf('\nplease input the thrust acceleration (meters/second^2)\n');

    tacc = input('? ');

    if (tacc > 0)
        break;
    end

end

% convert thrust acceleration to km/sec^2

thracc = tacc / 1000.0;

% convert inclinations to radians

inc1 = inc1 * dtr;

inc2 = inc2 * dtr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the orbit transfer problem %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate total inclination change

dinct = abs(inc2 - inc1);

% check for coplanar orbits

if (dinct == 0)
    
    dinct = 1.0e-8;
    
end

% compute geocentric radii of initial and final orbits (kilometers)

r1 = req + alt1;

r2 = req + alt2;

% compute "local circular velocity" of initial and final orbits

v1 = sqrt(mu / r1);

v2 = sqrt(mu / r2);

% initial yaw angle

beta0 = atan3(sin(0.5 * pi * dinct), (v1/v2) - cos(0.5 * pi * dinct));

% delta-v

dvt = v1 * cos(beta0) - v1 * sin(beta0) / tan(0.5 * pi * dinct + beta0);

% thrust duration

tdur = dvt / thracc;

thrdur_sec = tdur;

if (tdur < 3600.0)
    
    % minutes
    
    tdflag = 1;
    
    tdur = tdur / 60.0;
    
elseif (tdur < 86400.0)
    
    % hours
    
    tdflag = 2;
    
    tdur = tdur / 3600.0;
    
else
    
    % days
    
    tdflag = 3;
    
    tdur = tdur / 86400.0;
    
end

dtstep = tdur / 100.0;

tsim = -dtstep;

for i = 1:1:101
    
    tsim = tsim + dtstep;
    
    if (tdflag == 1)
        
        tsec = 60.0 * tsim;
        
    elseif (tdflag == 2)
        
        tsec = 3600.0 * tsim;
        
    else
        
        tsec = 86400.0 * tsim;
        
    end
    
    t(i) = tsim;
    
    beta(i) = rtd * atan3(v1 * sin(beta0), (v1 * cos(beta0) ...
        - thracc * tsec));
    
    tmp1 = atan((thracc * tsec - v1 * cos(beta0))/(v1 * sin(beta0)));
    
    dinc = rtd * (2/pi) * (tmp1 + 0.5 * pi - beta0);
    
    inc(i) = rtd * inc1 - dinc;
    
    v(i) = 1000.0 * sqrt(v1 * v1 - 2.0 * v1 * thracc * tsec * cos(beta0) ...
        + thracc * thracc * tsec * tsec);
    
    sma(i) = 1.0e6 * mu / (v(i) * v(i));
    
end

% print analytic results

clc; home;

fprintf('\n   Low-thrust Orbit Transfer Analysis \n\n');

fprintf('initial orbit altitude      %10.4f kilometers \n\n', alt1);

fprintf('initial orbit inclination   %10.4f degrees \n\n', inc1 * rtd);

fprintf('initial RAAN                %10.4f degrees \n\n', raan0);

fprintf('initial orbit velocity      %10.4f meters/second \n\n\n', 1000.0 * v1);

fprintf('final orbit altitude        %10.4f kilometers \n\n', alt2);

fprintf('final orbit inclination     %10.4f degrees \n\n', inc2 * rtd);

fprintf('final orbit velocity        %10.4f meters/second \n', 1000.0 * v2);

fprintf('\ntotal inclination change    %10.4f degrees\n\n', rtd * dinct);

fprintf('total delta-v               %10.4f meters/second \n\n', 1000.0 * dvt);

if (tdflag == 1)
    
    fprintf('thrust duration             %10.4f minutes \n\n', tdur);
    
elseif (tdflag == 2)
    
    fprintf('thrust duration             %10.4f hours \n\n', tdur);
    
else
    
    fprintf('thrust duration             %10.4f days \n\n', tdur);
    
end

fprintf('initial yaw angle           %10.4f degrees \n\n', rtd * beta0);

fprintf('thrust acceleration         %10.6f meters/second^2 \n\n', 1000.0 * tacc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create integrated solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial semimajor axis (kilometers)

oev_initial(1) = r1;

% initial orbital eccentricity (nd)

oev_initial(2) = 0.0;

% initial orbital inclination (radians)

oev_initial(3) = inc1;

% initial argument of periapsis (radians)

oev_initial(4) = 0.0 * dtr;

% initial RAAN (radians)

oev_initial(5) = raan0 * dtr;

% initial true anomaly (radians)

oev_initial(6) = 0.0;

% convert user-defined classical orbital elements to eci state vector

[reci_initial, veci_initial] = orb2eci(mu, oev_initial);

% convert eci state vector to modified equinoctial orbital elements

mee_initial = eci2mee(mu, reci_initial, veci_initial);

% data generation step size (seconds)

deltat = 1800.0d0;

% final simulation time (seconds)

tfinal = thrdur_sec;

% number of differential equations

neq = 6;

% truncation error tolerance

tetol = 1.0d-12;

tf = 0.0d0;

nplot = 1;

[reci, veci] = mee2eci(mu, mee_initial(1:6));

oev_wrk = eci2orb1 (mu, reci, veci);

%%%%%%%%%%%%%%%%%%%%%%%
% initial graphics data
%%%%%%%%%%%%%%%%%%%%%%%

% simulation time (seconds)

xplot(nplot) = 0.0;

% x position coordinate (kilometers)

yplot1(nplot) = reci(1) / req;

% y position coordinate (kilometers)

yplot2(nplot) = reci(2) / req;

% z position component (kilometers)

yplot3(nplot) = reci(3) / req;

% yaw angle (degrees)

yplot4(nplot) = rtd * beta0;

% semimajor axis (kilometers)

yplot5(nplot) = oev_wrk(1);

% orbital eccentricity (nd)

yplot6(nplot) = oev_wrk(2);

% orbital inclination (degrees)

yplot7(nplot) = rtd * oev_wrk(3);

% RAAN (degrees)

yplot8(nplot) = rtd * oev_wrk(5);

% integrate equations of motion and create data arrays

while (1)
    
    h = 10.0d0;
    
    ti = tf;
    
    tf = ti + deltat;
    
    if (tf > tfinal)
        
        tf = tfinal;
        
    end
    
    mee_final = rkf78 ('meeeqm', neq, ti, tf, h, tetol, mee_initial);
    
    [reci, veci] = mee2eci(mu, mee_final(1:6));
    
    oev_wrk = eci2orb1 (mu, reci, veci);
    
    pmee = mee_final(1);
    fmee = mee_final(2);
    gmee = mee_final(3);
    hmee = mee_final(4);
    xkmee = mee_final(5);
    xlmee = mee_final(6);
    
    % current yaw angle (radians)
    
    beta_tmp = atan3(v1 * sin(beta0), (v1 * cos(beta0) - thracc * tf));
    
    % current argument of latitude (radians)
    
    arglat = mod(oev_wrk(4) + oev_wrk(6), 2.0 * pi);
    
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
    
    % load data arrays for plotting
    
    nplot = nplot + 1;
    
    % simulation time (seconds)
    
    xplot(nplot) = ti;
    
    % x position coordinate (Earth radii)
    
    yplot1(nplot) = reci(1) / req;
    
    % y position coordinate (Earth radii)
    
    yplot2(nplot) = reci(2) / req;
    
    % z position component (Earth radii)
    
    yplot3(nplot) = reci(3) / req;
    
    % yaw angle (degrees)
    
    yplot4(nplot) = rtd * beta_wrk;
    
    % semimajor axis (kilometers)
    
    yplot5(nplot) = oev_wrk(1);
    
    % orbital eccentricity (nd)
    
    yplot6(nplot) = oev_wrk(2);
    
    % orbital inclination (degrees)
    
    yplot7(nplot) = rtd * oev_wrk(3);

    % RAAN (degrees)
    
    yplot8(nplot) = rtd * oev_wrk(5);
    
    if (tf == tfinal)
        
        break
        
    end
    
    % reload integration vector
    
    mee_initial = mee_final;
    
end

% print results

fprintf('\nfinal classical orbital elements and state vector - integrated solution');
fprintf('\n-----------------------------------------------------------------------\n');

oeprint1(mu, oev_wrk, 2);

[reci_final, veci_final] = orb2eci(mu, oev_wrk);

svprint(reci_final, veci_final);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display analytic solution graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);

plot(t, beta, 'LineWidth', 1.5);

title('Low-thrust Orbit Transfer - analytic solution', 'FontSize', 16);

if (tdflag == 1)
    
    xlabel('time since ignition (minutes)', 'FontSize', 12);
    
elseif (tdflag == 2)
    
    xlabel('time since ignition (hours)', 'FontSize', 12);
    
else
    
    xlabel('time since ignition (days)', 'FontSize', 12);
    
end

ylabel('|yaw angle| (degrees)', 'FontSize', 12);

grid;

print -depsc -tiff -r300 ltot1.eps

figure(2);

plot(t, inc, 'LineWidth', 1.5);

title('Low-thrust Orbit Transfer - analytic solution', 'FontSize', 16);

if (tdflag == 1)
    
    xlabel('time since ignition (minutes)', 'FontSize', 12);
    
elseif (tdflag == 2)
    
    xlabel('time since ignition (hours)', 'FontSize', 12);
    
else
    
    xlabel('time since ignition (days)', 'FontSize', 12);
    
end

ylabel('orbital inclination (degrees)', 'FontSize', 12);

grid;

print -depsc -tiff -r300 ltot2.eps

figure(3);

plot(t, v, 'LineWidth', 1.5);

title('Low-thrust Orbit Transfer - analytic solution', 'FontSize', 16);

if (tdflag == 1)
    
    xlabel('time since ignition (minutes)', 'FontSize', 12);
    
elseif (tdflag == 2)
    
    xlabel('time since ignition (hours)', 'FontSize', 12);
    
else
    
    xlabel('time since ignition (days)', 'FontSize', 12);
    
end

ylabel('velocity (meters/second)', 'FontSize', 12);

grid;

print -depsc -tiff -r300 ltot3.eps

figure(4);

plot(t, sma, 'LineWidth', 1.5);

title('Low-thrust Orbit Transfer - analytic solution', 'FontSize', 16);

if (tdflag == 1)
    
    xlabel('time since ignition (minutes)', 'FontSize', 12);
    
elseif (tdflag == 2)
    
    xlabel('time since ignition (hours)', 'FontSize', 12);
    
else
    
    xlabel('time since ignition (days)', 'FontSize', 12);
    
end

ylabel('semimajor axis (kilometers)', 'FontSize', 12);

grid;

print -depsc -tiff -r300 ltot4.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display integrated solution graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5);

plot(xplot / 86400.0, yplot4, 'LineWidth', 1.5);

title('Low-thrust Orbit Transfer - integrated solution', 'FontSize', 16);

xlabel('time since ignition (days)', 'FontSize', 12);

ylabel('yaw angle (degrees)', 'FontSize', 12);

grid;

print -depsc -tiff -r300 ltot5.eps;

figure(6);

plot(xplot / 86400.0, yplot5, 'LineWidth', 1.5);

title('Low-thrust Orbit Transfer - integrated solution', 'FontSize', 16);

xlabel('time since ignition (days)', 'FontSize', 12);

ylabel('semimajor axis (kilometers)', 'FontSize', 12);

grid;

print -depsc -tiff -r300 ltot6.eps;

figure(7);

plot(xplot / 86400.0, yplot6, 'LineWidth', 1.5);

title('Low-thrust Orbit Transfer - integrated solution', 'FontSize', 16);

xlabel('time since ignition (days)', 'FontSize', 12);

ylabel('orbital eccentricity', 'FontSize', 12);

grid;

print -depsc -tiff -r300 ltot7.eps;

figure(8);

plot(xplot / 86400.0, yplot7, 'LineWidth', 1.5);

title('Low-thrust Orbit Transfer - integrated solution', 'FontSize', 16);

xlabel('time since ignition (days)', 'FontSize', 12);

ylabel('orbital inclination (degrees)', 'FontSize', 12);

grid;

print -depsc -tiff -r300 ltot8.eps;

figure(9);

plot(xplot / 86400.0, yplot8, 'LineWidth', 1.5);

title('Low-thrust Orbit Transfer - integrated solution', 'FontSize', 16);

xlabel('time since ignition (days)', 'FontSize', 12);

ylabel('RAAN (degrees)', 'FontSize', 12);

grid;

print -depsc -tiff -r300 ltot9.eps;

figure(10);

% create axes vectors

xaxisx = [1 1.5];
xaxisy = [0 0];
xaxisz = [0 0];

yaxisx = [0 0];
yaxisy = [1 1.5];
yaxisz = [0 0];

zaxisx = [0 0];
zaxisy = [0 0];
zaxisz = [1 1.5];

hold on;

grid on;

% plot earth

[x y z] = sphere(24);

h = surf(x, y, z);

colormap([127/255 1 222/255]);

set (h, 'edgecolor', [1 1 1]);

% plot coordinate system axes

plot3(xaxisx, xaxisy, xaxisz, '-g', 'LineWidth', 1);

plot3(yaxisx, yaxisy, yaxisz, '-r', 'LineWidth', 1);

plot3(zaxisx, zaxisy, zaxisz, '-b', 'LineWidth', 1);

% plot transfer orbit

plot3(yplot1, yplot2, yplot3, '-b', 'LineWidth', 1);

xlabel('X coordinate (ER)', 'FontSize', 12);

ylabel('Y coordinate (ER)', 'FontSize', 12);

zlabel('Z coordinate (ER)', 'FontSize', 12);

title('Low-thrust Orbit Transfer', 'FontSize', 16);

axis equal;

view(50, 20);

rotate3d on;

print -depsc -tiff -r300 ltot10.eps;

