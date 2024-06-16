% Thomas Murphy
% 5/7/2023
% Project 2: Hitting a home run, with air resistance
% Phase 3: Exporting data and analyzing it in Excel

clear
clf


% ----- define given information -----

R0 = 446;   % measured HR distance, in ft
v0mph = 112;   % exit velocity, in mph
phi0deg = 32;   % launch angle, in deg

m = 0.1417;    % mass of a baseball, in kg
% from https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/forces-on-a-baseball/


x0 = 0; y0 = 0;  % initial x and y velocitiy
g = 10;   % gravitational constant, N/kg = m/s^2 

% time is in s, distance in m (speed in m/s, etc.)


% ----- set up more variables -----

mph2mps = 5280 * 12 * 2.54 / 100 / 3600;  % mph to m/s conversion
deg2rad = pi()/180;   % degrees to radians
m2ft = 3.281; % m to ft

v0 = v0mph*mph2mps;   % initial speed
phi0 = phi0deg*deg2rad;   % initial angle (in rad)

v0x = v0*cos(phi0);   % x-component of velocity
v0y = v0*sin(phi0);   % y-component of velocity


% ----- compute some useful quantities for the trajectory -----

tH = v0y/g;   % time to reach max. height
tLand = 2*tH;   % time to land (time of flight)

H = tH * v0y/2;   % max. height
R = v0x * tLand;   % range

H_ft = H*m2ft;  % max. height converted to feet
R_ft = R*m2ft;  % range converted to feet


% ---- set up a time array, compute x(t), y(t) analytically -----

tmin = 0; tmax = tLand; % start time, end time, in s
N = 2000;   % intervals

t = linspace(tmin, tmax, 1+N);   % time array, connects x(t) and y(t), in s

xt = x0 + v0x*t;   % analytic, x(t), no drag, ax = 0
yt = y0 + v0y*t - (1/2)*g*t.^2;   % analytic, y(t), no drag, ay = -g


% ----- add numeric solution -----
pair = 1.225; % density of air, in kg/m^3
C = input('Constant: '); % drag coef., between 0.2 and 0.3
A = 0.00426; % cross sectional area of a baseball, in m^2
% found on http://spiff.rit.edu/richmond/baseball/traj/traj.html
k = 0.5*C*pair*A; % constant used in Fdrag

dt = (tmax-tmin)/N; % time interval

y = zeros(1, N+1);   % initialize y(t)

y(1) = y0;   % initialize y(0)
vy = v0y;   % vy(1) = v0y
            
x = zeros(1, N+1);  % initialize x(t)

x(1) = x0;  % initialize x(0)
vx = v0x; % vx(1) = 0

for n = 1:N   % stop at N
    v = (vx^2 + vy^2)^0.5; % speed of the baseball
    
    Fnety = -g*m - k*v*vy;   % net force in y direction, N*m/s^2
    Fnetx = -k*v*vx;  % net force in x direction, N*m/s^2

    ay = Fnety/m;   % y acceleration, m/s^2
    ax = Fnetx/m;   % x acceleration, m/s^2

    y(n+1) = y(n) + vy*dt + (1/2)*ay*dt^2;   % updating height y at time t
    vy = vy + ay*dt;   % updating y velocity at time t
    
    x(n+1) = x(n) + vx*dt + (1/2)*ax*dt^2;  % updating distance x at time t
    vx = vx + ax*dt;    % updating x velocity at time t

end

checkSumy = sum(abs(y-yt))   % compare analytic to numeric, should be 0
checkSumx = sum(abs(x-xt))   % compare analytic to numeric, should be 0

xt = xt*m2ft; yt = yt*m2ft;   % convert analytic to mph
x = x*m2ft; y = y*m2ft;   % convert numeric to mph

plot(xt, yt, x, y, 'LineWidth', 2)  % plotting analytic and numeric functions

grid on; grid minor;

ax = gca; ax.FontSize = 14; ax.GridAlpha = 0.4; ax.MinorGridAlpha = 0.5;
% axis font size, 14 pt

xlabel('x (ft)', 'FontSize', 16)   % xlabel, 16 pt
ylabel('y (ft)', 'FontSize', 16)    % ylabel, 16 pt

title({'ECE 202, Project 2, Phase 3:',...
    'Trajectory of a baseball', ...
    ' with drag vs.  without drag'}, ...
    'FontSize', 18)   % title, 18 pt

legend({'without drag', sprintf('with drag (C = %g)',C)}, 'FontSize', 16)
% legend, 16 pt

ylim([-2 130])   % add a little space on the bottom, more on top for legend


% ----- make labels, export data -----

export = [t; x; y].'; % x and y in columns

labels = ["time t (s)", "x (ft)", "y (ft)"]; % labels of columns

export = [labels;export];
% putting the labels at the front, adding a heading

writematrix(export, 'Proj2Phase3.csv')
