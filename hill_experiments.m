%% Problem 1 Part 1: Fig 7A and 7B Simulations
close all;clear;clc;format short;format compact;

figure(1);

t = [linspace(0,5,1000)']; % set fixed time axis
l = [repmat(1,1000,1)]; % the first is isometric contraction 'A'
[p,h]=hill(l,t);

% Simulating Fig 7A
ctr=1; % array position counter
dec=-0.05;
d = 0.85:dec:0.3; % shortening length decrements

start=550; % shortening start time
short=100; % shortening duration
send=1000-start-short; % shortening end time

ex_hs=zeros(1,numel(d)); % preallocate

for d=d
    l = [repmat(1,start,1);linspace(1,d,short)';repmat(d,send,1)]; % different shortening distances (d)
    [pd,hd]=hill(l,t);
    ex_hs(ctr) = hd(end)-h(end); % excess heat = h of deflection - h of isometric contraction
    subplot(2,1,1);plot(t,h);plot(t,hd); title('Fig 7A Simulation'); 
    axis([0 5 0 60]); xlabel('Time (sec)'); ylabel('Heat'); grid on; hold on
    ctr=ctr+1;
end

% Simulating Fig 7B
ctr=1;
dinc=21;
dur=100:dinc:350; % duration increases
cdist=0.55; % constant distance to which we shorten

for dur=dur % vary durations for simulating shortening different loads (l)
    res=start+dur;  % time at which load change ends and fixed load resumes 
    l = [repmat(1,start,1);linspace(1,cdist,dur)';repmat(cdist,res,1)];
    [pl,hl]=hill(l,t);
    subplot(2,1,2);plot(t,h); plot(t,hl); title('Fig 7B Simulation'); 
    axis([0 5 0 60]); xlabel('Time (sec)'); ylabel('Heat'); grid on; hold on

end

%% Visualize load response of the hill.p system for 7A and 7B
figure(2);

after=start+short; % timepoint immediately after end of shortening

ctr=1;
d=0.85:dec:0.3;
pds=zeros(1,numel(d)); % preallocate

for d=d
    l = [repmat(1,start,1);linspace(1,d,short)';repmat(d,send,1)]; % the varied shortening lengths (d)
    [pd,hd]=hill(l,t);
    pds(ctr) = pd(after);
    subplot(2,1,1); plot(t,p);plot(t,pd); title('Load responses from Fig 7A Simulation'); 
    axis([0 5 -20 200]); xlabel('Time (sec)'); ylabel('Load'); grid on; hold on
    ctr=ctr+1;
end

dur=100:dinc:350; % duration increases
ctr=1;
for dur=dur % vary durations for simulating shortening different loads (l)
    res=start+dur;  % time at which load change ends and fixed load resumes 
    l = [repmat(1,start,1);linspace(1,cdist,dur)';repmat(cdist,res,1)];
    [pl,hl]=hill(l,t);
    pls(ctr)=pl(res);
    subplot(2,1,2);plot(t,p); plot(t,pl); title('Load Responses from Fig 7B Simulation'); 
    axis([0 5 0 200]); xlabel('Time (sec)'); ylabel('Load'); grid on; hold on
    ctr=ctr+1;
end

%% Problem 1 Part 2: Hill constants 'a' and 'b'
ctr=1;
for d=0.85:dec:0.3 % distances over which we shortened
    dists(ctr) = 1-d;
    ctr=ctr+1;
end

% Hill constant 'a' = excess heat / distance shortened

a = ex_hs/dists; % calculated directly from excess heat vs distance shortened

figure(3);
subplot(2,1,1); plot(dists, ex_hs, 'ro', 'MarkerSize', 8, 'LineWidth', 2); grid on;
xlabel('Distances Shortened');
ylabel('Excess Heat');
title('Linear Fit of Hill Constant a');

% Linear fit
a_fit = polyfit(dists, ex_hs, 1);
xFit = linspace(0, 1, 15);
yFit = polyval(a_fit, xFit);
hold on;

plot(xFit, yFit, 'b', 'MarkerSize', 15, 'LineWidth', 1);
legend('Simulated data', 'Fit', 'Location', 'Northwest');

%-----------------------------------------------------------%
% Hill constant 'b' = rate of excess energy liberated / load

% Calculate velocities from figure 7A
ctr=1;
for d=0.85:dec:0.3 % distances shortened
    v(ctr) = (1-d)/100; % 1-d = the distance over which we shortened, 100 = 0.5s (constant shortening time)
    ctr=ctr+1;
end

v_sec = v .* 200; % convert to seconds

ddtexe = (pds+a) .* v_sec; % calculating rate of excess energy liberated
p0 = p(end); % p0 is the load of the isometric contraction 'E'
load=p0-pds;

b = ddtexe/load; % calculated directly

subplot(2,1,2); plot(load, -ddtexe, 'ro', 'MarkerSize', 8, 'LineWidth', 2); grid on;
xlabel('Loads');
ylabel('d/dt Excess Energy Liberated');
title('Linear Fit of Hill Constant b');

% Linear fit
b_fit = polyfit(load, -ddtexe, 1); % plotting rate as negative to graphically represent a negative slope
xFit = linspace(60, 140, 15);
yFit = polyval(b_fit, xFit);
hold on;

plot(xFit, yFit, 'b', 'MarkerSize', 15, 'LineWidth', 1);
legend('Simulated data', 'Fit', 'Location', 'Northeast');
%% Problem 1 Part 3: Calculating system properties

p0 % p0 is the load of the isometric contraction 'E' 

vmax_graph = max(v_sec) % vmax calculated graphically from 7A data

vmax_calc = p0*b / a % derived from the force-velocity hyperbolic relationship
                     % units muscle length/s

%% Problem 2 Part 1: Simulating constant-velocity releases
%  Force and Velocity data as in Fig 12, Hill

% Attempt 1: Fig 7A data
pds; % 7A forces
v_sec; % 7A velocities

% Attempt 2: 7B
pls; % 7B forces

ctr=1;
dur=100:dinc:350;
for dur=dur % 7A velocities
    v_f(ctr) = (1-cdist)/dur; % 1-cdist = the constant distance over which we shortened
    ctr=ctr+1;
end

v_f = v_f *200; % convert to seconds
