% Modified from seismic-live python version (web) and wave1d.m (Ruhul Firdaus)
% Writer: Yudha Styawan
% email : yudhastyawan26@gmail.com

% opening the session (close all)
clc; clear; close all;

% SWITCHERS
% ========

isAbsorb = 1;               % 1 if Absorbing boundary condition
isplotstf = 0;              % 1 if plot the source time function
isplotstfspec = 0;          % 1 if plot the source spectrum
ishetero = 1;               % 1 if heterogeneous medium (2 layers)
isrecordedspec = 1;         % 1 if show the recorded wave spectrum

% PARAMETER CONFIGURATION
% =======================

% 1. distance
xmin = 10;                  % meter
xmax = 1000;                % meter
dx = 10;                    % meter
x = xmin:dx:xmax;           % xmin dx+xmin 2dx+xmin 3dx+xmin ... xmax
nx = length(x);

% 2. velocity model (edit the model here)
% homogeneous
c0 = 1000.;                 % m/s
c = zeros(1,nx);
c = c + c0;

% heterogeneous (2 layers)
if ishetero == 1
    c1 = 2000;
    c2 = 1000;
    c0 = max(c1,c2);
    xcb = 700;                        % boundary, meter
    icb = round((xcb - xmin)/dx);     % point
    c = zeros(1,nx);
    c(1:icb) = c1;
    c(icb+1:end) = c2;
end

% 3. source position
xsrc = 500;                      % meter
isrc = round((xsrc - xmin)/dx);  % point

% 4. receiver position
xr = 200;                       % meter
ir = round((xr - xmin)/dx);     % point

% 5. time length
tmax = 2;                       % sec
CFL = 1;                        % CFL stability criterion = c0 * dt/dx
dt = CFL*dx/c0;
time = 0:dt:tmax;
nt = length(time);

% SOURCE TIME FUNCTION (GAUSSIAN)
% ===============================

% 1. initialization
f0 = 30.;                                       % max frequency band (see the freq. spectrum)
t0 = 4. / f0;
src = -2. * (time - t0) * (f0 ^ 2) .*...
    (exp(-1.0 * (f0 ^ 2) * (time - t0) .^ 2));  % 1st diff. of gaussian
src = (src(3:end)-src(1:end-2))/2/dt;           % 2nd diff. of gaussian
src = -[0 src 0];
src = src/max(abs(src));                        % normalization

% 2. plot stf
if isplotstf == 1
    figure(1)
    plot(time, src)
    xlabel('t (s)'); ylabel('amplitude'); title('source time function')
end

% 3. plot the frequency spectrum
if isplotstfspec == 1
    Af = abs(fft(src));
    Af = Af(1:round(end/2));
    ff = linspace(0,1/(2*dt),length(Af));
    figure(2)
    plot(ff,Af)
    xlabel('frequency (Hz)'); ylabel('amplitude'); title('source spectrum')
end

%%
% 1D WAVE SIMULATION
% ==================

% 1. initialize empty pressure
p = zeros(1,nx);                % the current p
pold = p;                       % p before the current p
pnew = p;                       % p after the current p
d2px = p;                       % second derivative (without /dx^2)

% 2. initialize empty seismogram
seis = zeros(1,nt);

% -> if recorded wave spectrum is shown -> configure subplot numbers
if isrecordedspec == 1
    nsubplot = 3;
else
    nsubplot = 2;
end

% 3. calculate the finite difference
for it=1:nt
    for i=2:nx-1
        d2px(i) = p(i+1)-2*p(i)+p(i-1);
    end

    % time extrapolation
    k = (c .* dt ./ dx).^2;
    pnew = 2 * p - pold + k .* d2px;

    % Absorbing Boundary Condition (one-way wave equation)
    if isAbsorb == 1
        pnew(1) = p(1) + (c(1)*dt/dx) * ( p(2) - p(1) ); % + menjalar hanya ke kiri
        pnew(nx) = p(nx) - (c(nx)*dt/dx) * ( p(nx) - p(nx - 1) ); % - menjalar hanya ke kanan
    end

    % add source term at isrc
    pnew(isrc) = pnew(isrc) + src(it);

    % remap time levels
    pold = p;
    p = pnew;

    % output seismogram
    seis(it) = p(ir);

    % plot
    idisp = 5;              % set the interval animation (points)
    if mod(it,idisp) == 0
        figure(3)
        
        % animate the propagation
        subplot(nsubplot,1,1)
        plot(x,p)
        hold on
        plot(x(isrc),0,'r+')
        plot(x(ir),0,'bv')
        if ishetero == 1
            xline(x(icb))
        end
        hold off
        ylim([-1,1])
        legend('wave','source position','receiver position')
        xlabel('distance (m)'); ylabel('amplitude'); title('wave propagation')

        % live record
        subplot(nsubplot,1,2)
        plot(time(1:it),seis(1:it))
        ylim([-1,1])
        xlim([0,tmax])
        xlabel('t (s)'); ylabel('amplitude'); title('recorded wave in the station (live)')

        % recorded wave spectrum
        if isrecordedspec == 1
            subplot(3,1,3)
            Sf = abs(fft(seis));
            Sf = Sf(1:round(end/2))/max(abs(Sf));
            fx = linspace(0,1/(2*dt),length(Sf));
            plot(fx,Sf)
            ylim([0,1])
            xlabel('frequency (Hz)'); ylabel('amplitude'); title('recorded wave spectrum')
        end
        
        % delay animation in sec
        pause(0.1)
    end
end