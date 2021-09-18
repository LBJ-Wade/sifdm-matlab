close all
clear all
clc
%% Philip Mocz (2020), Princeton University
% Merge solitons with SI -- setup similar to Mocz+(2017)
% adds Chavanis equations 33 and 34 ( https://arxiv.org/abs/1710.06268 )

% Internal units:
% [L] = kpc
% [M] = Msun
% [E] = Msun (km/s)^2

% scaling symmetry of the equations
% {x, t, rho, m} --> {a x, b t, b^-2 rho, a^-2 b m}

% Mtot ~ 4e9 Msun
% Msoli,max, initial ~ 1e9 Msun
% goal: form solitn above critical mass at ~ 2e9 Msun / see what happens

% Mmax = 1.012 hbar / sqrt(G m |a_s|)
% f = sqrt(hbar c^2 m / (32 pi |a_s|));
% Mmax = 10.1468 * f/m * sqrt(hbar/(G*c^3))
% 10.1468*10^15 GeV/(10^-22 eV/c^2) * sqrt(hbar/(G*c^3)) / mass of sun = 1.1e9

% v = (hbar / m) * grad(phase)

addpath('helpers/')

%% Parameters
m22      = 1;                             % (m/ 10^-22 eV)
Lbox     = 20;                            % kpc
N        = 128;  64;128;256;              % resolution
Tfinal   = 4;                             % kpc/(km/s) ~ 978 Myr
Nout     = 40;   400;                     % number of output
f15      = Inf;  0;1;2;4;Inf;             % (f/10^15 GeV)

output   = ['output/f' num2str(f15) 'L' num2str(Lbox) 'T' num2str(Tfinal) 'n' num2str(Nout) 'r' num2str(N) '/'];
makePlot = true; false;



%% Constants
hbar = 1.71818131e-87;       % hbar / (mass of sun * (km/s) * kpc)
m = m22 * 8.96215327e-89;    % 10^-22 eV / c^2 / mass of sun
G = 4.3022682e-6;            % G/((km/s)^2*kpc/mass of sun)
c = 299792.458;              % c / (km/s)
f = f15 * 8.05478166e-32;          % 10^15 GeV/((km/s)^2*mass of sun)

m_per_hbar = m/hbar;
psi4fac = (hbar*c/(f*m)^(2/3))^3/16;      % (hbar^3*c^3/(16*f^2*m^2));
psi6fac = hbar^2 * (hbar*c/(f*m))^4/288;  % hbar^6 * (c/(f*m))^4/288

%bfac1 = c*(m*f/hbar)^2/hbar;         % (m^2*c*f^2/hbar^3);
%bfac2 = hbar*c/2 * (hbar/(f*m))^2;   % (hbar^3*c/(2*f^2*m^2))
%bfac3 = hbar*sqrt(2*hbar*c)/(f*m);   % sqrt(2*hbar^3*c/(f^2*m^2));



%% setup


% spacing & grid
dx = Lbox / N;
xlin = ((0:N-1)' + 0.5) * dx;   % Note, x=0 and x=Lbox are the same point!

%assert(dx < 1); % Make sure we have better than 1 kpc resolution!

[x, y, z] = meshgrid(xlin, xlin, xlin);



% IC
rng(42);
Ncores = 2^randi([2,5]); % \in [4, 32]
rc = .2 + .8*rand(Ncores,1); % \in [.2,1]  % orig 8..50
ctr = Lbox * rand(Ncores,3);
rho = zeros(size(x));
for i = 1:Ncores
    for j = 1:N
        rho(:,:,j) = rho(:,:,j) + solitonProfile( sqrt( ...
            periodicDist(x(:,:,j),ctr(i,1),Lbox).^2 + ...
            periodicDist(y(:,:,j),ctr(i,2),Lbox).^2 + ...
            periodicDist(z(:,:,j),ctr(i,3),Lbox).^2 ) ...
            , rc(i), m22 );
    end
end


psi = sqrt(rho);
clear rho;


Msoliton =  2.2e8 ./ rc / m22^2;

Mtot = sum(abs(psi(:)).^2) * dx^3;       % total mass
rhobar = Mtot / Lbox^3;                  % average density

[max(Msoliton) Mtot]
%stop

cx = 0.5*Lbox;
cy = 0.5*Lbox;
cz = 0.5*Lbox;
%r = sqrt((x-cx).^2 + (y-cy).^2 + (z-cz).^2);
clear xlin;
clear x;
clear y;
clear z;



%%


% fourier space variables
fftw('planner','measure');
klin = (-N/2:N/2-1)' * (2*pi/Lbox);
[kx, ky, kz] = meshgrid(klin, klin, klin);
kSq = fftshift(kx.^2 + ky.^2 + kz.^2);
clear klin;
clear kx;
clear ky;
clear kz;


% initialize the potential
V = zeros(size(psi));

t = 0;
dt = 0;


% Load any previously saved data, if any
loadedPreviousSave = false;
for snapnum = Nout:-1:0
    filename = [output 'snap' sprintf('%.04d',snapnum) '.h5'];
    if exist(filename,'file')
        t = hdf5read(filename, '/time');
        psi = hdf5read(filename, '/psiRe') + 1.i * hdf5read(filename, '/psiIm');
        loadedPreviousSave = true;
        break
    end
end

% else, save first snapshot
if ~loadedPreviousSave
    if ~exist(output,'dir')
        mkdir(output);
    end
    snapnum = 0;
    filename = [output 'snap' sprintf('%.04d',snapnum) '.h5'];
    % if exist(filename,'file')
    %     delete(filename);
    % end
    hdf5write(filename, '/time', double(t))
    hdf5write(filename, '/m22', double(m22), 'WriteMode', 'append')
    hdf5write(filename, '/Lbox', double(Lbox), 'WriteMode', 'append')
    hdf5write(filename, '/psiRe', double(real(psi)), 'WriteMode', 'append')
    hdf5write(filename, '/psiIm', double(imag(psi)), 'WriteMode', 'append')
end

snapnum = snapnum + 1;
saveNextTurn = 0;

clear Ncores;
clear rc;
clear ctr;

if makePlot
    fh = figure(1);
    set(fh,'position',[0,0,400,400]);
end

%stop

%% simulation
tic;
while t < Tfinal
    
    % potential - (1/2) kick
    psi = exp(-1.i * dt/2 * m_per_hbar * V).*psi;
    
    % kinetic - drift
    psi = fftn(psi);
    psi = exp(dt * (1/m_per_hbar)/2 * -1.i * kSq) .*psi;
    psi = ifftn(psi);
    
    % potential - (1/2) kick
    V = 4*pi*G * (abs(psi).^2 - rhobar);
    V = -fftn(V);
    V = V ./ ( kSq  + (kSq==0) );
    V = ifftn(V);
    %add SI
    V = V - 2*psi4fac.*abs(psi).^2 + 3*psi6fac.*abs(psi).^4;
    
    psi = exp(-1.i * dt/2 * m_per_hbar * V).*psi;
    t = t + dt
    
    dt = min( m_per_hbar/6*dx^2, 1./(m_per_hbar*max(abs(V(:)))) );
    
    
    if t+dt >= snapnum * Tfinal/Nout
        dt = snapnum * Tfinal/Nout - t;
        assert(dt > 0);
        saveNextTurn = 1;
    end
    
    % save snapshot
    if saveNextTurn %t > snapnum * Tfinal / Nout
        filename = [output 'snap' sprintf('%.04d',snapnum) '.h5'];
        %         if exist(filename,'file')
        %             delete(filename);
        %         end
        hdf5write(filename, '/time', double(t))
        hdf5write(filename, '/m22', double(m22), 'WriteMode', 'append')
        hdf5write(filename, '/Lbox', double(Lbox), 'WriteMode', 'append')
        hdf5write(filename, '/psiRe', double(real(psi)), 'WriteMode', 'append')
        hdf5write(filename, '/psiIm', double(imag(psi)), 'WriteMode', 'append')
        
        snapnum = snapnum + 1;
        saveNextTurn = 0;
        
        %% plot
        if makePlot
            figure(fh);
            imagesc(log10(mean(abs(psi).^2,3)))
            caxis([3 9])
            axis off
            axis square
            colormap(inferno)
            pause(0.0001);
        end
    end
    
    
end
toc;

