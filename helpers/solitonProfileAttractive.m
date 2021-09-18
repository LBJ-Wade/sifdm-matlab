function rhoc = solitonProfileAttractive( r, rc, m22, f15 )
%SOLITONPROFILE returns soliton density at radius r
%   units: kpc, Msun, km/s
%   input: sample radius r (kpc), core radius rc (kpc), m22
%   total mass is ~ 2.20e8 /rc / m22^2


hbar = 1.71818131e-87;       % hbar / (mass of sun * (km/s) * kpc)
m = m22 * 8.96215327e-89;    % 10^-22 eV / c^2 / mass of sun
G = 4.3022682e-6;            % G/((km/s)^2*kpc/mass of sun)
c = 299792.458;              % c / (km/s)
f = f15 * 8.05478166e-32;          % 10^15 GeV/((km/s)^2*mass of sun)


beta0 = hbar*c^5/(32*pi*G*f^2) * 6.9e-9; % lambda^2 \propto rho0^(1/4) \propto (rc/kpc)^-2
beta = beta0 * rc.^-2;



rho0 = 1.9e7 * m22^-2 * rc^-4;    % Msun/kpc^3
%rhoc = rho0 ./ (1 + 0.091 * (r/rc).^2).^8;
i1 = tanh(beta/5);
i2 = tanh(beta);
i3 = tanh(sqrt(beta))^2;
rhoc = rho0 * (1 + (1+2.6*i1)*0.091 * (r*(sqrt(1+beta))/rc).^(2-i2/5) ).^(-8+22/5*i3);

end

