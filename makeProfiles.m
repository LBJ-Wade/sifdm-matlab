close all
clear all
clc
% Philip Mocz (2020), Princeton University
% Make Radial Profiles

% Internal units:
% [L] = kpc
% [M] = Msun
% [E] = Msun (km/s)^2


%stop

%% simulation ID
m22      = 1;                            % (m/ 10^-22 eV)
Lbox     = 20;                           % kpc
Ns        = [256 400];                           % resolution
Tfinal   = 4;                            % kpc/(km/s) ~ 978 Myr
Nout     = 400;                            % number of output

f15s = [1 1.2 2 Inf];  % (f/10^15 GeV)


cc = 1;
for N = Ns
    ff = 1;
    for f15  = f15s
        
        
        
        
        snapdir   = [ 'output/f' num2str(f15) 'L' num2str(Lbox) 'T' num2str(Tfinal) 'n' num2str(Nout) 'r' num2str(N) '/'];
        
        clim = [5 9];
        
        
        addpath('helpers/')
        
        % constants
        hbar = 1.71818131e-87;
        G = 4.3022682e-6;
        
        
        %% Plot Radial Profiles
        snapnum = 400;
        %fh = figure;
        %set(fh,'position',[0 0 600 600],'PaperPosition',[0 0 6 6]);
        
        
        
        [ t, m22, Lbox, N, psi ] = readsnap( snapdir, snapnum );
        m = m22 * 8.96215327e-89;
        
        
        
        
        %%
        colormap(inferno);
        
        figure;
        imagesc([0 Lbox],[0 Lbox],log10(mean(abs(psi).^2,3)))
        caxis(clim)
        axis off
        axis square
        set(gca,'ydir','normal')
        
        set(gca, 'unit', 'normalize')
        set(gca, 'position', [0 0 1 1]);
        
        
        %% Grid
        
        dx = Lbox / N;
        xlin = ((1:N)-0.5)*dx;
        [xx,yy, zz] = meshgrid(xlin,xlin,xlin);
        
        %% Get radial profile
        
        maxrho = max(abs(psi(:)).^2);
        
        idx = find(abs(psi).^2 == maxrho);
        [ic,jc,kc] = ind2sub([N N N], idx);
        [abs(psi(ic,jc,kc)).^2 abs(psi(idx)).^2 maxrho]
        
        %%
        xc = xx(idx);
        yc = yy(idx);
        zc = zz(idx);
        
        % use sub-pixel centering
        norm = 0;
        xc = 0;
        yc = 0;
        zc = 0;
        for ii = -1:1
            for jj = -1:1
                for kk = -1:1
                    norm = norm + abs(psi(ic+ii,jc+jj,kc+kk)).^2;
                    xc = xc + abs(psi(ic+ii,jc+jj,kc+kk)).^2     * xx(ic+ii,jc+jj,kc+kk);
                    yc = yc + abs(psi(ic+ii,jc+jj,kc+kk)).^2     * yy(ic+ii,jc+jj,kc+kk);
                    zc = zc + abs(psi(ic+ii,jc+jj,kc+kk)).^2     * zz(ic+ii,jc+jj,kc+kk);
                end
            end
        end
        xc = xc / norm;
        yc = yc / norm;
        zc = zc / norm;

        ddx = xx - xc;
        ddy = yy - yc;
        ddz = zz - zc;
        
        ddx(ddx > Lbox/2) = ddx(ddx > Lbox/2) - Lbox;
        ddx(ddx < -Lbox/2) = ddx(ddx < -Lbox/2) + Lbox;
        ddy(ddy > Lbox/2) = ddy(ddy > Lbox/2) - Lbox;
        ddy(ddy < -Lbox/2) = ddy(ddy < -Lbox/2) + Lbox;
        ddz(ddz > Lbox/2) = ddz(ddz > Lbox/2) - Lbox;
        ddz(ddz < -Lbox/2) = ddz(ddz < -Lbox/2) + Lbox;
        
        rr = sqrt(ddx.^2 + ddy.^2 + ddz.^2);
        
        %%
        H = myBinMean1D(rr(:), abs(psi(:)).^2, xlin);
        
        
        r_profile{cc,ff} = xlin;
        rho_profile{cc,ff} = H;
        
        %%
        figure(101);
        loglog(xlin,H)
        hold on
        
        ff = ff + 1;
        
    end
    cc = cc + 1;
end

close all

%% Make a plot

figure(10);

%
m22 = 1;

f15 = 1.2;Inf;1.2;Inf;1.2;Inf;2;
rc = 0.2;

for rc = 0.136;[0.181 0.14];[.1 .15 .2];0.14;0.1.^(0.1:0.1:2); % [0.1 0.2 0.4 0.8]
    
    rrr = dx/16:dx/16:Lbox;
    
    rhoc = solitonProfileAttractive( rrr, rc, m22, f15 );
    loglog(rrr,rhoc,'color',[1 1 1]*rc*3,'linewidth',6)
    hold on
    
end

f15 = 2;
rc = 0.165;
    rhoc = solitonProfileAttractive( rrr, rc, m22, f15 );
    loglog(rrr,rhoc,'color',[1 1 1]*rc*3,'linewidth',6)


    f15 = Inf;
rc = 0.18;
    rhoc = solitonProfileAttractive( rrr, rc, m22, f15 );
    loglog(rrr,rhoc,'color',[1 1 1]*rc*3,'linewidth',6)
    
    
loglog(rrr,8e8*exp(-rrr/.2))
loglog(rrr,5e8*exp(-rrr/.3))
%loglog(rrr,1.8e10*exp(-rrr.^2/.05))

loglog(rrr, 3e7*rrr.^-2)

cc = 1;
for N = Ns
    ff = 1;
    for f15  = f15s
        
        loglog(r_profile{cc,ff},rho_profile{cc,ff},'linewidth',N/200,'color',[1-(ff-1)/3 0 (ff-1)/3])
        hold on
        
        ff = ff + 1;
    end
    cc = cc +1;
end

hold off
axis([2e-2 1e1 1e5 5e12])




