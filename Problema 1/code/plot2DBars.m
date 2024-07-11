function plot2DBars(data,x,Tn,u,sig,scale,units)

% Precomputations
X = x(:,1);
Y = x(:,2);
Ux = scale*u(1:data.ni:end,1);
Uy = scale*u(2:data.ni:end,1);
smin = min(sig);
smax = max(sig);

% Open plot window
figure; box on; hold on; axis equal;
% Plot undeformed structure
plot(X(Tn'),Y(Tn'),'-', ...
    'Color',0.85*[1,1,1],'LineWidth',1);
% Plot deformed structure with colorbar for stresses 
patch(X(Tn')+Ux(Tn'),Y(Tn')+Uy(Tn'),[sig';sig'], ...
    'EdgeColor','interp','LineWidth',2);
% Colorbar settings
clims = get(gca,'clim');
clim(max(abs(clims))*[-1,1]);
n = 128; % Number of rows
c1 = 2/3; % Blue
c2 = 0; % Red
s = 0.85; % Saturation
c = hsv2rgb([c1*ones(1,n),c2*ones(1,n);1:-1/(n-1):0,1/n:1/n:1;s*ones(1,2*n)]');
colormap(c); 
cb = colorbar;

% Add labels
title(sprintf('scale = %g (\\sigma_{min} = %.3g %s | \\sigma_{max} = %.3g %s)',scale,smin,units,smax,units)); 
xlabel('x (m)'); 
ylabel('y (m)');
cb.Label.String = sprintf('Stress (%s)',units); 

end