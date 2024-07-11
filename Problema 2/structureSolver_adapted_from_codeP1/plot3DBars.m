function plot3DBars(x,Tn,scale,t,u,sig,units)

% Units
if nargin < 7
    units = 'Pa';
end

% DOFs per node
ni = 3;

% Precomputations
X = x(:,1);
Y = x(:,2);
Z = x(:,3);
smin = min(sig);
smax = max(sig);

% Determine how many timesteps
if exist('t','var') && size(sig,2)==length(t) && size(u,2)==length(t) && length(t)>1
    ntimes = length(t);
else
    ntimes = 1;
end

% Plot first (or unique) timestep
k = 1;
Ux = scale*u(1:ni:end,k);
Uy = scale*u(2:ni:end,k);
Uz = scale*u(3:ni:end,k);
% Open plot window
figure; box on; hold on; axis equal; view(45,30); axis tight;
if ntimes > 1 % Add slider and play buttons
    sl = uicontrol('Style','slider','Position',[25 22.5 335 15],...
                   'Min',1,'Max',length(t),'Value',1,'Callback',@plotSelectTime);
    pl = uicontrol('Style','pushbutton','Position',[370 20 50 20],...
                   'String','Play','Callback',@plotAllTimes);
end

% Plot undeformed structure
plot3(X(Tn'),Y(Tn'),Z(Tn'),'-', ...
    'Color',0.5*[1,1,1],'LineWidth',1);

% Plot deformed structure with colorbar for stresses 
p = patch(X(Tn')+Ux(Tn'),Y(Tn')+Uy(Tn'),Z(Tn')+Uz(Tn'),[sig(:,k)';sig(:,k)'], ...
    'EdgeColor','interp','LineWidth',2);

% Colorbar settings
clim(max(abs([smin,smax]))*[-1,1]);
n = 128; % Number of rows
c1 = 2/3; % Blue
c2 = 0; % Red
s = 0.85; % Saturation
c = hsv2rgb([c1*ones(1,n),c2*ones(1,n);1:-1/(n-1):0,1/n:1/n:1;s*ones(1,2*n)]');
colormap(c); 
cb = colorbar;

% Add labels
if ntimes > 1 || exist('t','var')
    title(sprintf('Time = %g sec | scale = %g (\\sigma_{min} = %.3g %s | \\sigma_{max} = %.3g %s)',t(k),scale,smin,units,smax,units)); 
else
    title(sprintf('scale = %g (\\sigma_{min} = %.3g %s | \\sigma_{max} = %.3g %s)',scale,smin,units,smax,units)); 
end
xlabel('x (m)'); 
ylabel('y (m)');
cb.Label.String = sprintf('Stress (%s)',units); 
set(gca,'color','none','xcolor','none','ycolor','none','zcolor','none');

% Plot all times
if ntimes > 1
    plotAllTimes
end

function plotSelectTime(src,~)
    if ntimes > 1
        k_ = round(src.Value);
    end
    Ux = scale*u(1:ni:end,k_);
    Uy = scale*u(2:ni:end,k_);
    Uz = scale*u(3:ni:end,k_);
    p.XData = X(Tn')+Ux(Tn');
    p.YData = Y(Tn')+Uy(Tn');
    p.ZData = Z(Tn')+Uz(Tn');
    p.CData = [sig(:,k_)';sig(:,k_)'];
    title(sprintf('Time = %g sec | scale = %g',t(k_),scale)); 
    drawnow;
end

function plotAllTimes(src,~)
    for k_ = 1:ntimes
        sl.Value = k_;
        plotSelectTime(sl);
    end
end

end