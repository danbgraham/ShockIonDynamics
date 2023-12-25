Vsw = -320e3; 
n0 = 4*1e6; 
n1 = 13*1e6; 
B0 = 13.5*1e-9; 
B1 = 45.5*1e-9; 
nsw = n1-n0;

P0 = 0.0425e-9;
P1 = 0.0575e-9;

l = 10e3; 

xvec = [-10e2:1:5e2]*1e3;
[xx,yy] = meshgrid(xvec,xvec);
Byxx = -B0*tanh(xx/l)+B1;

By = -B0*tanh(xvec/l)+B1;
ni = -n0*tanh(xvec/l)+n1;

Units = irf_units;
mu0 = Units.mu0;
e = Units.e;
mi = Units.mp;

Jz = -B0*sech(xvec/l).^2/(mu0*l);
Vx = Vsw*nsw./ni;

Emultiplier = 1;

Ex = -Jz.*By./(e*ni) + P0*sech(xvec/l).^2./(e*ni*l);
Ez = -Vx.*By;

dx = 1e3;
phi = cumsum(Ex)*dx;

opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
%%
if 1
dv = 10e3;
vxvec = -1300e3:dv:700e3;
vzvec = -1000e3:dv:1000e3;
xpos = (-2000:10:200)*1e3;

plotpos1 = -200e3;
plotpos2 = 30e3;
[~,idx1] = min(abs(xpos-plotpos1));
[~,idx2] = min(abs(xpos-plotpos2));

[VX,VZ] = meshgrid(vxvec, vzvec);
fxz = zeros(size(VX));
fxvxvz = repmat(fxz,1,1,length(xpos));
fxvxvz = squeeze(permute(fxvxvz,[3 1 2]));

fxvxvztemp = squeeze(zeros(size(fxvxvz(1,:,:))));

Ti = 12;
vi = sqrt(2*e*Ti/mi/4);
nsw = 9e6;
vsw = -320e3;
spos = 200e3;
y0 = 0;
z0 = 0;
vy0 = 0;

%%

for kk = 1:length(xpos)

x0 = xpos(kk);  
tic
for ii = 1:length(vxvec)
	parfor jj = 1:length(vzvec)
    [t1b,d1b] = ode45(@temp,[20 0],[x0; vxvec(ii); y0; vy0; z0; vzvec(jj)],opts);
    if d1b(end,1) > spos
      [a,~] = find(d1b(:,1) > spos);
      idx = a(1)-1;
      fxvxvztemp(jj,ii) = nsw/(pi*vi^2)*exp((-(d1b(idx,2) - vsw)^2-d1b(idx,6)^2)/vi^2);
    else
      fxvxvztemp(jj,ii) = 0;
    end
  end

fxvxvz(kk,:,:) = fxvxvztemp;
end
toc

kk
end

fdistalpha = struct('xpos',xpos,'vxvec',vxvec,'vzvec',vzvec,'fxvxvz',fxvxvz,'Emultiplier',Emultiplier,'dv',dv);
save('fdistalpha.mat','fdistalpha')

end

%%
function dydt = temp(t,q)

l = 10e3;
vsw = -320e3;

n0 = 4*1e6; 
n1 = 13*1e6; 
B0 = 13.5*1e-9; 
B1 = 45.5*1e-9; 

P0 = 0.0425e-9;
P1 = 0.0575e-9;

mu0 = 1.2566e-06;
eomi = 9.5788e+07/2;
e = 1.6022e-19;

By = -B0*tanh(q(1)/l)+B1;
ni = -n0*tanh(q(1)/l)+n1;
Jz = -B0*sech(q(1)/l).^2/(mu0*l);
Ex = -Jz.*By./(e*ni) + P0*sech(q(1)/l).^2./(e*ni*l);
Vx = vsw*(n1-n0)./ni;
Ey = 0;
Ez = -Vx*By;

dydt = [q(2) ; eomi*Ex - eomi*q(6)*By; ...
        q(4) ; eomi*Ey; ...
        q(6) ; eomi*Ez + eomi*q(2)*By];
end