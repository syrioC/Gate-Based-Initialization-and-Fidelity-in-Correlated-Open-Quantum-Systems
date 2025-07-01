function [T,bcf]=BCFkT(tmax,nmax,omegac,s,kT);
%generates bcf vs time, in time stepts tmax/N and range [0,tmax];
%N=2^15;
n0=32;
n0=256;
N=nmax;
% %increases Nyquist frequency
tmax=tmax*n0;N1=N*n0;%leeways
dt=tmax/N1;
T=(0:N1)'*dt;
omegamax=pi/dt;  %Nyquist frequency
domega=omegamax/N1;
omega=(1:N1)'*domega;

% Nw=N1/16;
% windw=ones(N1+1,1);
% windw(N1+1:-1:N1-Nw+2)=(1:Nw)/Nw;
omeg=omega/omegac;
tic
J=2*pi*omegac*omeg.^s.*exp(-omeg); % power law SD
J=[0;J];
omega=[0;omega];
toc
if kT == 0
    J(N1+2:2*N1)=0;
else
    Eom=exp(-omega(2:N1+1)/kT);
    J(2:N1+1)=J(2:N1+1)./(1-Eom);    
    J(N1+2:2*N1)=J(N1:-1:2).*Eom(N1:-1:2);
    J(1)=4*(pi*kT*2^(-s)/s)*(domega/omegac)^(s-1);%divided by 2 because of fft algorithm. 
end
clear Eom
toc
bcf=fft(J)/tmax;
toc
%bcf=bcf(1:numel(T));    %removes negative times
T=T(1:N+1);
bcf=bcf(1:N+1);
toc
end


