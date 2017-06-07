clear all 
close all
clc


% A well with 3000 m long
% A segement is 1 m long
% dt = 10s
% Pr = 2800 psi = 1.9305*10^7 Pascal (N/m^2)
% Assume vertical well
% At the beginning, hydrostatic pressure = Pr
Lwell = 3000 % 3000 m
N = 3001 % 3001 points. 1m per segment
dt = 10 % s 
dL = 1
verd = (0:N-1)';
Pr = 2800*6894.75729 % Change to Pascal
Pwf = Pr % Balanced, no flow

% We assume that liquid is incompressable so rhoL doesn't change
% Gas is compressable so rhoG can change
rhoL = 1050 % kg/m^3 rho water
rhoG = zeros(N,1); % At the beginning there is no gas in the tube.

g = 9.8 % m/s^2 
Pth = 101325 % Atmospher pressure in Pascal
disp('current level of water in tubing');
tllev = 3000 - (Pwf-Pth)/(g*rhoL)

% At the beginning there is only liquid in the tube
HG = zeros(N,1); 
HL = zeros(N,1); 

for iter = floor(tllev)+2:N
    HL(iter) = 1.0;
end    

% At the beginning there is no flow in the tube
Vsl = zeros(N,1); 
Vsg = zeros(N,1); 



% Hydrostatic pressure
P=zeros(N,1);
for iter = 1:floor(tllev)+1
    P(iter)=Pth;
end

for iter = floor(tllev)+2:N
    P(iter) = P(iter-1) + rhoL*g;
end

figure(1)
plot(verd,P/6894.75729)
legend('Hydrostatic pressure in Psi');

figure(2)
plot(verd,HL);
legend('Liquid Holdup');

figure(3);
plot(verd,HG);
legend('Void fraction'); 












   

    











