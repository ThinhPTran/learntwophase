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
Nstep = 1000; 
N = 3001 % 3001 points. 1m per segment
dt = 10 % s 
dL = 1
verd = (0:N-1)';
Pr = 2800*6894.75729 % Change to Pascal
Pwf = Pr*ones(1,Nstep); % Balanced, no flow

% We assume that liquid is incompressable so rhoL doesn't change
% Gas is compressable so rhoG can change
rhoL = 1050; % kg/m^3 rho water
rhoG = zeros(N,Nstep); % At the beginning there is no gas in the tube.

g = 9.8 % m/s^2 
Pth = 101325 % Atmospher pressure in Pascal
disp('current level of water in tubing');
tllev = 3000 - (Pwf(1)-Pth)/(g*rhoL(1))

% At the beginning there is only liquid in the tube
HG = zeros(N,1); 
HL = zeros(N,1); 

for iter = floor(tllev)+2:N
    HL(iter) = 1.0;
end

for iter = 1:N
   HG(iter) = 1.0 - HL(iter);  
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
    P(iter) = P(iter-1) + rhoL(1)*g;
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

rhoLHL=zeros(N,Nstep); 
rhoGHG=zeros(N,Nstep); 
rhoLVsl=zeros(N,Nstep); 
rhoGVsg=zeros(N,Nstep); 
rhoLVslrhoGVsg=zeros(N,Nstep); 

% Initialize for all these information 
for iter = 1:N
    rhoLHL(iter,1) = rhoL*HL(iter,1);
    rhoGHG(iter,1) = rhoG(iter,1)*HG(iter,1);
    rhoLVsl(iter,1) = rhoL*Vsl(iter,1); 
    rhoGVsg(iter,1) = rhoG(iter,1)*Vsg(iter,1); 
    rhoLVslrhoGVsg(iter,1) = rhoL*Vsl(iter,1) + rhoG(iter,1)*Vsg(iter,1); 
end


fldirect=1;


for iter = 2:Nstep
    
    if (Pwf<=Pr)
        fldirect=1;
    else
        fldirect=-1; 
    end
    
    % Solve Mass conservation equations
    rhoGHG(N,iter) = rhoGHG(N,iter-1);
    rhoLHL(N,iter) = rhoLHL(N,iter-1);
    HL(N,iter) = rhoLHL(N,iter)/rhoL;
    HG(N,iter) = 1.0 - HL(N,iter); 
    if (HG(N,iter)==0)
      rhoG(N,iter) = 0;
    else 
      rhoG(N,iter) = rhoGHG(N,iter)/HG(N,iter);
    end
    
    for jter = N-1:-1:1
        rhoGHG(jter,iter) = rhoGHG(jter,iter-1) - (rhoGVsg(jter,iter-1) - rhoGVsg(jter+1,iter-1))*dt/dL; 
        rhoLHL(jter,iter) = rhoLHL(jter,iter-1) - (rhoLVsl(jter,iter-1) - rhoLVsl(jter+1,iter-1))*fldirect*dt/dL;
        HL(jter,iter) = rhoLHL(jter,iter)/rhoL;
        HG(jter,iter) = 1.0 - HL(jter,iter);
        if (HG(jter,iter)==0)
           rhoG(jter,iter)=0;
        else
           rhoG(jter,iter)=rhoGHG(jter,iter)/HG(jter,iter);  
        end
    end
    
    % Solve Momentum equation
    rhoLVslrhoGVsg(N,iter) = rhoLVslrhoGVsg(N,iter-1);
    
    for jter = N-1:-1:1
       if (HG(jter,iter-1)==0) 
         rhoLVslrhoGVsg(jter,iter) = rhoLVslrhoGVsg(jter,iter-1) + (rhoL*Vsl(jter,iter-1)*Vsl(jter,iter-1)-rhoL*Vsl(jter+1,iter-1)*Vsl(jter+1,iter-1))*fldirect*dt/dL;
       elseif (HL(jter,iter-1)==0)   
         
       else 
           
       end
    end
    
    
end






















   

    











