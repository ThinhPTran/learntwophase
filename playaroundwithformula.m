clear all 
close all
clc

N=10000;

Vsg = 1000*rand(N,1);
Vsl = 10*ones(N,1); 

for i = 1:N
  Hg(i) = 0.5249*Vsg(i)^(0.2089+0.0844*Vsl(i))*exp((-0.0175 - 0.0103*Vsl(i))*(log(Vsg(i)))*(log(Vsg(i)))-0.179*Vsl(i));
end

max(Hg)
min(Hg)
hist(Hg)






  





