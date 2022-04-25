function [chan] = create_rayleigh(M_paths,delay,Tb)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    cA_gains=zeros(1,M_paths);
%     cE_gains=zeros(1,M);
    
%     A=exp(-[0:M_paths-1]*Tb/2/delay);
    A=zeros(1,M_paths);
    A(1)=1-exp(-Tb/delay);
    A=A(1)*exp(-[0:M_paths-1]*Tb/delay);
    
    
%     A1=exp(-[0:M-1]*Ts/delay);
    
    %exp(t/tau)

    cA_gains=sqrt(A/2).*(randn(1,M_paths)+1i*randn(1,M_paths));
%     cE_gains=sqrt(1/2)*(randn(1,M)+1i*randn(1,M));
       


    
%     chan=cA_gains.*A;
    chan=cA_gains;
    
    
    
    %pdf x(n)/sigma2 * exp(- x(n)^2/(2*sigma2)
%     cB=cA;
%     cE=cE_gains.*A;
    
%     t_A=abs(cA).^2; %Power delay profile? eq (4.15)
%     t_E=abs(cE).^2;
end

