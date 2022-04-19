function [c_hat] = chan_estimator(Y,pilot,L)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    %%% Directly compute channel estimate using 11.54 from Dr. Farhang's book %%%
    %Generate S matrix
%     pilot_L=conv(expander([pilot; pilot],L),pT);
%     pilot_L=pilot_L((length(pT)-1)/2+1:(length(pT)-1)/2+N*L);
%     pilot_L2=conv(beacon(L*N:2*L*N-1),pR);
%     pilot_L2=conv(beacon,pR);
% %     cutoff=(length(pR)-1)/2;
%     pilot_L2=pilot_L2(L*N:2*L*N-1);
    
%     LPF=r_cos_p(10*L,L,alpha_srrc);
%     pilot_L3=conv(beacon,LPF);
%         cutoff=(length(LPF)-1)/2;
%     pilot_L3=pilot_L3(cutoff+1:end-cutoff);
%     pilot_L3=pilot_L3(1:L*N);
%     s=circshift(pilot_L,1);
%     S=zeros(N*L,length(s));
%     for i=1:N*L
%         S(i,:)=s;
%         s=circshift(s,1);
%     end

    N=length(pilot);
    s=circshift(pilot,1);
    S=zeros(N,N);
    for i=1:N
        S(i,:)=s;
        s=circshift(s,1);
    end
    
%     s_L=circshift(pilot_L3,1);
%     S_L=zeros(L*N,L*N);
%     for i=1:L*N
%         S_L(i,:)=s_L;
%         s_L=circshift(s_L,1);
%     end
    %average signal over several periods, slightly more accurate
    %Compute channel estimate C = (S*S^H)^-1 * (S*conj(Y))
    %     y_avg=reshape(yA(2*N*L:length(yA))-2*N,N,length(yA(2*N:length(yA)))/N);
    
    %4 decimations experiment
    for i=1:L
%         S_t=S_L(i:L:end,i:L:end);
        y=Y(i:L:end);
        y_avg=reshape(y(2*N:22*N-1),N,20);
        y_avg=mean(y_avg,2);
%         plot(abs(y_avg));hold on;
%         c_hat_dec(i,:)=inv(S_t*S_t')*(S_t.'*conj(y_avg));
        cA_t(i,:)=inv(S*S')*(S.'*conj(y_avg));
        
%         y=yB_r(i:L:end);
%         y_avg=mean(reshape(y(2*N:22*N-1),N,20),2);
%         cB_t(i,:)=inv(S*S')*(S.'*conj(y_avg));
%         
%         y=yE_r(i:L:end);
%         y_avg=mean(reshape(y(2*N:22*N-1),N,20),2);
%         cE_t(i,:)=inv(S*S')*(S.'*conj(y_avg));
%         ce_hat_dec(i,:)=inv(S_t*S_t')*(S_t.'*conj(y_avg));
    end
%     hold off;
    c_hat=cA_t(:)';
%     cB_hat=cB_t(:)';
%     cE_hat=cE_t(:)';
   
    
%     cA_hat_es=c_hat_dec(:)';
%     cE_hat_es=ce_hat_dec(:)';
    
    
    %Long experiment
    
    
%     y_avg_long=reshape(yA_r(2*L*N:22*L*N-1),L*N,20);
%     y_avg_long=mean(y_avg_long,2);
%     c_hat_L=inv(S_L*S_L')*(S_L.'*conj(y_avg_long));
%     
%     %1 decimation experiment
%     yA_avg=reshape(yA_dec(2*N:22*N-1),N,20);
%     yA_avg=mean(yA_avg,2);
%     cA_hat=inv(S*S')*(S.'*conj(yA_avg)); %Using equation (11.54) from software radio textbook
    
%     figure(1),plot(abs(cA_hat)),title('4 decimations, original pilot'),shg;
%     figure(2),plot(abs(c_hat_L)),title('Long pilot'),shg;
%     figure(3),plot(abs(cA_hat_dec)),title('4 decimations on pilot'),shg;
%     
    
    
%     
%     yB_avg=reshape(yB(2*N:22*N-1),N,20);
%     yB_avg=mean(yB_avg,2);
%     cB_hat=inv(S*S')*(S.'*conj(yB_avg));
%     
%     yE_avg=reshape(yE(2*N:22*N-1),N,20);
%     yE_avg=mean(yE_avg,2);
%     cE_hat=inv(S*S')*(S.'*conj(yE_avg));
%     
    
%         figure(1)
%         spec_analysis(cA_hat,fs)
%         title('A')
    %     figure(2)
    %     spec_analysis(cB_hat,fs)
    %     title('B')
    %     figure(3)
    %     spec_analysis(cE_hat,fs)
    %     title('E')

end

