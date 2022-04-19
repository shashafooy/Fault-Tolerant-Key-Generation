

clear all

% cA=[1 zeros(1,91) 0.3];
% cB=[1 zeros(1,91) 0.3]+randn(1,length(cA))*sqrt(0.00001);
% cE=[1 zeros(1,71) 0.4];
%%%% Channel Creation %%%%
L=4;
M=16;
fc=100000;
fs=130e6;
Ts=1/fs;
Tb=Ts*L;
alpha_srrc=0.5;
sigmanuc=0.05;
epsilon=1e-6;
iterations=10000;
iteratoins=500;

Rho_R=zeros(iterations,1);
Rho_AB=zeros(iterations/2,1);
Rho_AB_SPC=zeros(iterations/2,1);
Rho_AE=zeros(iterations/2,1);
Rho_AE_SPC=zeros(iterations/2,1);
SPC_counter=1;

% SNRo_B=zeros(iterations,1);
B_errors=zeros(iterations,1);
E_errors=zeros(iterations,1);
SNRo_B=zeros(iterations,1);


for iter=1:iterations
    pT=sr_cos_p(10*L,L,alpha_srrc);
%     pR=pT;
    pR=conj(pT(end:-1:1)); %more accurate in cases where Pr is slightly different
    
    %ZC sequence
    N=64;
    L=4;
    pilot=CycPilot(N-1);
    s=pilot;
    ss=[s;s;s;s;s];
    beacon=[ss;ss;ss;ss;ss];
    
    
    % cA=expander(cA,L);
%     beacon=conv(expander(beacon,L),pT);
    beacon=expander(beacon,L); %will get LP filtered when convolved with channel
    p=conv(pT,pR);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Channel creation %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Es=((beacon)'*beacon)/Ls;
    delay=50e-9;
    
    cA=create_rayleigh(M,delay,Tb);
    cB=cA;
    cE=create_rayleigh(M,delay,Tb);
    
    Rho_R(iter)=abs(cA*cE'/(norm(cA,2)*norm(cE,2)))^2;
%     continue
      
    A=exp(-[0:M-1]*Tb/delay/2);

    
    Rho_R(iter)=abs(cA*cE'/(norm(cA,2)*norm(cE,2)))^2;
%     Rho_R_gains(iter)=abs(t_A*t_E'/(norm(t_A,2)*norm(t_E,2)))^2;
    
%     continue
% cA=1;
% cE=1;  
    %%% Baseband
    cA_bb=conv(cA,pT); %cA_bb=cA_bb(1:L:end);
    cB_bb=cA_bb;
%     cE=cE.*exp(-1i*2*pi*[0:length(cE)-1]*Ts*fc);
    cE_bb=conv(cE,pT);% cE_bb=cE_bb(1:L:end);
    
%     pR=sqrt(L/2)*pR(1:L/2:end);

    %%% Used to compare to the estimate
    cA_true=conv(cA_bb,pR);
    cA_true=[cA_true;zeros(256-96,1)];
    cE_true=conv(cE_bb,pR);
    cE_true=[cE_true;zeros(256-96,1)];
    
    %%%%%%%%temporary%%%%%%%
%     Rho_R(iter)=abs(cA*cE'/(norm(cA,2)*norm(cE,2)));
%     corr=abs(corrcoef([cA_gains.' cE_gains.']));
%     Rho_R(iter)=corr(1,2);
%     continue
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Transmit ZC over Rayleigh channel %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %sigma = power/snr
    yA=conv(beacon,cA_bb);
    vA=sigmanuc*(randn(2*length(yA),1)+1i*randn(2*length(yA),1))/sqrt(2); %noise
    vA=conv(pT,vA);
    yA=yA+vA(1:2:2*length(yA));
    
    yB=conv(beacon,cB_bb);
    vB=sigmanuc*(randn(2*length(yB),1)+1i*randn(2*length(yB),1))/sqrt(2); %noise
    vB=conv(pT,vB);
    yB=yB+vB(1:2:2*length(yB));
    
    yE=conv(beacon,cE_bb);
    vE=sigmanuc*(randn(2*length(yE),1)+1i*randn(2*length(yE),1))/sqrt(2); %noise
    vE=conv(pT,vE);
    yE=yE+vE(1:2:2*length(yE));

    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% Channel Estimate %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    %Decimate Y by L
    yA_r=conv(yA,pR);
    yB_r=conv(yB,pR);
    yE_r=conv(yE,pR);
    yA_dec=decimator(conv(yA,pR),L);
    
    
    cA_hat=chan_estimator(yA_r,pilot,L);
    cB_hat=chan_estimator(yB_r,pilot,L);
    cE_hat=chan_estimator(yE_r,pilot,L);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Timing Alignment %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    L2=10;
    rc=r_cos_p(10*L2,L2,alpha_srrc);
    
%     cA_L=conv(expander(cA_hat,L),p);
%     cutoff=(length(p)-1)/2;
    cA_L=conv(expander(cA_hat,L2),rc);
    cutoff=(length(rc)-1)/2;
    
    
    cA_L=cA_L(cutoff+1:end-cutoff);
    Nc=length(cA_L);
%     cB_L=conv(expander(cB_hat,L),p);
    cB_L=conv(expander(cB_hat,L2),rc);
    cB_L=cB_L(cutoff+1:end-cutoff);
%     cB_L=conv(expander(c_hat_L,L),p);
    
%     cE_L=conv(expander(cE_hat,L),p);
    cE_L=conv(expander(cE_hat,L2),rc);
    cE_L=cE_L(cutoff+1:end-cutoff);
    
    
    %%% Find Path Candidates %%%
    
%     p_L=conv(expander(p,L),p); %Expand srrc by L
    p_L=conv(expander(p,L2),rc); %Expand srrc by L
    p_L=p_L(cutoff+1:end-cutoff);
%     p_L=p;
%     Nz=(Nc-length(srrc)/2)/2;
    z=zeros(1,Nc-length(p_L)+1);
    z(ceil(length(z)/2))=1;
    % long_p=conv(expander(p,L2),p);
    long_p=conv(z,p_L);
%     long_p=conv(z,srrc);
    long_p=circshift(long_p,-Nc/2+2);
%     long_p=[p_L; zeros(length(cA_L)-length(p_L),1)]';
    
    alpha=zeros(1,M);
    k=zeros(1,M);
    
    
    [~,idx_inc]=max(p_L);

    
    [alpha_A,k_A]=path_candidates(cA_L,M,long_p);
    [alpha_B,k_B]=path_candidates(cB_L,M,long_p);
    [alpha_E,k_E]=path_candidates(cE_L,M,long_p);
    
    %%%shift channels to be centered
    [~,idx]=max(abs(alpha_A));
%     idx=1;
    cA_L=circshift(cA_L,length(cA_L)/2-k_A(idx));
    k_A=mod(k_A-k_A(idx),Nc);
    
    [~,idx]=max(abs(alpha_B));
    cB_L=circshift(cB_L,length(cB_L)/2-k_B(idx));
    k_B=mod(k_B-k_B(idx),Nc);
    
    [~,idx]=max(abs(alpha_E));
%     idx=1;
    cE_L=circshift(cE_L,-k_E(idx)+1); %just shift eve to n=0
    k_E=mod(k_E-k_E(idx),Nc);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Alice and Bob alignment  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    kA_bar=round(sum(k_A.*abs(alpha_A).^2)/sum(k_A));
    kB_bar=round(sum(k_B.*abs(alpha_B).^2)/sum(k_B));
    
    k_D=Nc/2-kA_bar; %Sent to Bob
    
    k_ref=kB_bar+k_D;
    oldmin=1e6;
    ref_p=circshift(long_p,k_ref-1);
    scaling=norm(alpha_B,2);
    for i=1:M
        newmin=norm(abs(alpha_B(i)/scaling)*circshift(long_p,k_B(i)+Nc/2-1) -ref_p ,2)^2;
        if(newmin<oldmin)
            oldmin=newmin;
            alpha_B_sp=alpha_B(i);
            k_B_sp=k_B(i);
        end
    end
    cB_L=circshift(cB_L,k_B_sp); %align bob to alice
    
    %%% align Alice and Bob to n=0
    cA_L=circshift(cA_L,-Nc/2+1);
    cB_L=circshift(cB_L,-Nc/2+1);
    
    %%% Decimate
%     cutoff = (length(p)-1)/(2*L);
%     cA_hat_dec=decimator(cA_L,L);
%     cB_hat_dec=decimator(cB_L,L);
%     cE_hat_dec=decimator(cE_L,L);

    cA_hat_dec=decimator(cA_L,L2);
    cB_hat_dec=decimator(cB_L,L2);
    cE_hat_dec=decimator(cE_L,L2);

    %0 phase at n=0
    cA_hat_dec=cA_hat_dec*exp(-1i*angle(cA_hat_dec(1)));
    cB_hat_dec=cB_hat_dec*exp(-1i*angle(cB_hat_dec(1)));
    cE_hat_dec=cE_hat_dec*exp(-1i*angle(cE_hat_dec(1)));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Strongest Path Cancelation %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(mod(iter,2))
        Ns=length(cA_hat_dec);
        Nz=(Ns-length(p)/2)/2;
        z=zeros(Ns-length(p)+1,1);
        z(ceil(length(z)/2))=1;
        % long_p=conv(expander(p,L2),p);
        long_p_SPC=conv(z,p);
        long_p_SPC=circshift(long_p_SPC,-floor(Ns/2)+1).';
%         long_p_SPC=[p;zeros(length(cA_hat_dec)-length(p),1)]';
        
        
        
        cA_bar=cA_hat_dec-abs(alpha_A(1))*long_p_SPC;
        [~,idx]=min(k_B); %shift of 0 may not be at idx==1
        cB_bar=cB_hat_dec-abs(alpha_B(idx))*long_p_SPC;
        cE_bar=cE_hat_dec-abs(alpha_E(1))*long_p_SPC;
    else
        cA_bar=cA_hat_dec;
        cB_bar=cB_hat_dec;
        cE_bar=cE_hat_dec;
    end
    %%%%%%%%%%%%%%%%%%%%%%
    %%% Key Generation %%%
    %%%%%%%%%%%%%%%%%%%%%%

    CA=fftshift(fft(cA_bar,length(cA_bar)));
    CB=fftshift(fft(cB_bar,length(cB_bar)));
    CE=fftshift(fft(cE_bar,length(cE_bar)));
 
    passband=[length(CA)/2-N/2+1:length(CA)/2+N/2];
    
    gamma_A=(CA(passband)/norm(CA(passband),2)).';
    gamma_B=(CB(passband)/norm(CB(passband),2)).';
    gamma_E=(CE(passband)/norm(CE(passband),2)).';

    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% Calculate Rhos %%%
    %%%%%%%%%%%%%%%%%%%%%%
    
    if(mod(iter,2)) %used SPC

        Rho_AB_SPC(SPC_counter)=abs(gamma_A'*gamma_B)^2;
        Rho_AE_SPC(SPC_counter)=abs(gamma_A'*gamma_E)^2;

        
    else %no SPC
        Rho_AB(SPC_counter)=abs(gamma_A'*gamma_B)^2;
        Rho_AE(SPC_counter)=abs(gamma_A'*gamma_E)^2;
        
        SPC_counter=SPC_counter+1; %%%%%%TEMPORARY, later used below for Rho snr%%%%%%%%%%%%%%
    end
%     Rho_R(iter)=abs(cA_gains*cE_gains')^2;%/(norm(cA_gains,2)*norm(cE_gains,2)))^2;
    Rho_R(iter)=abs(cA*cE'/(norm(cA,2)*norm(cE,2)))^2;

    continue %%%%%%%%%%%%TEMPORARY %%%%%%%%%%%%%%

    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% Artificial Noise %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    % gamma_A=abs(randn(M,1))+1i*abs(randn(M,1));
    % gamma_B=sqrt(0.01/2/M)*(randn(M,1)+1i*randn(M,1))+gamma_A;
    % gamma_A=gamma_A/norm(gamma_A,2);
    % gamma_B=gamma_B/norm(gamma_B,2);
    % gamma_E=sqrt(10/2/M)*(abs(randn(M,1))+1i*abs(randn(M,1)));
    % gamma_E=gamma_E/norm(gamma_E,2);
    
    K=10000; %num symbols
    power_cap=100;
    phi=1/N;
    % phi=0.25;
    sigma_s=phi*power_cap;
    sigma_w=(1-phi)*N*power_cap/(N-1);
    
    SNRi_B=10;
    sigma_n=power_cap/SNRi_B;
    sigma_e=0; %0 noise added to eve
    v=zeros(K,N);
    w=zeros(K,N);
    x=zeros(K,N);
    y=zeros(K,N);
    z=zeros(K,N);
    
    s=sqrt(sigma_s/2)*bits2QPSK(randi([0 1],1,2*K));
    
    for k=1:K
        %%% generate noise vector 'v' and add it to the spreaded 's'
        w(k,:)=sqrt(sigma_w/2/N)*(randn(N,1)+1i*randn(N,1));
        v(k,:)=w(k,:)-((gamma_A'*w(k,:).')*gamma_A).';
        x(k,:)=gamma_A*s(k)+v(k,:).';
        
        P(k)=x(k,:)*x(k,:)'; %Power of signal
        o_w(k)=w(k,:)*w(k,:)'; %variance of noise
        o_s=P(k)-(N-1)/N*o_w(k); %calculated variance of data
        calc_phi(k)=1-o_w(k)*(N-1)/(N*P(k)); %phi for this iteration
        
        
        
        
        %%% Transmit with noise
        sigma_n=P(k)/SNRi_B;
        sigma_e=0; %no noise for eve
        y(k,:)=x(k,:)+sqrt(sigma_n/2/N)*(randn(1,N)+1i*randn(1,N));
        z(k,:)=x(k,:)+sqrt(sigma_e/2/N)*(randn(1,N)+1i*randn(1,N));
        
        
        
        
        
        %Despread
        s_b(k)=gamma_B'*y(k,:).';
        s_e(k)=gamma_E'*z(k,:).';
        
    end
    % mean(P)
    % mean(calc_phi)
    % mean(o_w)
%     rho_AB=abs(gamma_B'*gamma_A)
%     rho_AE=abs(gamma_E'*gamma_A)
    if mod(iter,2)
        SNRo_B(iter)=N*phi*Rho_AB_SPC(SPC_counter)*SNRi_B/(N/(N-1)*(1-phi)*(1-Rho_AB_SPC(SPC_counter))*SNRi_B+1)+1;
    else
       SNRo_B(iter)=N*phi*Rho_AB(SPC_counter)*SNRi_B/(N/(N-1)*(1-phi)*(1-Rho_AB(SPC_counter))*SNRi_B+1)+1; 
       SPC_counter=SPC_counter+1;
    end
    
    differences=sign(real(s))-sign(real(s_b).') | sign(imag(s))-sign(imag(s_b)).';
    B_errors(iter)=nnz(differences);
    
    differences=sign(real(s))-sign(real(s_e).') | sign(imag(s))-sign(imag(s_e)).';
    E_errors(iter)=nnz(differences);
    
end

%Plot CDF of each rho
figure(1);
h=cdfplot(Rho_AB_SPC);
h.Color = 'blue';
h.LineStyle = '--';
hold on;

h=cdfplot(Rho_AB);
h.Color = 'blue';
h.LineStyle = '-';

h=cdfplot(Rho_AE_SPC);
h.Color = 'red';
h.LineStyle = '--';

h=cdfplot(Rho_AE);
h.Color = 'red';
h.LineStyle = '-';

h=cdfplot(Rho_R);
h.Color = 'green';
h.LineStyle = '-';
hold off
legend('\rho_{AB} SPC','\rho_{AB}','\rho_{AE} SPC','\rho_{AE}','\rho_R')
title('Key Correlation');
xlabel('\rho')
ylabel('F(\rho)')


% tx_power=x(k,:)*x(k,:)';
% temp_sigma2_w=w(k,:)*w(k,:)';
%
% % phi=1/N;
%
%
% y(k,:)=x(k,:)+sigma_n*(randn(N,1)+1i*randn(N,1));
% z(k,:)=x(k,:)+sigma_e*(randn(N,1)+1i*randn(N,1));