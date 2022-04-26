load('keys.mat');
iterations=length(store_gamma_A(:,1));

for iter=1:iterations
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