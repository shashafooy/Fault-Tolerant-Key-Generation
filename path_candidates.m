function [alpha,k] = path_candidates(c_in,M,long_p)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

alpha=zeros(1,M);
k=zeros(1,M);

% Nc=length(c_in);
% z=zeros(1,Nc-length(p)+1);
% z(ceil(length(z)/2))=1;
% 
% long_p=conv(z,p);
% %     long_p=conv(z,srrc);
% long_p=circshift(long_p,-Nc/2+2);

p_L=circshift(long_p,324/2);
p_L=p_L(1:324);

c(1,:)=c_in;
%Get M path candidates -  k and alpha
for i=1:M

    [t,lags]=xcorr(c(i,:),p_L);
%     t=t(1:length(long_p));
    [~,idx]=max(t);
    k(i)=lags(idx)+length(p_L)/2+1;

%     k(i)=k(i)-5;
    
%     temp=long_p;
    %cross correlation

%     t=zeros(1,length(c_in));
%     for j=1:length(c_in)/2
%         t(j)=(c(i,:))*circshift(long_p,j-1)';
%     end
%     [~,k(i)]=max(abs(t));


%     [~,idx]=max(abs(t));
%     
%     oldmin=10000;
%     for j=idx-100:idx+100
%         newmin=norm(c(i,:)-c(i,j)*circshift(long_p,j-1));
%         if(newmin<oldmin)
%             oldmin=newmin;
%             best=j;
%         end
%     end
%     k(i)=best;
    
    alpha(i)=c(i,k(i));
    c(i+1,:)=c(i,:)-alpha(i)*circshift(long_p,k(i)-1);
    
%     figure(1),plot(abs(c(i,:))),hold on,plot(abs(abs(alpha(i))*circshift(long_p,k(i)-1))),plot(abs(c(i+1,:))),hold off,shg
% %     figure(2),plot(abs(c(i+1,:)))
%     pause(0.5)
    
    
    
%     [~, idx]=maxk(abs(c(i,:)),M*4); %Get indicies of M items %%%%%%%%%%%%USED TO BE 4%%%%%%%%%%%%%%%%%
%     oldmin=norm(c(i,:),2)^2;
%     %             for j=1:length(c(i,:))
%     for j=idx
%         %           for j=1:M
%         %                 idx=mod(j + idx_inc-2,length(cA_L))+1;
%         %                 newmin=norm(c(i,:)-c(i,idx)*circshift(long_p,j-1),2)^2;
%         newmin=norm(c(i,:)-c(i,j)*circshift(long_p,j-1),2)^2;
%         %                 newmin=norm(c(i,:)-c(i,idx(j))*circshift(long_p,idx(j)-1),2)^2;
%         if(newmin<oldmin)
%             
%             %                     best_a=c(i,idx);
%             best_a=c(i,j);
%             best_idx=j;
%             
%             %                     best_a = c(i,idx(j));
%             %                     best_idx = idx(j);
%             oldmin=newmin;
%         end
%     end
%     %save alpha and k that give the minimum
%     %             [~,idx]=max(abs(c(i,:)))
%     %             best_idx
%     alpha(i)=best_a;
%     k(i)=best_idx;
%     c(i+1,:)=c(i,:)-alpha(i)*circshift(long_p,k(i)-1);
    %                     figure(2);plot(abs(c(i+1,:)));title('c(i+1)')
    %                     figure(3);plot(abs(c(i,:)));hold on;plot(abs(best_a*circshift(long_p,best_idx)));hold off;title('c and strongest path');shg
    %                     figure(1);plot(abs(c(i,:)));hold on;plot(abs(c(i,:)-best_a*circshift(long_p,best_idx)));hold off;title('removed strongest path');shg
    %                     figure(4);plot(abs(cA_L));
    %             figure(1),plot(abs(c(i,:))),hold on,plot(abs(alpha(i))*circshift(long_p,best_idx-2)),hold off,shg
    %             figure(2),plot(abs(c(i+1,:)))
end

end

%     for num=1:3
%     for num=1:3
%         clear c
%         if(num==1)
%             c(1,:)=cA_L;
%         elseif(num==2)
%             c(1,:)=cB_L;
%         else
%             c(1,:)=cE_L;
%         end
%         %Get M path candidates -  k and alpha
%         for i=1:M
%             [~, idx]=maxk(abs(c(i,:)),M*L); %Get indicies of M items
%             oldmin=norm(c(i,:),2)^2;
% %             for j=1:length(c(i,:))
%             for j=idx
% %           for j=1:M
% %                 idx=mod(j + idx_inc-2,length(cA_L))+1;
% %                 newmin=norm(c(i,:)-c(i,idx)*circshift(long_p,j-1),2)^2;
%                 newmin=norm(c(i,:)-c(i,j)*circshift(long_p,j-1),2)^2;
% %                 newmin=norm(c(i,:)-c(i,idx(j))*circshift(long_p,idx(j)-1),2)^2;
%                 if(newmin<oldmin)
%                     
% %                     best_a=c(i,idx);
%                     best_a=c(i,j);
%                     best_idx=j;
%                     
% %                     best_a = c(i,idx(j));
% %                     best_idx = idx(j);
%                     oldmin=newmin;
%                 end
%             end
%             %save alpha and k that give the minimum
% %             [~,idx]=max(abs(c(i,:)))
% %             best_idx
%             alpha(i)=best_a;
%             k(i)=best_idx;
%             c(i+1,:)=c(i,:)-alpha(i)*circshift(long_p,k(i)-1);
% %                     figure(2);plot(abs(c(i+1,:)));title('c(i+1)')
% %                     figure(3);plot(abs(c(i,:)));hold on;plot(abs(best_a*circshift(long_p,best_idx)));hold off;title('c and strongest path');shg
% %                     figure(1);plot(abs(c(i,:)));hold on;plot(abs(c(i,:)-best_a*circshift(long_p,best_idx)));hold off;title('removed strongest path');shg
% %                     figure(4);plot(abs(cA_L));
% %             figure(1),plot(abs(c(i,:))),hold on,plot(abs(alpha(i))*circshift(long_p,best_idx-2)),hold off,shg
% %             figure(2),plot(abs(c(i+1,:)))
%         end
%         if(num==1)
%             alpha_A=alpha;
%             k_A=k;
%         elseif(num==2)
%             alpha_B=alpha;
%             k_B=k;
%         else
%             alpha_E=alpha;
%             k_E=k;
%         end
%         %     [k(i) alpha(i)]=argmin(norm(c(1)-alpha*p(n-k),2));
%         %     c(i+1)=c(1)-alpha(i)*p(n-k(i));
%         
%     end

