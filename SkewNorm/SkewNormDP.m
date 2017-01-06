function [opts,Gamma, Phi, E_inv_sigma,E_s] = SkewNormDP(Y, U, V, Cov_U, Cov_V,Gamma,Phi,E_inv_sigma,E_s,Mask,opts)


%%Initialize 
[M, N] = size(Y);

lamda_s     = opts.lamda_s;
alpha0   = opts.alpha0;
a0       = opts.a0;
b0       = opts.b0;
s_mu        = opts.s_mu;
s_sigma     = opts.s_sigma;
K        = opts.K;
R        = opts.r;       % ��ʼrank
Pai      = zeros(1,K);
err      = zeros(M,N);   % E(y-uv) = err
err2     = zeros(M,N);   % E(y-uv)^2 = err2
temp0    = (U*V').*Y;

E_s2      = ones(K,1);   % ��E��s2����ʼ������Ϊһ��ʼ�ڼ���sigma�ķֲ���ʱ��s�Ĺ��ƻ���δ֪,���������ʼ����

% ����err, err2
err = Y - U*V'; 
N_1 = reshape(V',R,1,N);
N_2 = reshape(V',1,R,N);
temp_N=  bsxfun(@times,N_1,N_2 )+Cov_V ;

M_1 = reshape(U',R,1,M);
M_2 = reshape(U',1,R,M);
temp_M=  bsxfun(@times,M_1,M_2 )+Cov_U ;

MatU =reshape(temp_M, R*R,M);
MatV =reshape(temp_N, R*R,N);
err2= Y.^2+MatU'*MatV-2*temp0;   %  size M*N
clear temp0 temp_M temp_N;
clear MatU;
clear MatV;
for i =1:K
    phi(i) =  sum(sum(Phi(:,:,i)));     % ��M,N���
end
Sum_to_k = sum(phi);    % ��K���
%% Update
% try  

for iter = 1 : 1
%     for k=1:K
%          Gamma(k,1) = 1 + phi(k);
%         Sum_to_k   = Sum_to_k - phi(k);
%         Gamma(k,2) = alpha0+Sum_to_k;
%         psi_vec(k) = psi(Gamma(k,2))-psi(Gamma(k,1)+Gamma(k,2));
%         pai_vec(k) = Gamma(k,2)/(Gamma(k,1)+Gamma(k,2));
%     end

    for x=1:K
        Gamma(x,1) =1+phi(x);   
        Sum_to_k = Sum_to_k -phi(x);
        Gamma(x,2) = alpha0+Sum_to_k;
        if(x==1)
            psi_vec_sum(x) =0;
            pai_vec_prod(x)=1;
        elseif(x==2)
            psi_vec_sum(x) =psi(Gamma(x-1,2))-psi(Gamma(x-1,1)+Gamma(x-1,2));
            pai_vec_prod(x)=Gamma(x-1,2)/(Gamma(x-1,1)+Gamma(x-1,2));
        else
            psi_vec_sum(x) = psi_vec_sum(x-1)+psi(Gamma(x-1,2))-psi(Gamma(x-1,1)+Gamma(x-1,2));
            pai_vec_prod(x) = pai_vec_prod(x-1)*Gamma(x-1,2)/(Gamma(x-1,1)+Gamma(x-1,2));
        end
    end

    for k = 1 : K
        Sum_of_psi =psi_vec_sum(k);  
        Prod_of_pai =pai_vec_prod(k);
%         Tao(k,1)   = a0 + 0.5*phi(k);                           % change guoxin
%         Tao(k,2)   = b0 + 0.5*sum(sum(Phi(:,:,k).* err2));       % change guoxin
         A_sigma = 2*a0+2+phi(k);
         B_sigma = lamda_s*sqrt(1+lamda_s^2)* E_s(k) * sum(sum(Phi(:,:,k).* err));
         C_sigma = b0 + 0.5*(1+lamda_s^2)* sum(sum(Phi(:,:,k).* err2));
         J_sigma = J_plus(A_sigma-3,B_sigma,C_sigma);
         E_inv_sigma(k,1) = exp(J_plus(A_sigma+1-3,B_sigma,C_sigma)-J_sigma);
         E_inv_sigma(k,2) = exp(J_plus(A_sigma+2-3,B_sigma,C_sigma)-J_sigma);
         E_log_sigma = J_log(A_sigma-3,B_sigma,C_sigma,J_sigma);
%          [tem0,tem1] = J_log(A_sigma-3,B_sigma,C_sigma);
%          if isnan(tem1)
%              E_log_sigma = -exp(tem0-J_sigma);
%          elseif tem0 =='no need'
%              E_log_sigma = tem1/exp(J_sigma);
%          else
%              E_log_sigma = -exp(tem0-J_sigma)+exp(tem1-J_sigma);
%          end
          %��s���й���       
         A_s = 1 + phi(k)*lamda_s^2;
         B_s = lamda_s*sqrt(1+lamda_s^2)*E_inv_sigma(k,1)*sum(sum(Phi(:,:,k).* err));
         E_s(k) = B_s/A_s;
         E_s2(k) = E_s(k)^2 + 1/A_s;
         
        if k == K
            Gamma(K,2) = 0;
        end
        
        Phi(:,:,k) = exp((psi(Gamma(k,1)) - psi(Gamma(k,1)+Gamma(k,2))+Sum_of_psi) *ones(M,N) -0.5*lamda_s^2* E_s2(k)^2*ones(M,N) + E_log_sigma*ones(M,N) - 0.5*(1+lamda_s^2)*E_inv_sigma(k,2)*err2 + lamda_s*sqrt(1+lamda_s^2)*E_s(k)*E_inv_sigma(k,1)*err);   %%change guoxin���˴���Ҫ�����޸�
        Pai(k)     = Prod_of_pai*Gamma(k,1)/(Gamma(k,1)+Gamma(k,2));

    end 
 
%      if opts.verbose
%          disp(Pai);
%      end
    % catch 
     Pai   = Pai/sum(Pai);
     flag  = (Pai>1e-4);
     Pai            = Pai(flag);
     K              = size(Pai,2);
     opts.K         = K;
     E_s = E_s(flag);
     E_s2 = E_s2(flag);
     Gamma          = Gamma(flag,:);
     E_inv_sigma    = E_inv_sigma(flag,:);
     Phi            = Phi(:,:,flag);
     Phi = bsxfun(@rdivide, Phi,sum(Phi,3));
     Phi = bsxfun(@times, Phi, Mask);
end

end
