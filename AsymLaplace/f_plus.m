function f_plus = f_plus(x,mu0,sigma0,A,B,C)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% function for J_plus
    f_plus =  exp( (B-2*C*mu0)*sqrt(2)*sigma0*x - 2*C*sigma0^2 *x.^2 + A*log(1+ sqrt(2)*sigma0/mu0 * x));
%      f_plus = (1+sqrt(2)*sigma0/mu0.*x).^A .*exp(sqrt(2)*B*sigma0.*x - 2*sigma0^2*C.*x.^2 -2*sqrt(2)*mu0*sigma0*C.*x);  
end

