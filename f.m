function f = f(x,mu0,sigma0,A,B,C)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    f =  exp( (B-2*C*mu0)*sqrt(2)*sigma0*x - 2*C*sigma0^2 *x.^2 + A*log(1+ sqrt(2)*sigma0/mu0 * x));
end

