function f = f(x,mu0,sigma0,A,B,C)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    f =  exp( (B-2*C*mu0)*sqrt(2)*sigma0*x - 2*C*sigma0^2 *x.^2 + A*log(1+ sqrt(2)*sigma0/mu0 * x));
end

