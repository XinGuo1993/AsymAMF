function I = J_plus(A,B,C)
%UNTITLED8 �˴���ʾ�йش˺�����ժҪ
% A xָ���� B exp��xָ���� C exp��x^2ָ��
%%%%% ע���������logֵ��������ԭ����
%   �˴���ʾ��ϸ˵��
if A>0
mu0 = (B+sqrt(B^2 +8*A*C))/(4*C);
sigma0 = (B+sqrt(B^2 +8*A*C))/sqrt(4*C*(8*A*C+B^2+B*sqrt(B^2 +8*A*C) ));
low_bound = -mu0/(sqrt(2)*sigma0); 
%quadl('f(x,0.6076,0.3389,1,2,3)',1,2)
%quad(['f(x,',num2str(mu0),',',num2str(sigma0),',',num2str(A),',',num2str(B),',',num2str(C),')'],-1,1)
I0 = quadToInf(['f_plus(x,',num2str(mu0),',',num2str(sigma0),',',num2str(A),',',num2str(B),',',num2str(C),')'],low_bound);
% if isreal(I0)
%     I = sqrt(2)* real(I0) * sigma0 * mu0^A*exp(B*mu0-C*mu0^2);
if imag(I0)<0.001
     I = log(sqrt(2)*sigma0) + log(real(I0)) + A*log(mu0)+B*mu0-C*mu0^2;
end

else     
     I = log(quadToInf(['f_plus2(x,',num2str(A),',',num2str(B),',',num2str(C),')']));
end
end

