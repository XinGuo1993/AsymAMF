function f_log2 = f_log2(x,A,B,C)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    f_log2 =  log(x).*exp(B*x - C*x.^2 + A*log(x));
end
