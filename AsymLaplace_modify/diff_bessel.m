function result = diff_bessel(mu,z)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
  h = 0.00001;
%   result = (log(besselk(mu+h,z)) - log(besselk(mu-h,z)))/(2*h);
  result = log( div_bessel(mu-2*h,z,2*h) )/(2*h);
end

