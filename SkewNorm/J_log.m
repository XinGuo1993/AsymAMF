function result = J_log(A,B,C,J_sigma)
%UNTITLED8 �˴���ʾ�йش˺�����ժҪ
% A xָ���� B exp��xָ���� C exp��x^2ָ��
%%%%% ע���������logֵ��������ԭ����
%   �˴���ʾ��ϸ˵��
flag = 0;
if A>0.25
    flag =1;
    %% ȷ��mu0
    f=@(x)1./(x.*log(x))+A./x+B-2*C.*x;
    a = 0.001;
    b = 1-0.001;
    c = 1+0.001;
    d = 10;
    threshold = 0.0001;
    %ȷ�� ��ʼ�㣬���ں�������������ö��ַ��������׼ȷ��ֱ��ʹ��fzero�����������޷�����������һ��������
    fa=f(a);
    fb=f(b);
    fc=f(c);
    fd=f(d);
    while fa<0
        a = a/2;
        fa=f(a);
    end
    while fb>0
        b = (1+b)*0.5;
        fb=f(b);
    end
    while fc<0
        c = (1+c)*0.5;
        fc=f(c);
    end
    while fd>0
        d = 2*d;
        fd=f(d);
    end
    while min(abs(fa),abs(fb))>threshold
        tem = (a+b)*0.5;
        f_tem = f(tem);
        if f_tem>0
            a = tem;
            fa = f_tem;
        else
            b = tem;
            fb = f_tem; 
        end

    end

    while min(abs(fc),abs(fd))>threshold
        tem = (c+d)*0.5;
        f_tem = f(tem);
        if f_tem>0
            c = tem;
            fc = f_tem;
        else
            d = tem;
            fd = f_tem;
        end
    end

    if abs(fa) < abs(fb)
        tem1 = a;
    else
        tem1 = b;
    end
    if abs(fc) < abs(fd)
        tem2 = c;
    else
        tem2 = d;
    end
    % fprintf('erfen ending\n');

    %%
    mu1 = tem1;
    mu2 = tem2;
    sigma1 = 1/sqrt(1/mu1^2*(1/log(mu1)+1/log(mu1)^2+A) +2*C);
    sigma2 = 1/sqrt(1/mu2^2*(1/log(mu2)+1/log(mu2)^2+A) +2*C);
  if isreal(sigma1)&&isreal(sigma1)
        upper_bound1 =  (1 - mu1)/(sqrt(2)*sigma1);
        low_bound1 = (0 - mu1)/(sqrt(2)*sigma1);
        low_bound2 = (1 - mu2)/(sqrt(2)*sigma2);

        I1 = quadl(['f_log(x,',num2str(mu1),',',num2str(sigma1),',',num2str(A),',',num2str(B),',',num2str(C),')'],low_bound1,upper_bound1);
        I2 = quadToInf(['f_log(x,',num2str(mu2),',',num2str(sigma2),',',num2str(A),',',num2str(B),',',num2str(C),')'],low_bound2);

%         if isnan(I1)
%             I11 = NaN;
%         end
%         if isnan(I2)
%             I22 = NaN;
%         end

        if imag(I1) <0.01
            I11 = log(sqrt(2)*sigma1) + log(real(I1)) + log(-log(mu1))+A*log(mu1)+B*mu1-C*mu1^2;
        end
        if imag(I2) <0.01
            I22 = log(sqrt(2)*sigma2) + log(real(I2)) + log(log(mu2))+A*log(mu2)+B*mu2-C*mu2^2;
        end
        result = -exp(I11-J_sigma)+exp(I22-J_sigma);
  else
      flag =0;
  end
end
if flag ==0
    I22=quadToInf(['f_log2(x,',num2str(A),',',num2str(B),',',num2str(C),')']);
    result = I22/exp(J_sigma);
end


end

