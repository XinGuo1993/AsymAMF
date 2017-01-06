%%测试 J_plus 和J_log的稳定性
A_sigma = 10;
B_sigma = 0;
C_sigma = 10;
% for A_sigma = 5:100
%     for C_sigma = 5:100
%         J_sigma = J_plus(A_sigma-3,0,C_sigma);
%         A1 = exp(J_plus(A_sigma+2-3,0,C_sigma)-J_sigma);
%         A2 = (A_sigma-2)/(2.0*C_sigma);
%         (A1-A2)<0.0001
%     end
% end
 for A_sigma = 1000:2000
     for C_sigma = 1000:2000
        J_sigma = J_plus(A_sigma-3,0,C_sigma);
        A3 = -0.5*(log(C_sigma)-psi((A_sigma-2)/2.0));
        A4 = J_log(A_sigma-3,0,C_sigma,J_sigma);
        (A1-A2)<0.0000001
     end
 end

% I22=quadToInf(['f_log2(x,',num2str(A_sigma-3),',',num2str(B_sigma),',',num2str(C_sigma),')'])
% I222=quadl(['f_log2(x,',num2str(A_sigma-3),',',num2str(B_sigma),',',num2str(C_sigma),')'],0,1)
% I23=quadToInf(['f_log2(x,',num2str(A_sigma-3),',',num2str(B_sigma),',',num2str(C_sigma),')'],1)
% A5 = I22/exp(J_sigma)
% A6 = J_log(A_sigma-3,0,C_sigma,1)
% A7 = I222/exp(J_sigma)