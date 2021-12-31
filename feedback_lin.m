function u= feedback_lin(v,gamma,E50)
a=v*((1-v)^-1);
% x = sprintf("this is %d ",a);
% disp(x);
%fprintf("this is %f \n",a);
u=(a^(1/gamma))*E50;
end