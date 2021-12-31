clc;
% clear all;
k10=0.006;
ke0=0.0063;
k12=11.0;
k21=14.04;
k13=10.02;
k31=283.50;
Ce50=6.10;
gamma=1;
Bisrf=0.5;
%Bis=100 (maximum alert state), Bis=0 zero electrical activity,
%Bis_desired=0.45, NORMALIZED VALUE = 50/100 = 0.5 
%________________________________TRANSFER FUNCTION THAT MAPS U to C1_________________________________________________________________________
a11=(k10+k12+k13);
A=[-a11 k21 k31;
  k12 -k21 0 ;
  k13 0 -k31 ]
B=[1 0 0];
B = transpose(B);
C=[1 0 0];
D=0;
sys_State_sp = ss(A,B,C,D);  %state space model
A1=[-ke0]; B1=[ke0]; C1=[1]; sys_State_sp2 = ss(A1,B1,C1,D);
[Num,Den] = ss2tf(A,B,C,D);
transf=tf(Num,Den); %input is U, output is C_1

%C_1=U*transf

%________________________2nd transfer function that maps C1 to Ce__________
transf_ce=tf([ke0],[1,ke0]); %input is C_1 , output is C_e

%C_e=C_1*transf_ce

%________________________1st part of hill function__ Paper 2 method________
% transf_v =1/Ce50; %input is C_e , output is v
%v=C_e*transf_ce

psi = series(transf,transf_ce) %___linear map of u to effect site C_e;
% PSI MAPS DRUG TO EFFECT SITE LINEARLLY
%create v with pid control
Bisref=0.5;
Kp=0.86;
Kd=0.00;
Ki=0;
v=pid(Kp,Ki,Kd);
%____ Feedback linearization________
u_0=feedback_lin(v,gamma,Ce50); % write v in terms of u so it includes non linear cancelling terms
u=Bisref*u_0/psi;
bis=hill_function(u,gamma,psi,Ce50);
fdbck_loop=feedback(bis*u,1);

%                    ___ verify how systems evolves with no control____
t_vec=[0:0.1:1400];

for i=1:length(t_vec)
    U_vec(i)=0;
end  
init=[2 0 0];
y=lsim(sys_State_sp,U_vec,t_vec); 
yce=lsim(sys_State_sp2,y,t_vec,65090);
bisvec=[];
contr=[];
for i=1:length(t_vec)
    bisvec(i)=hill_0(1,yce(i),Ce50);
end 
step(fdbck_loop);
% hold on;
