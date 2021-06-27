%Numerical example of Fedotenkov (2016, Econ Letters)

clear
clc
close all

%Calibration
alfa = 0.4;
rho = 0.4166;
n = 0;
tau = 0.20;
tau1 = 0.15;
SP_discount = 0.5;

%Model
k0 = 1;
k_e = k0;
k_ue = k0;
T =505;
%T = 5000 %for welfare analysis

for t=1:T
    
    theta_e = 1/ ( 1 + (1+rho)*(1+(1-alfa)*tau1/alfa) );
    theta_ue = 1/ ( 1 + (1+rho)*(1+(1-alfa)*tau/alfa) );
    
    if t==1
    cy_e(t) = (1-tau)*(1-alfa)*(1-theta_ue)*k0^alfa;
    co_e(t) = alfa*(theta_ue*(1-tau)*(1-alfa)*k0^alfa)^alfa*(1+n)^(1-alfa)*(1 +(1-alfa)*tau/alfa);
    Ue(t) = log(cy_e(t)) + (1/(1+rho))*log(co_e(t));

    cy_ue(t) = (1-tau)*(1-alfa)*(1-theta_ue)*k0^alfa;
    co_ue(t) = alfa*(theta_ue*(1-tau)*(1-alfa)*k0^alfa)^alfa*(1+n)^(1-alfa)*(1 +(1-alfa)*tau/alfa);
    Uue(t) = log(cy_ue(t)) + (1/(1+rho))*log(co_ue(t));
    end
    
    if t > 1
    k_e(t) = theta_ue*(1-alfa)*(1-tau)*k_e(t-1)^alfa;     
    cy_e(t) = (1-tau)*(1-alfa)*(1-theta_ue)*k_e(t-1)^alfa;
    co_e(t) = alfa*(theta_ue*(1-tau)*(1-alfa)*k_e(t-1)^alfa)^alfa*(1+n)^(1-alfa)*(1 +(1-alfa)*tau/alfa);
    Ue(t) = log(cy_e(t)) + (1/(1+rho))*log(co_e(t));
    c_old_e(t) = alfa*k_e(t)^alfa*(1+n)^(1-alfa)*(1 +(1-alfa)*tau/alfa);
    U_old_e(t) = log(cy_e(t-1)) + (1/(1+rho))*log(c_old_e(t));
    
    k_ue(t) = theta_ue*(1-alfa)*(1-tau)*k_ue(t-1)^alfa; 
    cy_ue(t) = (1-tau)*(1-alfa)*(1-theta_ue)*k_ue(t-1)^alfa;
    co_ue(t) = alfa*(theta_ue*(1-tau)*(1-alfa)*k_ue(t-1)^alfa)^alfa*(1+n)^(1-alfa)*(1 +(1-alfa)*tau/alfa);
    Uue(t) = log(cy_ue(t)) + (1/(1+rho))*log(co_ue(t));
    c_old_ue(t) = alfa*k_ue(t)^alfa*(1+n)^(1-alfa)*(1 +(1-alfa)*tau/alfa);
    U_old_ue(t) = log(cy_ue(t-1)) + (1/(1+rho))*log(c_old_ue(t));
    end
    
    if t == 499
    k_e(t) = theta_e*(1-alfa)*(1-tau)*k_e(t-1)^alfa;     
    cy_e(t) = (1-tau)*(1-alfa)*(1-theta_e)*k_e(t-1)^alfa;
    co_e(t) = alfa*(theta_e*(1-tau)*(1-alfa)*k_e(t-1)^alfa)^alfa*(1+n)^(1-alfa)*(1 +(1-alfa)*tau1/alfa);
    Ue(t) = log(cy_e(t)) + (1/(1+rho))*log(co_e(t));
    c_old_e(t) = alfa*k_e(t-1)^alfa*(1+n)^(1-alfa)*(1 +(1-alfa)*tau/alfa);
    U_old_e(t) = log(cy_e(t-1)) + (1/(1+rho))*log(c_old_e(t));

    k_ue(t) = theta_ue*(1-alfa)*(1-tau)*k_ue(t-1)^alfa; 
    cy_ue(t) = (1-tau)*(1-alfa)*(1-theta_ue)*k_ue(t-1)^alfa;
    co_ue(t) = alfa*(theta_ue*(1-tau)*(1-alfa)*k_ue(t-1)^alfa)^alfa*(1+n)^(1-alfa)*(1 +(1-alfa)*tau/alfa);
    Uue(t) = log(cy_ue(t)) + (1/(1+rho))*log(co_ue(t));
    k_old_ue(t) = k_ue(t-1);
    c_old_ue(t) = alfa*k_old_ue(t)^alfa*(1+n)^(1-alfa)*(1 +(1-alfa)*tau/alfa);
    U_old_ue(t) = log(cy_ue(t-1)) + (1/(1+rho))*log(c_old_ue(t));
    end
    
    if t >= 500
    k_e(t) = theta_e*(1-alfa)*(1-tau1)*k_e(t-1)^alfa;     
    cy_e(t) = (1-tau1)*(1-alfa)*(1-theta_e)*k_e(t-1)^alfa;
    co_e(t) = alfa*(theta_e*(1-tau1)*(1-alfa)*k_e(t-1)^alfa)^alfa*(1+n)^(1-alfa)*(1 +(1-alfa)*tau1/alfa);
    Ue(t) = log(cy_e(t)) + (1/(1+rho))*log(co_e(t));
    c_old_e(t) = alfa*k_e(t-1)^alfa*(1+n)^(1-alfa)*(1 +(1-alfa)*tau1/alfa);
    U_old_e(t) = log(cy_e(t-1)) + (1/(1+rho))*log(c_old_e(t));

    k_ue(t) = theta_e*(1-alfa)*(1-tau1)*k_ue(t-1)^alfa; 
    cy_ue(t) = (1-tau1)*(1-alfa)*(1-theta_e)*k_ue(t-1)^alfa;
    co_ue(t) = alfa*(theta_e*(1-tau1)*(1-alfa)*k_ue(t-1)^alfa)^alfa*(1+n)^(1-alfa)*(1 +(1-alfa)*tau1/alfa);
    Uue(t) = log(cy_ue(t)) + (1/(1+rho))*log(co_ue(t));
    k_old_ue(t) = k_ue(t-1);
    c_old_ue(t) = alfa*k_old_ue(t)^alfa*(1+n)^(1-alfa)*(1 +(1-alfa)*tau1/alfa);
    U_old_ue(t) = log(cy_ue(t-1)) + (1/(1+rho))*log(c_old_ue(t));
    end

    X(t) = t-499;
    
end

for t=100:T
    %EV calcs for cy (so no need to multiply by (1+rho)/(2+rho))
    Lambda(t) = 100*( exp((Ue(t)-Ue(498))) -1);
    Lambda1(t) = 100*( exp((Uue(t)-Uue(498))) -1);
    Lambda2(t) = 100*( exp((U_old_e(t)-U_old_e(498))) -1);
    Lambda3(t) = 100*( exp((U_old_ue(t)-U_old_ue(498))) -1);
    
end

for t=499:T
 
 if t==499
  Welfare_e1(t) = SP_discount^(-1)*U_old_e(t) + Ue(t); 
  Welfare_ue1(t) = SP_discount^(-1)*U_old_ue(t) + Uue(t);
 end
 
 if t>499
 Welfare_e1(t) = SP_discount^(t-499)*Ue(t);
 Welfare_ue1(t) = SP_discount^(t-499)*Uue(t);
 end
 
end

for t=500:T
 
    if t==500
 Welfare_e11(t) = SP_discount^(-1)*U_old_e(t) + Ue(t);
 Welfare_ue11(t) = SP_discount^(-1)*U_old_ue(t) + Uue(t);
    end
 
 if t>500
  Welfare_e11(t) = SP_discount^(t-500)*Ue(t);
 Welfare_ue11(t) = SP_discount^(t-500)*Uue(t);
 end
 
end

Welfare_e = sum(Welfare_e1);
Welfare_ue = sum(Welfare_ue1);

%Welfare_e = sum(Ue(499:T));
%Welfare_ue = sum(Uue(499:T));

%EV0
Lambda4 = 100*( exp( (Welfare_ue - Welfare_e)*(1-SP_discount)/(1-SP_discount^(T+1)) ) -1 )

Welfare_e1 = sum(Welfare_e11);
Welfare_ue1 = sum(Welfare_ue11);

%Welfare_e1 = U_old_e(500) + sum(Ue(500:T));
%Welfare_ue1 = U_old_ue(500) + sum(Uue(500:T));

%EV1
Lambda5 = 100*( exp( (Welfare_ue1 - Welfare_e1)*(1-SP_discount)/(1-SP_discount^(T+1)) ) -1 )

%figure(1)
%hold on, plot(X(497:end), Ue(497:end)), hold on, plot(X(497:end), Uue(497:end),'r')

%figure(2)
%hold on, plot(X(497:end), U_old_e(497:end)), hold on, plot(X(497:end), U_old_ue(497:end), 'r')

%figure(3)
%hold on, plot(X(497:end), Lambda(497:end)), hold on, plot(X(497:end), Lambda1(497:end),'r')

%figure(4)
%hold on, plot(X(497:end), Lambda2(497:end)), hold on, plot(X(497:end), Lambda3(497:end),'r')

figure(1)
hold on, subplot(2,1,1),  plot(X(497:end), cy_e(497:end)/cy_e(400)), hold on, plot(X(497:end), cy_ue(497:end)/cy_e(400),'r'), 
%hold on, subplot(2,2,2), plot(X(497:end), Lambda(497:end)/100), hold on, plot(X(497:end), Lambda1(497:end)/100,'r')
xlabel('Time')
title('Young-age consumption: Announced (blue), Unnannounced (red)') 
hold on, subplot(2,1,2), plot(X(497:end), c_old_e(497:end)/c_old_e(400)), hold on, plot(X(497:end), c_old_ue(497:end)/c_old_e(400),'r')
%hold on, subplot(2,2,4), plot(X(497:end), Lambda2(497:end)/100), hold on, plot(X(497:end), Lambda3(497:end)/100,'r')
xlabel('Time')
title('Realized old-age consumption: Announced (blue), Unnannounced (red)') 