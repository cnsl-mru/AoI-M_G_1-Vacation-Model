%% Initialization
clear ; close all; clc

%% ==================== Part 1: Basic Function ====================
% Complete warmUpExercise.m
fprintf('Average AoI Running ... \n');


fprintf('Average AoI: \n');

fprintf('Program paused. Press enter to continue.\n');
pause;
%% ======================= Part 2: Base Line Model with respect lambda1=======================

fprintf('Plotting Data ...\n')
p=[];
lambda=[];
lambda
rho=[];
AoI=[];
AoI1=[];

ES=0.5% expected service time
ES2=2*(ES)^2  %2nd order moment of service time
for M=2:5
        for i=2:M
            lambda(i)=0.12;
         end
        l1=sum(lambda);
       
        m1=1
        for k=0.1:0.02:1.4
            lambda(1)=k;
              for i=1:M
                 rho(i)=lambda(i)*ES;
             end
                 r1=sum(rho);
            AoI11=(lambda(1))*(1-r1)/(1/r1-rho(1))+ES*((1/rho(1))+(r1/(1-r1))+((2*rho(2)-1)*(r1-1))/(1-rho(2))^2+(2*rho(1)*rho(2)*(r1-1))/(1-rho(2))^3)
            AoI(m1)=AoI11;
            m1=m1+1
        end
        x=(0.1:0.02:1.4)
        p(M)=plot(x,AoI,'*--')
        hold on;
end
legend([p(2),p(3),p(4),p(5)],'M=2','M=3','M=4','M=5')
%% ======================= Part3: Base Line Model with respect  lambda1=======================

fprintf('Plotting Data ...\n')
p=[];
lambda=[];
rho=[];
AoI=[];
AoI1=[];
ES=0.5% expected service time
ES2=2*(ES)^2  %2nd order moment of service time
for M=2:5
% %Exponential
%ES=1;% expected service time

%Erlang-2
% ES=0.5% expected service time
% ES2=1*(ES)^2  %2nd order moment of service time
% EV=1/2 %expected vacation time
% EV2=1*EV^2 % second order moment of vacation
%HyperExponential with p
% p=0.4;
% ES=p*3/4+(1-p)*1/(2^2)% expected service time
% ES2=p*2*3^2/(4^2)+(1-p)*2/(2^4)  %2nd order moment of service time
m1=1;
for k=0.1:0.02:1.4
   
     lambda(1)=k;
        for i=2:M
        lambda(i)=0.1;
         end
        l1=sum(lambda);
         for i=1:M
        rho(i)=lambda(i)*ES;
        end
        rho(1)
        r1=sum(rho);
ES2=2*ES^2;  %2nd order moment of service time
EW=l1*ES2/(2*(1-r1)); % Expected waiting Time
LS1=(1/ES)/((1/ES)+lambda(1));    %Exponential S*(lambda(1)))
LS2=-(1/ES)/((1/ES)+lambda(1))^2; %Exponential S*'(lambda(1))) 
LS3=(2/ES)/((1/ES)+lambda(1))^3;  %Exponential S*''(lambda(1))) 
WS1=((1-r1)*lambda(1)*LS1)/((lambda(1)-l1*(1-LS1)))
WS2=(1-r1)*(l1*LS1^2+(lambda(1)^2-lambda(1)*l1)*LS2-l1*LS1)/(lambda(1)-l1*(1-LS1))^2
WS3=0;
for k1=2:M
LS11=(1/ES)/((1/ES)+lambda(k1));    %Exponential S*(lambda(1)))
LS21=-(1/ES)/((1/ES)+lambda(k1))^2; %Exponential S*'(lambda(1))) 
LS31=(2/ES)/((1/ES)+lambda(k1))^3;  %Exponential S*''(lambda(1))) 
WS11=((1-r1)*lambda(1)*LS11)/((lambda(1)-l1*(1-LS11)))
WS21=(1-r1)*(l1*LS11^2+(lambda(1)^2-lambda(1)*l1)*LS21-l1*LS11)/(lambda(1)-l1*(1-LS11))^2;
WS3=WS3+lambda(k1)*((2/lambda(1))+WS21-(2.0*WS11/lambda(1))) 
end;
AoI(m1)=(EW+2.0*ES+(2.0*WS1/lambda(1))-WS2-(1/lambda(1))+ES*WS3)
m1=m1+1;

end;
x=(0.1:0.02:1.4)
p(M)=plot(x,AoI,'*--')
hold on;
end;
legend([p(2),p(3),p(4),p(5)],'M=2','M=3','M=4','M=5')
%% ======================= Part 4: Plotting with respect  lambda1 for vacation model=======================
fprintf('Plotting Data ...\n')
 
for M=2:5% number of sources
    lambda=[];
    rho=[];
    AoI=[];
    AoI1=[];
   % %Exponential
    ES=0.4% expected service time
    ES2=2*(ES)^2  %2nd order moment of service time
    EV=2%expected vacation time
    EV2=2*ES^2 % second order moment of vacation
    %Erlang-2
    % ES=0.4% expected service time
    % ES2=1*(ES)^2  %2nd order moment of service time
    % EV=1/2 %expected vacation time
    % EV2=1*EV^2 % second order moment of vacation
    %HyperExponential with p
    % p=0.4;
    % ES=p*3/4+(1-p)*1/(2^2)% expected service time
    % ES2=p*2*3^2/(4^2)+(1-p)*2/(2^4)  %2nd order moment of service time
    % EV=1/2 %expected vacation time
    % EV2=2*EV^2 % second order moment of vacation
    
    m1=1;
    for k=0.1:0.02:1.2
    
        lambda(1)=k;
        for i=2:M
        lambda(i)=0.4;
         end
        l1=sum(lambda);
         for i=1:M
        rho(i)=lambda(i)*ES;
        end
        rho(1)
        r1=sum(rho);

    EW=(l1*ES2/2*(1-r1))+(EV2)/(2*EV); % Expected waiting Time
    LS1=(1/ES)/((1/ES)+lambda(1));    %Exponential S*(lambda(1)))
    LS2=-(1/ES)/((1/ES)+lambda(1))^2; %Exponential S*'(lambda(1))) 
    LS3=(2/ES)/((1/ES)+lambda(1))^3;  %Exponential S*''(lambda(1))) 
    LV1=(1/EV)/((1/EV)+lambda(1));    %%Exponential V*(lambda(1)))
    LV2=-(1/EV)/((1/EV)+lambda(1))^2;  %Exponential V*'(lambda(1)))
    LV3=(2/EV)/((1/EV)+lambda(1))^3;  %Exponential V*''(lambda(1)))
    WS1=((1-r1)*LS1)*(1-LV1)/((lambda(1)-l1*(1-LS1))*EV)
    WS2=((1-r1)/EV)*((lambda(1)-l1*(1-LS1))*(LS2*(1-LV1)-LS1*LV2)-(LS1*(1-LV1)*(1+l1*LS2)))/(lambda(1)-l1*(1-LS1))^2;
    WS3=0;
    for k1=2:M
        LS11=(1/ES)/((1/ES)+lambda(k1));    %Exponential S*(lambda(1)))
        LS21=-(1/ES)/((1/ES)+lambda(k1))^2; %Exponential S*'(lambda(1))) 
        LS31=(2/ES)/((1/ES)+lambda(k1))^3;  %Exponential S*''(lambda(1))) 
        LV11=(1/EV)/((1/EV)+lambda(k1));    %%Exponential V*(lambda(1)))
        LV21=-(1/EV)/((1/EV)+lambda(k1))^2;  %Exponential V*'(lambda(1)))
        LV31=(2/EV)/((1/EV)+lambda(k1))^3;  %Exponential V*''(lambda(1)))
        WS11=((1-r1)*LS11)*(1-LV11)/((lambda(k1)-l1*(1-LS11))*EV)
        WS21=((1-r1)/EV)*((lambda(k1)-l1*(1-LS11))*(LS21*(1-LV11)-LS11*LV21)-(LS11*(1-LV11)*(1+l1*LS21)))/(lambda(k1)-l1*(1-LS11))^2;
        WS3=WS3+lambda(k1)*(2/lambda(1)+WS21-2.0*WS11/lambda(1)) 
     end;
    AoI(m1)=EW+2.0*(ES+EV)+(1/lambda(1))*2.0*WS1-WS2-(1/lambda(1))+(ES+EV)*WS3
    m1=m1+1;
    
    end;
     x=(0.1:0.02:1.2)
%     plot(x,AoI,'-*');
     p(M)=plot(x,AoI,'*-');
     hold on;
end;
legend([p(2),p(3),p(4),p(5)],'M=2','M=3','M=4','M=5')
   
%% ======================= Part5: Base Line Model with respect  ES =======================

fprintf('Plotting Data ...\n')
p=[];
lambda=[];
rho=[];
AoI=[];
AoI1=[];
ES=0.5% expected service time
ES2=2*(ES)^2  %2nd order moment of service time
for M=2:5
% %Exponential
%ES=1;% expected service time

%Erlang-2
% ES=0.2% expected service time
% ES2=1*(ES)^2  %2nd order moment of service time
% EV=1/2 %expected vacation time
% EV2=1*EV^2 % second order moment of vacation
%HyperExponential with p
% p=0.4;
% ES=p*1/3+(1-p)*1/(3^2)% expected service time
% ES2=p*2/(3^2)+(1-p)*2/(3^4)  %2nd order moment of service time
% EV=1/2 %expected vacation time
% EV2=2*EV^2 % second order moment of vacation
m1=1;lambda(1)=1;
for k=0.1:0.02:0.6
   
        ES=k% expected service time
        ES2=2*(ES)^2  %2nd order moment of service time
%Erlang-2
% ES=k% expected service time
% ES2=1*(ES)^2  %2nd order moment of service time
%HyperExponential with p
% p=0.4;
% ES=p*k+(1-p)*k^2% expected service time
% ES2=p*2*ES^2+(1-p)*2*ES^2  %2nd order moment of service time
        for i=2:M
        lambda(i)=0.12;
         end
        l1=sum(lambda);
         for i=1:M
        rho(i)=lambda(i)*ES;
        end
        rho(1)
        r1=sum(rho);
ES2=2*ES^2;  %2nd order moment of service time
EW=l1*ES2/(2*(1-r1)); % Expected waiting Time
LS1=(1/ES)/((1/ES)+lambda(1));    %Exponential S*(lambda(1)))
LS2=-(1/ES)/((1/ES)+lambda(1))^2; %Exponential S*'(lambda(1))) 
LS3=(2/ES)/((1/ES)+lambda(1))^3;  %Exponential S*''(lambda(1))) 
WS1=((1-r1)*lambda(1)*LS1)/((lambda(1)-l1*(1-LS1)))
WS2=(1-r1)*(l1*LS1^2+(lambda(1)^2-lambda(1)*l1)*LS2-l1*LS1)/(lambda(1)-l1*(1-LS1))^2
WS3=0;
for k1=2:M
LS11=(1/ES)/((1/ES)+lambda(k1));    %Exponential S*(lambda(1)))
LS21=-(1/ES)/((1/ES)+lambda(k1))^2; %Exponential S*'(lambda(1))) 
LS31=(2/ES)/((1/ES)+lambda(k1))^3;  %Exponential S*''(lambda(1))) 
WS11=((1-r1)*lambda(1)*LS11)/((lambda(1)-l1*(1-LS11)))
WS21=(1-r1)*(l1*LS11^2+(lambda(1)^2-lambda(1)*l1)*LS21-l1*LS11)/(lambda(1)-l1*(1-LS11))^2;
WS3=WS3+lambda(k1)*((2/lambda(1))+WS21-(2.0*WS11/lambda(1))) 
end;
AoI(m1)=(EW+2.0*ES+(2.0*WS1/lambda(1))-WS2-(1/lambda(1))+ES*WS3)
m1=m1+1;

end;
x=(0.1:0.02:0.6)
p(M)=plot(x,AoI,'*--')
hold on;
end;
legend([p(2),p(3),p(4),p(5)],'M=2','M=3','M=4','M=5')




%% ======================= Part 6: Plotting with respect to ES for Vacation Model =======================
fprintf('Plotting Data ...\n')
% number of sources
for M=2:5
   lambda=[];
    rho=[];
    AoI=[];
    AoI1=[];
% %Exponential
    ES=0.2% expected service time
    ES2=2*(ES)^2  %2nd order moment of    service time
% EV=1/2 %expected vacation time
% EV2=1/2 % second order moment of vacation
%Erlang-2
% ES=0.2% expected service time
% ES2=1*(ES)^2  %2nd order moment of service time
% EV=1/2 %expected vacation time
% EV2=1*EV^2 % second order moment of vacation
%HyperExponential with p
% p=0.4;
% ES=p*1/3+(1-p)*1/(3^2)% expected service time
% ES2=p*2/(3^2)+(1-p)*2/(3^4)  %2nd order moment of service time
% EV=1/2 %expected vacation time
% EV2=2*EV^2 % second order moment of vacation
    lambda(1)=0.2;
    m1=1;EV=0.5;EV2=2*EV^2;
    for k=0.1:0.02:0.6
        %%Exponential
        ES=k;
        ES2=2*ES^2;
        %Erlang-2
    % ES=k% expected service time
    % ES2=1*(ES)^2  %2nd order moment of service time
    %HyperExponential with p
    % p=0.4;
    % ES=p*k+(1-p)*k^2% expected service time
    % ES2=p*2*k^2+(1-p)*2*k^4  %2nd order moment of service time
    for i=2:M
        lambda(i)=0.12;
    end
     l1=sum(lambda);
    for i=1:M
        rho(i)=lambda(i)*ES/M;
    end
    rho(1)
    r1=sum(rho);

    EW=(l1*ES2/2*(1-r1))+(EV2)/(2*EV); % Expected waiting Time
    LS1=(1/ES)/((1/ES)+lambda(1));    %Exponential S*(lambda(1)))
    LS2=-(1/ES)/((1/ES)+lambda(1))^2; %Exponential S*'(lambda(1))) 
    LS3=(2/ES)/((1/ES)+lambda(1))^3;  %Exponential S*''(lambda(1))) 
    LV1=(1/EV)/((1/EV)+lambda(1));    %%Exponential V*(lambda(1)))
    LV2=-(1/EV)/((1/EV)+lambda(1))^2;  %Exponential V*'(lambda(1)))
    LV3=(2/EV)/((1/EV)+lambda(1))^3;  %Exponential V*''(lambda(1)))
    WS1=((1-r1)*LS1)*(1-LV1)/((lambda(1)-l1*(1-LS1))*EV)
    WS2=((1-r1)/EV)*((lambda(1)-l1*(1-LS1))*(LS2*(1-LV1)-LS1*LV2)-(LS1*(1-LV1)*(1+l1*LS2)))/(lambda(1)-l1*(1-LS1))^2;
    WS3=0;
    for k1=2:M
        LS11=(1/ES)/((1/ES)+lambda(k1));    %Exponential S*(lambda(1)))
        LS21=-(1/ES)/((1/ES)+lambda(k1))^2; %Exponential S*'(lambda(1))) 
        LS31=(2/ES)/((1/ES)+lambda(k1))^3;  %Exponential S*''(lambda(1))) 
        LV11=(1/EV)/((1/EV)+lambda(k1));    %%Exponential V*(lambda(1)))
        LV21=-(1/EV)/((1/EV)+lambda(k1))^2;  %Exponential V*'(lambda(1)))
        LV31=(2/EV)/((1/EV)+lambda(k1))^3;  %Exponential V*''(lambda(1)))
        WS11=((1-r1)*LS11)*(1-LV11)/((lambda(k1)-l1*(1-LS11))*EV)
        WS21=((1-r1)/EV)*((lambda(k1)-l1*(1-LS11))*(LS21*(1-LV11)-LS11*LV21)-(LS11*(1-LV11)*(1+l1*LS21)))/(lambda(k1)-l1*(1-LS11))^2;
        WS3=WS3+lambda(k1)*(2/lambda(1)+WS21-2.0*WS11/lambda(1)) 
    end;
    AoI(m1)=EW+2.0*(ES+EV)+(1/lambda(1))*2.0*WS1-WS2-(1/lambda(1))+(ES+EV)*WS3
    m1=m1+1;

    end;
    x=(0.1:0.02:0.6)
%     plot(x,AoI,'-*');
     p(M)=plot(x,AoI,'*-');
    hold on;
end;
legend([p(2),p(3),p(4),p(5)],'M=2(Exp)','M=3(Exp)','M=4(Exp)','M=5(Exp)')


%% ======================= Part 7: Plotting with respect to EV=======================
fprintf('Plotting Data ...\n')
% number of sources
for M=2:5
   lambda=[];
    rho=[];
    AoI=[];
    AoI1=[];
% %Exponential
%     ES=0.2% expected service time
%     ES2=2*(ES)^2  %2nd order moment of    service time
% EV=1/2 %expected vacation time
% EV2=1/2 % second order moment of vacation
%Erlang-2
ES=0.2% expected service time
ES2=1*(ES)^2  %2nd order moment of service time

    lambda(1)=0.1;
    m1=1;EV=0.05;EV2=2*EV^2;
    for k=0.1:0.01:0.6
        %Exponential0
        EV=k;
        EV2=1*ES^2;
         %Erlang-2
         % EV=k% expected service time
        % EV2=1*(EV)^2  %2nd order moment of service time
        %HyperExponential with p
        % p=0.4;
         % EV=p*k+(1-p)*k^2% expected service time
        % EV2=p*2*k^2+(1-p)*2*k^4  %2nd order moment of service time
        
    for i=2:M
        lambda(i)=0.1;
    end
     l1=sum(lambda);
    for i=1:M
        rho(i)=lambda(i)*ES;
    end
    rho(1)
    r1=sum(rho);

    EW=(l1*ES2/2*(1-r1))+(EV2)/(2*EV); % Expected waiting Time
    LS1=(1/ES)/((1/ES)+lambda(1));    %Exponential S*(lambda(1)))
    LS2=-(1/ES)/((1/ES)+lambda(1))^2; %Exponential S*'(lambda(1))) 
    LS3=(2/ES)/((1/ES)+lambda(1))^3;  %Exponential S*''(lambda(1))) 
    LV1=(1/EV)/((1/EV)+lambda(1));    %%Exponential V*(lambda(1)))
    LV2=-(1/EV)/((1/EV)+lambda(1))^2;  %Exponential V*'(lambda(1)))
    LV3=(2/EV)/((1/EV)+lambda(1))^3;  %Exponential V*''(lambda(1)))
    WS1=((1-r1)*LS1)*(1-LV1)/((lambda(1)-l1*(1-LS1))*EV)
    WS2=((1-r1)/EV)*((lambda(1)-l1*(1-LS1))*(LS2*(1-LV1)-LS1*LV2)-(LS1*(1-LV1)*(1+l1*LS2)))/(lambda(1)-l1*(1-LS1))^2;
    WS3=0;
    for k1=2:M
        LS11=(1/ES)/((1/ES)+lambda(k1));    %Exponential S*(lambda(1)))
        LS21=-(1/ES)/((1/ES)+lambda(k1))^2; %Exponential S*'(lambda(1))) 
        LS31=(2/ES)/((1/ES)+lambda(k1))^3;  %Exponential S*''(lambda(1))) 
        LV11=(1/EV)/((1/EV)+lambda(k1));    %%Exponential V*(lambda(1)))
        LV21=-(1/EV)/((1/EV)+lambda(k1))^2;  %Exponential V*'(lambda(1)))
        LV31=(2/EV)/((1/EV)+lambda(k1))^3;  %Exponential V*''(lambda(1)))
        WS11=((1-r1)*LS11)*(1-LV11)/((lambda(k1)-l1*(1-LS11))*EV)
        WS21=((1-r1)/EV)*((lambda(k1)-l1*(1-LS11))*(LS21*(1-LV11)-LS11*LV21)-(LS11*(1-LV11)*(1+l1*LS21)))/(lambda(k1)-l1*(1-LS11))^2;
        WS3=WS3+lambda(k1)*(2/lambda(1)+WS21-2.0*WS11/lambda(1)) 
    end;
    AoI(m1)=EW+2.0*(ES+EV)+(1/lambda(1))*2.0*WS1-WS2-(1/lambda(1))+(ES+EV)*WS3
    m1=m1+1;

    end;
    x=(0.1:0.01:0.6)
%     plot(x,AoI,'-*');
    p(M)=plot(x,AoI,'+-');
    hold on;
end;
legend([p(2),p(3),p(4),p(5)],'M=2(Exp)','M=3(Exp','M=4(Exp)','M=5(Exp)')

