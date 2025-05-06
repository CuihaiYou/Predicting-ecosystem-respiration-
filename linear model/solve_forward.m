%this function solves the forward problem
function y = solve_forward(tau,GPP,B,a,c,x0,Nt)
A = [-1     0       0         0          0          0
    0      -1      0         0          0          0
    1      1       -1        0          0          0
    0      0       a(1)     -1         a(5)       a(7)
    0      0       a(2)     a(3)       -1         0
    0      0       0        a(4)       a(6)       -1 ];

C=[c(1)        0           0           0           0            0
    0            c(2)       0           0           0            0
    0            0           c(3)       0           0            0
    0            0           0          c(4)        0            0
    0            0           0           0           c(5)        0
    0            0           0           0           0           c(6)]; %the constant matrix A

x_last=x0';
co2=zeros(1,5);% 
x  = zeros(12,365*9); %
n_record = 0;
 for i_yr = 2009:2017
    i_d=(i_yr-2008-1)*365+1;
    i_d2=(i_yr-2008)*365;
    for i=i_d:i_d2
        %     display(i);
        tau_used = diag(squeeze(tau(:,i)));
        x_present = [eye(6)+A*C*tau_used]*x_last + B'.*GPP(i);
        co2(1)= tau(4,i)*c(4)*(1-a(3)-a(4))*x_last(4);
        co2(2)= tau(5,i)*c(5)*(1-a(5)-a(6))*x_last(5);
        co2(3)= tau(6,i)*c(6)*(1-a(7))*x_last(6);
        co2(4)=tau(3,i)*c(3)*(1-a(1)-a(2))*x_last(3);
        co2(5)= 0.45*GPP(i);
        x_last = x_present;
       if i== 302||i== 667||i== 1032||i== 1397||i== 1762||i== 2127||i== 2492||i== 2857||i== 3222||i== 3587||i== 3952||i== 4317||i== 4682;
            x_last(3) = x_last(3) + x_last(1);
            x_last(1) = 0;
        end
        n_record=n_record+1;
        x(1:6,n_record) = x_present;
        x(7:11,n_record) = co2; %x(7:9,n_record) = co2
        x(12,n_record) = sum(co2); %x(10,n_record) = sum(co2)
    end
 end

y = x;
