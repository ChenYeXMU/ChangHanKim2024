function F=CHK2024_ReturnFn(h,aprime,a,z,r,w,sigma,B,gamma,tau,Tr,Reform)

F = -Inf;
c = 0;

if Reform == 0
    c = w*z*h+(1+r)*a -aprime;
elseif Reform == 1
    c = (1-tau)*(w*z*h+r*a)+a+Tr-aprime; % Budget Constraint
end

if c>0
    if sigma == 1
        F = log(c)-B*(h^(1+gamma))/(1+gamma);
    else
        F = (c^(1-sigma))/(1-sigma)-B*(h^(1+gamma))/(1+gamma);
    end
end

end
