function F=CHK2024v2_ReturnFn(h,aprime,a,z,r,w,sigma,B,gamma,tau,Tr,Reform)

F = -Inf;
c = 0;
y = w*z*h+r*a;

if Reform == 0
    c = y -aprime;
elseif Reform == 1
    c = (1-tau)*y+a+Tr-aprime; % Budget Constraint
elseif Reform ==2
    c = y-tau*(y-Tr)+a-aprime;
end

if c>0
    if sigma == 1
        F = log(c)-B*(h^(1+gamma))/(1+gamma);
    else
        F = (c^(1-sigma))/(1-sigma)-B*(h^(1+gamma))/(1+gamma);
    end
end

end
