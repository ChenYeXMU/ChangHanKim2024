function F=benchmark_ReturnFn(h,aprime,a,z,r,w,sigma,B,gamma,tau,Tr)

F = -Inf;
y = w*z*h+r*a
c = (1-tau)*y+a+Tr-aprime; % Budget Constraint

if c>0
    if sigma == 1
        F = log(c)-B*(h^(1+gamma))/(1+gamma);
    else
        F = (c^(1-sigma))/(1-sigma)-B*(h^(1+gamma))/(1+gamma);
    end
end

end
