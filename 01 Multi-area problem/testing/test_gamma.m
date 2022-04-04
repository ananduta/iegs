gamma(1) = 0;
r_max = 100;

gammaUpper = inf;
gammaLower = gamma(1);
d =0.6;
for r = 1:r_max
    er(r) = error_func(gamma(r));
    
    [gamma(r+1),gammaLower,gammaUpper]=gammaRules(er(r),gamma(r),gammaLower,gammaUpper);
   
   
   
   
end

figure
subplot(2,1,1)
plot(gamma)
subplot(2,1,2)
plot(er)