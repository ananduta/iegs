function y = error_func(x)


c = 400;

b = -sqrt(c*4);

x0 = -b/2;

%x=[0:0.1:x0+10];

%for t = 1:length(x)
    if x < x0
        y = x^2 + b*x+c;
    else
        y = 0;
    end
%end
end
