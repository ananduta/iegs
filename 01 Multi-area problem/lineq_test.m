% test

A = [1 -1 0;
     0 1 -1];

b = rand(size(A,1),1);

x0 = A\b;

x1 = pinv(A)*b;

xd = x0-x1;