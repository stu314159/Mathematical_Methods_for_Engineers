function y =  myThreePointFifthOrder(f,xMin,xMax,N)
xI = NaN(3,1);
xI(1) = 0.1127016653792583;
xI(2) = 0.5;
xI(3) = 0.887298334620755;

wI = NaN(3,1);
wI(1) = 0.277777777777777786;
wI(2) = 0.444444444444444493;
wI(3) = wI(1);

xSplit = linspace(xMin,xMax,N+1);
Jac = xSplit(2) - xSplit(1);

y = 0;
for d = 1:N
    a = xSplit(d); b = xSplit(d+1);
    xT = a + (b-a)*xI;
    y = y + f(xT)'*wI*Jac;
end