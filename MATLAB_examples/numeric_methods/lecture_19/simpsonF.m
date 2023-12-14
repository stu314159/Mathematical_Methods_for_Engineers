function w = simpsonF(xW)

A = [1 1 1;
    0 xW(1) 1;
    0 xW(1)^2 1;
    0 xW(1)^3 1];

w = A*[xW(2);xW(3);xW(4)] - [1; 1/2; 1/3; 1/4];

end