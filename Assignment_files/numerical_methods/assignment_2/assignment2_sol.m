%% Assignment #2 Solution
clear
clc
close 'all'



%% Problem #2
Fun = @(x) x.^4 - 200;
FunDer = @(x) 4*(x.^3);
Xest = 8;

Xs = NewtonSol(Fun,FunDer,Xest);
fprintf('Root from NewtonSol is: %2.15f \n',Xs);

%% Problem #3
fprintf('\n\n\n Problem #3 \n\n');
Fun = @(x) x - 2*exp(-x);
Xest = 2;
Xs = SteffensenRoot(Fun,Xest);
fprintf('Root from SteffensenRoot is: %2.15f \n',Xs);

%% Problem #4
fprintf('\n\n\n Problem #4 \n\n');
L = 4; % m, length of beam
w = 20e3; % N/m, distributed load
E = 70e9; % N/m^2, elastic modulus
I = 52.9e-6; % m^4, moment of inertia of the beam
k = w/(360*L*E*I);
Fun = @(x) k*(7*(L^4)*x - 10*(L^2)*(x.^3)+ 3*(x.^5));

FunDer = @(x) k*(7*(L.^4)-30*(L.^2)*(x.^2) + ...
    15*(x.^4));

FunDDer = @(x) k*(-60*(L.^2)*x + 60*(x.^3));

imax = 100;
Err = 1e-9;
Xa = 1.5;
Xb = 2.5;


Xs = SecantRoot(FunDer,Xa,Xb,Err,imax);
fprintf('Root from SecantRoot is: %2.15f \n',Xs);
fprintf('Position of maximum deflection is x = %2.15f \n',...
    Xs);
fprintf('Maximum deflection y = %g \n',Fun(Xs));

% Use NewtonSol
Xs = NewtonSol(FunDer,FunDDer,Xa);
fprintf('Root from NewtonSol is: %2.15f \n',Xs);
fprintf('Position of maximum deflection is x = %2.15f \n',...
    Xs);
fprintf('Maximum deflection y = %g \n',Fun(Xs));

% Use Steffensen Root
Xs = SteffensenRoot(FunDer,Xa);
fprintf('Root from SteffensenRoot is: %2.15f \n',Xs);
fprintf('Position of maximum deflection is x = %2.15f \n',...
    Xs);
fprintf('Maximum deflection y = %g \n',Fun(Xs));

%% Local Functions
function Xs = NewtonSol(Fun,FunDer,Xest)
Err = 1e-9;
imax = 100;

for k = 1:imax
   Xk = Xest - Fun(Xest)./FunDer(Xest);
   
   if abs((Xk - Xest)./Xest) < Err
      % relative error tolarance met
      Xs = Xk;
      return;
   end
   % if error tolerance not met, update Xest
   Xest = Xk;
   
   if k == imax
      error('Error! NewtonSol failed to converge in %d iterations.\n',...
          imax);
   end
    
end
end

function Xs = SteffensenRoot(Fun,Xest)
Err = 1e-9;
imax = 100;
for k = 1:imax
   Xk = Xest - (Fun(Xest).^2)./(Fun(Xest + Fun(Xest))-Fun(Xest));
     
   if abs((Xk - Xest)./Xest) < Err
      % relative error tolarance met
      Xs = Xk;
      return;
   end
   % if error tolerance not met, update Xest
   Xest = Xk;
   
   if k == imax
      error('Error! SteffenRoot failed to converge in %d iterations.\n',...
          imax);
   end
    
end

end

function Xs = SecantRoot(Fun,Xa,Xb,Err,imax)

for i = 1:imax
   FXb = Fun(Xb);
   Xi = Xb - FXb*(Xa - Xb)/(Fun(Xa) - FXb);
   
   if abs((Xi - Xb)/Xb) < Err
       Xs = Xi;
       break;
   end
   Xa = Xb;
   Xb = Xi;
     
        
end

end