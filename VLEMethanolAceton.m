% Dieser Skript berechnet das VLE von dem binären Gemisch Aceton-Methanol
% bei 55°C
clc
clear

% Array für experimentelle Daten
% Innere Struktur ist P|x1|y1

X_exp = readmatrix("VLE Aceton-Methanol.xlsm","Sheet","Experimentelle Daten","Range","A2:C23");
l = length(X_exp(:,1));


% Now the experimental activity coefficients are calculated:
% Inner structure x1 | gamma1 | gamma2 | x2
% Formula: gamma(i) = (y(i)*P)/(x*PLV(i))

PLV1 = X_exp(l,1);
PLV2 = X_exp(1,1);
display(PLV1)
gammaExp = zeros(l-2,4);
for i=2:(l-1)
    gammaExp(i-1,1) = X_exp(i,2);
    gammaExp(i-1,2) = (X_exp(i,3)*X_exp(i,1))/(X_exp(i,2)*PLV1);
    gammaExp(i-1,3) = ((1-X_exp(i,3))*X_exp(i,1))/((1-X_exp(i,2))*PLV2);
    gammaExp(i-1,4) = 1-X_exp(i,2);
end

% Calculating the experimental Excess Energy
% Equation: gE/RT = sum(xi*ln(gamma(i))

gE_Exp = zeros(l,2);    %Internal Structure: x1 | gE/RT
deltaG = zeros(l,2);
for i=1:l
    if i==1
        gE_Exp(i,2) = 0;
        gE_Exp(i,1) = 0;
        deltaG(i,1) = 0;
        deltaG(i,2) = 0;
    elseif i==l
        gE_Exp(i,2) = 0;
        gE_Exp(i,1) = 1;
        deltaG(i,1) = 0;
        deltaG(i,2) = 0;
    else 
        gE_Exp(i,2) = gammaExp(i-1,1).*log(gammaExp(i-1,2)) + gammaExp(i-1,4).*log(gammaExp(i-1,3));
        gE_Exp(i,1) = gammaExp(i-1,1);
        deltaG(i,1) = gammaExp(i-1,1);
        deltaG(i,2) = gammaExp(i-1,1).*log(gammaExp(i-1,1)) + gammaExp(i-1,4).*log(gammaExp(i-1,4)) + gE_Exp(i,2);
    end
     
end

%¿ Least square regression:
% X contains the x1 values
% y contains the corresponding experimental gE/RT values
% ModelFunction Implements the model function

x = gE_Exp(:,1);
y = gE_Exp(:,2);
ft = fittype('ModelFunction(x, A12, A21)');
f = fit(x,y,ft,"StartPoint", [0.5, 0.5]);
A = coeffvalues(f);     % Structure is A12 | A21
gamma_calc = ones(length(x),4);    % Structure is x1 | gamma1 | gamma2 | x2

h = 0.00001;
x = linspace(0,1,100);
% Now calculate the "Siede- und Taulinie"
% First The Siedelinie:
SiedeLinie = zeros(length(x),2);       % Structure is (x,P) point pairs
SiedeLinie(:,1) = x;
TauLinie = zeros(length(x),2);       % Structure is (x,P) point pairs
TauLinie(:,1) = x;

for i=1:(length(x))

    % First calculate the activity coefficients
    DerivativeX1 = (ModelFunction(x(i)+h,A(1,1),A(1,2)) - ModelFunction(x(i),A(1,1),A(1,2)))/h;
    DerivativeX2 = (-1).*DerivativeX1;
    gamma_calc(i,1) = x(i);
    gamma_calc(i,2) = exp(DerivativeX1.*(1-x(i)) + ModelFunction(x(i),A(1,1),A(1,2)));
    gamma_calc(i,3) = exp(DerivativeX2.*x(i) + ModelFunction(x(i),A(1,1),A(1,2)));
    gamma_calc(i,4) = 1 - x(i);

    % Now calculate the pressures of the Siedelinie:
    syms p
    eqn1 = 1 == (gamma_calc(i,1).*PLV1.*gamma_calc(i,2))/p + (gamma_calc(i,3).*PLV2.*gamma_calc(i,4))/p;
    SiedeLinie(i,2) = solve(eqn1,p);

    % Now calculate the corresponding pressures of the Taulinie:
    TauLinie(i,1) = (SiedeLinie(i,1).*PLV1.*gamma_calc(i,2))/(SiedeLinie(i,2));
    TauLinie(i,2) = SiedeLinie(i,2);
    
end

scatter(X_exp(:,2), X_exp(:,1), "red"); hold on
title("Aceton-Methanol VLE at T=55°C")
xlabel("x1 [-]")
ylabel("Pressure [kPa]")
scatter(X_exp(:,3), X_exp(:,1), "red"); hold on
p1 = plot(SiedeLinie(:,1), SiedeLinie(:,2)); hold on
p1.Color="black";
p2 = plot(TauLinie(:,1), TauLinie(:,2));
p2.Color="black";






    


             