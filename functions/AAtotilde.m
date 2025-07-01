function tilde = AAtotilde(AA,isResonant)

tilde = zeros(6,1);

if isResonant
    AA(3) = AA(3)-AA(2);
    AA(5) = AA(5)+AA(6);
    ind = abs(AA) < 1e-15; %avoid issue of getting extremely small negative I3 when I2 and I3 are almost equal
    AA(ind) = 0;            %Come back for a better fix later
end

tilde(1) = sqrt(AA(1))*exp(real(AA(4)));

tilde(2) = sqrt(2*AA(2))*cos(AA(5));
tilde(3) = sqrt(2*AA(3))*cos(AA(6));
tilde(4) = sqrt(AA(1))*exp(-real(AA(4)));

tilde(5) = -sqrt(2*AA(2))*sin(AA(5));
tilde(6) = -sqrt(2*AA(3))*sin(AA(6));
temp = imag(AA(4));
switch temp
    case 0.5*pi
        tilde(1) = -tilde(1);
    case pi
        tilde(1) = -tilde(1);
        tilde(4) = -tilde(4);
    case 1.5*pi
        tilde(4) = -tilde(4);
end