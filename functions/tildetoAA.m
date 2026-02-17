function AA = tildetoAA(tilde,isResonant)

AA = zeros(6,1);

AA(1) = abs(tilde(1)*tilde(4));
AA(2) = 0.5*(tilde(2)^2+tilde(5)^2);
AA(3) = 0.5*(tilde(3)^2+tilde(6)^2);
AA(4) = 0.5*log(abs(tilde(1)/tilde(4)));

if ~tilde(1) || ~tilde(4)
    AA(4) = 0;
else
    if tilde(1) < 0
        if tilde(4) < 0
            AA(4) = AA(4) + 1i*pi;
        else
            AA(4) = AA(4) + 1i*0.5*pi;
        end
    else
        if tilde(4) < 0
            AA(4) = AA(4) + 1i*1.5*pi;
        end
    end
end

% AA(5) = -atan2(tilde(5),tilde(2)); %OLD
% AA(6) = -atan2(tilde(6),tilde(3)); %OLD
AA(5) = atan2(-tilde(2),-tilde(5)); %NEW
AA(6) = atan2(-tilde(3),-tilde(6)); %NEW

if isResonant
    AA(3) = AA(2)+AA(3);
    AA(5) = AA(5)-AA(6);
end


end