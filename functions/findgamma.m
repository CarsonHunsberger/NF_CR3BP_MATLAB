function gamma = findgamma(Lpt,mu)
if Lpt==1
    p = [1 -(3-mu) (3-2*mu) -mu 2*mu -mu];
    g = roots(p);
    gamma = g(g==real(g));
end

if Lpt==2
    p = [1 (3-mu) (3-2*mu) -mu -2*mu -mu];
    g = roots(p);
    gamma = g(g==real(g));
end

if Lpt==3
    p = [1 2+mu 1+2*mu -(1-mu) -2*(1-mu) -(1-mu)];
    g = roots(p);
    gamma = g(g==real(g));
end