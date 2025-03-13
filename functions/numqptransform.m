function  x = numqptransform(x,GenFuncEOMs,GenFuncEOMseval,flag,indarr,opts)
for n=1:length(indarr)
    m = indarr(n);
    [~,rv] = ode113(@(t,x)GenFuncEOMseval(t,x,GenFuncEOMs{m}),[0 flag],x,opts);
    x= rv(end,:)';
end
end