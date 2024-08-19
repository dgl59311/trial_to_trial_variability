function gmd=GlobalMapDiss(x,y)
    n=length(x);
    u=mean(x);
    v=mean(y);
    den_1=sqrt(sum(((x-u).^2)/n));
    den_2=sqrt(sum(((y-v).^2)/n));
    gmd=sqrt(sum((((x-u)/den_1)-((y-v)/den_2)).^2)/n);
end