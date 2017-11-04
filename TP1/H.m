function Himp=H(nu)
if ((0 <= nu) & (nu<=0.1)) 
    Himp=1;
end    
if ((0.15<=nu) & (nu<=0.5)) 
    Himp=0;
end
