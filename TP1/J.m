function critere=J(H) %Fonction définissant le critère J à minimiser
K=0;
nu_j=[0:0.1/14:0.1,0.15:0.35/14:0.5];
for i=1:30
    Hj=H0(nu_j(i));
    Hnuj=H_nu(H,nu_j(i));
    if abs(Hj-Hnuj) >= K 
        K=abs(Hj-Hnuj);
    end
end    
critere=K;

function Hi=H_nu(H,nu) % Fonction qui définie la fonction H(nu) comme la réponse fréquentielle
Vecteur=[1:1:30];
Hi=cos(2*pi*Vecteur*nu)*H;
end

function Himp=H0(nu) % Fonction qui définit H0, la réponse fréquentielle idéale
if ((0 <= nu) & (nu<=0.1)) 
    Himp=1;
end    
if ((0.15<=nu) & (nu<=0.5)) 
    Himp=0;
end
end

end