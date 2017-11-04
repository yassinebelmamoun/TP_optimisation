function GradResultsAdaptatif = gradient_rho_adaptatif(han_f, han_df, U0, pas, tol)
% Fonction permettant de minimiser la fonction f(U) par rapport au vecteur U 
% M�thode : gradient � pas fixe
% INPUTS :
% - han_f   : handle vers la fonction � minimiser
% - han_df  : handle vers le gradient de la fonction � minimiser
% - U0      : vecteur initial 
% - rho0     : Valeur initiale du pas 
% - tol     : tol�rance pour d�finir le crit�re d'arr�t
% OUTPUT : 
% - GradResultsAdaptatif : structure d�crivant la solution 

itermax=10000;  % nombre maximal d'it�rations 

xn=U0; f=han_f(xn); % point initial de l'algorithme
it=0;          % compteur pour les it�rations
rhon=pas; 
converged = false;

while ~converged && it < itermax  
    dfx=han_df(xn);       % valeur courante de la fonction � minimiser
    xnp1=xn-rhon*dfx;      % nouveau point courant (x_{n+1})
    fnp1=han_f(xnp1);
    if fnp1<han_f(xn)    
        rhonp1=2*rhon;
        if abs(fnp1-f)<tol
            converged = true;
        end
        xn=xnp1; f=fnp1;     it=it+1;
    else
        rhonp1=rhon/2;
    end    
    rhon=rhonp1;     % xnp1 : nouveau point courant
end

GradResultsAdaptatif.initial_x=U0;         % vecteur initial
GradResultsAdaptatif.minimum=xnp1;         % vecteur apr�s optimisation
GradResultsAdaptatif.f_minimum=fnp1;       % valeur optimale de la fonction
GradResultsAdaptatif.cout=it;        % nombre d'appels de la fonction
GradResultsAdaptatif.converged=converged;  % true si l'algorithme a converg�