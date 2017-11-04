function GradResults = gradient_rho_constant(han_f, han_df, U0, rho, tol)
% Fonction permettant de minimiser la fonction f(U) par rapport au vecteur U 
% M�thode : gradient � pas fixe
% INPUTS :
% - han_f   : handle vers la fonction � minimiser
% - han_df  : handle vers le gradient de la fonction � minimiser
% - U0      : vecteur initial 
% - rho     : param�tre g�rant l'amplitude des d�placement 
% - tol     : tol�rance pour d�finir le crit�re d'arr�t
% OUTPUT : 
% - GradResults : structure d�crivant la solution 

itermax=10000;  % nombre maximal d'it�rations 

xn=U0; f=han_f(xn); % point initial de l'algorithme
it=0;          % compteur pour les it�rations
converged = false;

while ~converged && it < itermax  
    it=it+1;
    dfx=han_df(xn);       % valeur courante de la fonction � minimiser
    xnp1=xn-rho*dfx;      % nouveau point courant (x_{n+1})
    fnp1=han_f(xnp1);
    if abs(fnp1-f)<tol
        converged = true;
    end
    xn=xnp1; f=fnp1;      % xnp1 : nouveau point courant
end

GradResults.initial_x=U0;         % vecteur initial
GradResults.minimum=xnp1;         % vecteur apr�s optimisation
GradResults.f_minimum=fnp1;       % valeur optimale de la fonction
GradResults.iterations=it;        % nombre d'it�rations
GradResults.converged=converged;  % true si l'algorithme a converg�