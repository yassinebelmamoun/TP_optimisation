close all, clear all
[Aineq,B,S] = definition_constantes;



%% gradient avec rho constant

%rho = 0.1
d=5; U0 = zeros(d,1);
tolerance=1e-6; pas1=0.01; 
han_f1  = @(u) f1(u,B,S);
han_df1 = @(u) df1(u,B,S);
tic
GradResults=gradient_rho_constant(han_f1,han_df1,U0,pas1,tolerance)
toc
tic
options=optimoptions('fminunc','algorithm','quasi-newton','TolFun',tolerance)
[X,FVAL,EXITFLAG,OUTPUT]=fminunc(han_f1,U0,options)
toc
disp(['gradient_rho_constant pour rho=0.01: minimum=',num2str(GradResults.f_minimum)])

%rho = 0.05
d=5; U0 = zeros(d,1);
tolerance=1e-6; pas2=0.05; 
han_f1  = @(u) f1(u,B,S);
han_df1 = @(u) df1(u,B,S);
tic
GradResults=gradient_rho_constant(han_f1,han_df1,U0,pas2,tolerance)
toc
tic
options=optimoptions('fminunc','algorithm','quasi-newton','TolFun',tolerance)
[X,FVAL,EXITFLAG,OUTPUT]=fminunc(han_f1,U0,options)
toc
disp(['gradient_rho_constant pour rho=0.05: minimum=',num2str(GradResults.f_minimum)])

%Rho=0.1
d=5; U0 = zeros(d,1);
tolerance=1e-6; pas3=0.1; 
han_f1  = @(u) f1(u,B,S);
han_df1 = @(u) df1(u,B,S);
tic
GradResults=gradient_rho_constant(han_f1,han_df1,U0,pas3,tolerance)
toc
tic
options=optimoptions('fminunc','algorithm','quasi-newton','TolFun',tolerance)
[X,FVAL,EXITFLAG,OUTPUT]=fminunc(han_f1,U0,options)
toc
disp(['gradient_rho_constant pour rho=0.1: minimum=',num2str(GradResults.f_minimum)])



%% gradient avec rho adaptatif

%rhoinitial=0.01
tic
GradResultsAdaptatif=gradient_rho_adaptatif(han_f1,han_df1,U0,pas1,tolerance)
toc
disp(['gradient_rho_constant pour rho_initial=0.01: minimum=',num2str(GradResultsAdaptatif.f_minimum)])


%rhoinitial=0.05
tic
GradResultsAdaptatif=gradient_rho_adaptatif(han_f1,han_df1,U0,pas2,tolerance)
toc
disp(['gradient_rho_constant pour rho_initial=0.05: minimum=',num2str(GradResultsAdaptatif.f_minimum)])


%rhoinitial=0.1
tic
GradResultsAdaptatif=gradient_rho_adaptatif(han_f1,han_df1,U0,pas3,tolerance)
toc
disp(['gradient_rho_constant pour rho_initial=0.1: minimum=',num2str(GradResultsAdaptatif.f_minimum)])



%% 1.1.2 Utilisation des routines d'optimisation : Méthode de Quasi Newton

%Sans renseigner le gradient de f
tolerance=1e-6;
tic
options1=optimoptions('fminunc','algorithm','quasi-newton','TolFun',tolerance)
[X,FVAL,EXITFLAG,OUTPUT]= fminunc(han_f1,U0,options1)
toc

%En renseignant le gradient de f
han_f1_df1=@(u) f_df_quasi_newton(u,B,S);

tic
options2=optimoptions('fminunc','SpecifyObjectiveGradient',true)
[X,FVAL,EXITFLAG,OUTPUT]= fminunc(han_f1_df1,U0,options2)
toc


%% Partie 1.2 : Optimisation sous contraintes
close all, clear all
[Aineq,B,S] = definition_constantes;

%1.2.1 Utilisation de l'algorithme SQP pour f1 
d=5; U0 = zeros(d,1);
A=[]; b=[]; Aeq=[]; beq=[]; lb=[0,0,0,0,0]; ub=[1,1,1,1,1];
han_f1  = @(u) f1(u,B,S);
tic
options=optimoptions('fmincon','Display','iter','Algorithm','sqp')
[X,FVAL,EXITFLAG,OUTPUT]=fmincon(han_f1,U0,A,b,Aeq,beq,lb,ub)
toc

%% Utilisation de l'algorithme SQP pour f2 

d=5; U0 = zeros(d,1);
A=[]; b=[]; Aeq=[]; beq=[]; lb=[0,0,0,0,0]; ub=[1,1,1,1,1];
han_f2  = @(u) f2(u,B,S);
tic
options=optimoptions('fmincon','Display','iter','Algorithm','sqp')
[X,FVAL,EXITFLAG,OUTPUT]=fmincon(han_f2,U0,A,b,Aeq,beq,lb,ub)
toc

%% 1.2.2  Optimisation sous contrainte et penalisation
epsilon=1e-15;
tolerance=1e-20;
han_f1penal  = @(u) f1penal(u,B,S,epsilon);
han_f2penal  = @(u) f2penal(u,B,S,epsilon);

tic
options1 = optimoptions('fminunc','algorithm','quasi-newton','TolFun',tolerance)
[X,FVAL,EXITFLAG,OUTPUT] = fminunc(han_f1penal,U0,options1)
toc

tic
options2 = optimoptions('fminunc','algorithm','quasi-newton','TolFun',tolerance)
[X,FVAL,EXITFLAG,OUTPUT] = fminunc(han_f2penal,U0,options2)
toc


%% 1.2.3 Methodes duales pour l'optimisation sous contraintes Uzawa et lagrangien
lambdan=zeros(10,1);
converged=false;
[Aineq,B,S] = definition_constantes;
U0=zeros(5,1);
pas=0.01;
tolerance=1e-6;
    while ~converged
        han_Lagrangien =@(u) Lagrangien(u,B,S,lambdan);
        han_dLagrangien_U =@(u) dLagrangien_U(u,B,S,lambdan);
        GradResultsAdaptatif=gradient_rho_adaptatif(han_Lagrangien,han_dLagrangien_U,U0,pas,tolerance);
        xn=GradResultsAdaptatif.minimum;
        lambdanp1=max(0,lambdan + dLagrangien_lambda(xn,lambdan));
        if abs(lambdanp1-lambdan)<1e-3
            converged=true;
        else
            lambdan=lambdanp1;
            xnp1=xn;
        end
    end    
disp(han_f1(xn))


%% 1.3 Optimisation non convexe - Recit simulé

%U0 = [0;0;0;0;0]
close all, clear all
[Aineq,B,S] = definition_constantes;
tolerance=1e-6;
U01 = [0;0;0;0;0];
lb=[]; ub=[];
han_f1  = @(u) f1(u,B,S);
han_f4  = @(u) f4(u,B,S);
tic
options = optimoptions('fminunc','algorithm','quasi-newton','TolFun',tolerance)
[X,FVAL,EXITFLAG,OUTPUT] = fminunc(han_f4,U01,options)
toc
tic
options = saoptimset ('InitialTemperature',50, 'ReannealInterval',10)
[X,FVAL,EXITFLAG,OUTPUT] = simulannealbnd(han_f4,U01,lb,ub,options)
toc



%U0 = [0.5;0.5;0.5;0.5;0.5]
close all, clear all
[Aineq,B,S] = definition_constantes;
tolerance=1e-6;
U02 = [0.5;0.5;0.5;0.5;0.5];
lb=[]; ub=[];
han_f1  = @(u) f1(u,B,S);
han_f4  = @(u) f4(u,B,S);
tic
options = optimoptions('fminunc','algorithm','quasi-newton','TolFun',tolerance)
[X,FVAL,EXITFLAG,OUTPUT] = fminunc(han_f4,U02,options)
toc
tic
options = saoptimset ('InitialTemperature',50, 'ReannealInterval',10)
[X,FVAL,EXITFLAG,OUTPUT] = simulannealbnd(han_f4,U02,lb,ub,options)
toc





%U0 = [20;5;5000;1000;50]
close all, clear all
[Aineq,B,S] = definition_constantes;
tolerance=1e-6;
U03 = [20;5;5000;1000;50];
lb=[]; ub=[];
han_f1  = @(u) f1(u,B,S);
han_f4  = @(u) f4(u,B,S);
tic
options = optimoptions('fminunc','algorithm','quasi-newton','TolFun',tolerance)
[X,FVAL,EXITFLAG,OUTPUT] = fminunc(han_f4,U03,options)
toc
tic
options = saoptimset ('InitialTemperature',50, 'ReannealInterval',10)
[X,FVAL,EXITFLAG,OUTPUT] = simulannealbnd(han_f4,U03,lb,ub,options)
toc

%% 1.4 Application Synthèse d’un filtre à réponse impulsionnelle finie
tolerance=1e-10;
nu_j=[0:0.1/14:0.1,0.15:0.35/14:0.5];
Hinitial=ones(30,1);
han_J = @(h) J(h);
tic
options=optimoptions('fminunc','algorithm','quasi-newton','TolFun',tolerance)
[X,FVAL,EXITFLAG,OUTPUT]= fminunc(han_J,Hinitial,options)
toc
 




