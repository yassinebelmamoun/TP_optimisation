n = 15;
load('DonneesEnginsChantier.mat', 'D');
n = length(D);

% Fonction à minimiser

d = [];
for i=1:n
    d(1,i)=200; d(2,i)=800; d(3,i)=1200;
end
d_vector = d(:);

% Contraintes de chantier: Il faut avoir le nombre d'engins nécessaire

A = kron(eye(n), [-1, 0, 0]);
b = (-1) * transpose(D);

ub = [];
lb = zeros(3*n,1);

% Contraintes sur la conservation du nombre d'engins de semaine en semaine

Aeq_1 = [zeros(1, 3*n); [kron(eye(n-1),[1,0,0]), zeros(n-1, 3)]];
Aeq_2 = kron(eye(n),[-1,1,-1]);
Aeq = Aeq_1 + Aeq_2;
beq = zeros(n, 1);

% Résolution avec la fonction intlinprog

[x, fval, exitflag, output] = intlinprog(d_vector, [1:3*n], A,b, Aeq, beq, lb, ub);

for i=1:length(x)
    if mod(i-1,3) == 0
        disp(x(i));
    end
end



