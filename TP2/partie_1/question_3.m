n = 15;


load('DonneesRangerObjets.mat');
n = length(PositionCasiers);
d = [];
for i=1:n
    for j=1:n
        d(i,j)=abs(PositionCasiers(j)-PositionObjets(i));
    end
end
d_vector = d(:);


% Define A and b
E = kron(eye(n),ones(1, n));
F = kron(ones(1,n), eye(n));

G1_bis = [1, zeros(1,n-1)];
G1 = kron(eye(n-1,n), G1_bis);
G2_bis = [0, -1, zeros(1,n-2)];
G2 = [zeros(n-1,n),kron(eye(n-1, n-1), G2_bis)];
G3 = [0, 1, zeros(1, n^2-2)];
G = [G1 + G2; G3];

A = [];
b = [];
Aeq = [E;F;G];
beq = [ones(2*n,1);zeros(n,1)];
lb = zeros(225,1);
ub = ones(225,1);

% Mathematical Problem
[x, fval, exitflag, output] = linprog(d_vector,A,b,Aeq, beq, lb, ub);




