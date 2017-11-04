n = 15;

% get O and B
load('DonneesRangerObjets.mat');
n = 15;
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
%------------------ Question 4

H_1 = [0,0,1,zeros(1,n-3)];
H_2 = [0,0,0,1,zeros(1,n-4)];
A_1 = [H_1,kron(ones(1,n-1),H_2)];
for i=1:(n-1)
    B = [zeros(1,i*n), H_1, kron(ones(1, n-1-i), H_2)];
    A_1 = [A_1; B];
end
%------------------ Question 5
% i = 1
simple_seven = [zeros(1, 6), 1, zeros(1, n-7)];
simple_nine = [zeros(1,8), 1, zeros(1,n-9)];
simple_zero = [zeros(1,n)];
for i=1:n
    for j=1:n
        if abs(j-i)>1
            A_9_aux(i,j) = 1;
        else 
            A_9_aux(i,j) = 0;
        end
    end
end
A_2 = kron(eye(n), simple_seven) + kron(A_9_aux, simple_nine)
A = [A_1;A_2]
b = ones(2*n,1);
Aeq = [E;F;G];
beq = [ones(2*n,1);zeros(n,1)];
lb = zeros(225,1);
ub = ones(225,1);

% Mathematical Problem
[x, fval, exitflag, output] = intlinprog(d_vector, [1:n^2], A,b,Aeq, beq, lb, ub);





