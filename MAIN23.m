%%Partie 2.3 : Methode gloutonne 

N=1000;
ainf=0.02;
asup=1;
binf=0;
b=[];
r=(asup-ainf)*rand(N,1)+ainf;
bsup=r-0.01*ones(N,1);
    for i=1:N
b(i)=(bsup(i)-binf)*rand(1,1)+binf;
    end 
p_values=[];
d_values=[];

for j=1:N
    p_values(j)=p(r(j),b(j));
    d_values(j)=d(r(j),b(j));
end

plot(p_values,d_values,'bx');




    p_values_pareto=[];
    d_values_pareto=[];
for i=1:N
    domine=0;
    for j = 1:N
        if ((p_values(i) <= p_values(j)) )
            domine=domine+1;
        end
    if domine == N-1
        p_values_pareto=[p_values_pareto;p_values(i)];
        d_values_pareto=[d_values_pareto;d_values(i)];
    end
    end
end    

%%Partie 2.3 : Methode plus sophistiquée


