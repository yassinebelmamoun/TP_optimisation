%% 2.2 Communication entre espions (optimisation combinatoire)

Donnees_probas=importdata('ProbaInterception.txt')
Donnees_probas=(log(Donnees_probas));
Donnees_probas(isnan(Donnees_probas))=0;
G=graph(Donnees_probas);
plot(G);

[Tree, pred]=minspantree(G);

p=plot(Tree);

display(['La probabilité minimale est donc : ', num2str(1-prod(1-exp((Tree.Edges.Weight))))])




