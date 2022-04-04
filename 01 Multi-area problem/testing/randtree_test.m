%%
noNodes = 25;
maxBranch = 4;
Adj = randTreeGraph(noNodes,maxBranch);
G = graph(Adj);
p = plot(G,'layout','force')