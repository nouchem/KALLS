function KNN=Learn( S)

l=size( S.coord_active_set,1);
KNN = fitcknn(S.coord_active_set(1:l,:),S.label_active_set(1:l),'NumNeighbors',1,'distance','euclidean');%,'Standardize',1);

end