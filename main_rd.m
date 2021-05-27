
%This code is an implementation of the K-nn Active Learning Under Local Smoothness (KALLS). For theoritical details see the publication: "K-NN active learning under local smoothness assumption" by  B. N. Njike and X. Siebert

%Input:
%  epsilon: accuracy parameter
%  delta: Confidence parameter
%  c: margin noise parameters
%  L: smoothness parameters
%  n:label budget   
%Output:
%  activeErrorf
%  passiveErrorf
close all
clear all
clc

synth_data=1;
reeldata=1;


switch reeldata
    case 1
        load('HTRU_2.csv');
        data_X=HTRU_2(:,1:8);
        data_Y=HTRU_2(:,9);
    case 2
         load('SUSY.mat');
        data_X=SUSY(:,2:19);
        data_Y=SUSY(:,1);
end

alpha=1/6:1/6:1;
beta=[0,1,2,3,4,6,12,33];

for alphai=1:length(alpha)
   for betai=1:length(beta)
        
for j=1:10
perm=randperm(size(data_X,1));
data_X=data_X(perm,:);
data_Y=data_Y(perm);
Train_X=data_X(1:(size(data_X,1)*0.8),:);
Train_Y=data_Y(1:(size(data_X,1)*0.8));
Test_X=data_X((size(data_X,1)*0.8)+1:size(data_X,1),:);
Test_Y=data_Y((size(data_X,1)*0.8)+1:size(data_X,1));

valCoord_tree = createns(Train_X,'nsmethod','kdtree');      

ii=1;      
c=2^(beta(betai))*7;        
L=5*log(1/epsilon);        
n=round(.1*size(data_X,1));        
d=size(Train_X,2);        

passiveError=[];
activeError=[];


Q={};
I={};
S={};
T=1;
r=0;
i=1;
epsilon=0.09;        
delta=0.1;  
    
deltachap=maxi(epsilon,beta(betai),c);
jj=1;
   t=n;
   %first_selected_point
   r=r+1;
   deltas=delta/(32*(ii^2));
        fprintf('Points queried by Active Learning algorithm : %03d\n',r);
        
        while size(S,1)==0
        Index=ii;
        ks=minimum(deltachap,deltas);
        ks=min(ks,t);
        format long g
        [ind,N] = knnsearch(Train_tree,Train_X(Index,:),'K',round(ks+1));
        N(N==0)=[]; 
        ind(ind==Index)=[];
        
        [rcoord,rlabel,Y]=confidentlabel(ii,ks,deltas,Train_X(ind,:),Train_Y(ind));
        t1=size(rcoord,1);
        dist(r)=pdist2(rcoord,Train_X(Index,:),'euclidean','Largest',1);
        
        Q.coord_N{r}=rcoord;
        Q.label_N{r}=rlabel;
        Q.indices_N{r}=ind(1:t1);
        
        Bs(r)=sqrt((2/t1)*(log(1/deltas)+log(log(1/deltas))+log(log(exp(1)*t1))));
        LBs(r)=abs(((1/t1)*sumlabel(rlabel))-(1/2))-(Bs(r));
        
        I.coord_Train(r,:)=Train_X(ii,:);
        I.LB(r)=LBs(r);
        I.size_N(r)=t1;
        I.indices_Train(r)=Index;
        I.dist(r)=dist(r);
        
        if LBs(r)>=0.1*Bs(r)
        S.coord_active_set(i,:)= I.coord_Train(r,:); 
        S.label_active_set(i)=Y;
        S.indices_active_set(i)=Index;
        S.LB(i)=LBs(r);
        S.dist(i)=dist(r);
        S.Neighbors_active_set{i}=Q.coord_N{r};%rcoord;
        S.label_Neighbors_active_set{i}=Q.label_N{r};
        S.Neighbors_indices_active_set{i}=Q.indices_N{r};
        S.size_N(i)=I.size_N(r);      
        
        KNN1=LearnF(S);
        
        KNN=fitcknn(Train_X(1:sum(S.size_N(1:i)),:),Train_Y(1:sum(S.size_N(1:i))),'NumNeighbors',1,'distance','euclidean');%classifier 1NN passive
        [testPredp, prob] = predict(KNN,Test_X);
        [testPreda, prob] = predict(KNN1,Test_X);
                                                                             
        passiveError(i,1)=error_calc(testPredp,Test_Y);
        passiveError(i,2)=sum(S.size_N(1:i));
        activeError(i,1)=error_calc(testPreda,Test_Y);
        activeError(i,2)=sum(S.size_N(1:i));
        
        i=i+1;
        end
        
        t=t-t1;
        ii=ii+1;
        r=r+1;      
        end
        while t>0 & ii< size(Train_X,1)
             fprintf('Points queried by Active Learning algorithm : %03d\n',r);
             [ind1,N1] = knnsearch(valCoord_tree,Train_X(ii,:),'K',size(Train_X,1));
             N1(N1==0)=[]; 
             ind(ind==Index)=[];
   
             T=reliable(Train_X(ii,:),alpha(alphai),L,S,d, N1,size(Train_X,1));
       
             if T==1
                 r=r+1;
                 Index=ii;
                 deltas=delta/(32*(ii^2));
                 ks=minimum(deltachap,deltas);
                 format long 
                 [ind,N] = knnsearch(Train_tree,Train_X(Index,:),'K',round(ks+1));
                 N(N==0)=[]; 
                 ind(ind==Index)=[];
                 coord_v=[];
         
                 for i_1=1:size(S.coord_active_set,1)
                     if norm(S.coord_active_set(i_1,:)-Train_X(ii,:))<=S.dist(i_1)
                         coord_v=S.Neighbors_active_set{i_1};
                         break;
                     end
                 end
                 
                 if isempty(coord_v)==1
         
                     [rcoord,rlabel,Y]=confidentlabelF(ii,ks,deltas,Train_X(ind,:),Train_Y(ind));
                     t1=size(rcoord,1);
                 else
                     [rcoord,rlabel,Y,t1,s1]=confidentlabel1(ii,ks,deltas,Train_X(ind,:),Train_Y(ind),S,i_1,ind);
                     N_used(jj,1)=ii;
                     N_used(jj,2)=s1;
                     jj=jj+1;
         
                 end
                 
                 dist=pdist2(rcoord,Train_X(Index,:),'euclidean','Largest',1);

                 Q.coord_N{r}=rcoord;
                 Q.label_N{r}=rlabel;
                 Q.indices_N{r}=ind(1:size(rcoord,1));
        
                 Bs(r)=sqrt((2/size(rcoord,1))*(log(1/deltas)+log(log(1/deltas))+log(log(exp(1)*size(rcoord,1)))));
                 LBs(r)=abs(((1/size(rcoord,1))*sumlabel(rlabel))-(1/2))-(Bs(r));
      
                 I.coord_Train(r,:)=Train_X(ii,:);
                 I.LB(r)=LBs(r);
                 I.size_N(r)=t1;
                 I.indices_Train(r)=Index;
                 I.dist(r)=dist;
                 
                 if LBs(r)>=0.1*Bs(r)
                     S.coord_active_set(i,:)= I.coord_Train(r,:); 
                     S.label_active_set(i)=Y;
                     S.indices_active_set(i)=Index;
                     S.LB(i)=LBs(r)
                     S.dist(i)=dist;  
                     S.size_N(i)=I.size_N(r);
                     S.Neighbors_active_set{i}=rcoord;
                     S.label_Neighbors_active_set{i}=Q.label_N{r};
                     S.Neighbors_indices_active_set{i}=Q.indices_N{r};
        
                     KNN1=Learn(S);
                     KNN=fitcknn(Train_X(1:sum(S.size_N(1:i)),:),Train_Y(1:sum(S.size_N(1:i))),'NumNeighbors',1,'distance','euclidean');
                     [testPredp, prob] = predict(KNN,Test_X);
                     [testPreda, prob] = predict(KNN1,Test_X);

                     passiveError(i,1)=error_calc(testPredp,Test_Y);
                     passiveError(i,2)=sum(S.size_N(1:i));
                     activeError(i,1)=error_calc(testPreda,Test_Y);
                     activeError(i,2)=sum(S.size_N(1:i));
                     
                     i=i+1;   
        
                 end
                 
                 t=t-t1;        
                 ii=ii+1;
             else
                 ii=ii+1;
             end
        end
        
        KNN=fitcknn(Train_X((1:n),:),Train_Y(1:n),'NumNeighbors',1,'distance','euclidean');%classifier 1NN passive
        [testPredp, prob] = predict(KNN,Test_X);%prédire l'erreur d'entrainement
        passiveError(i-1,1)=error_calc(testPredp,Test_Y);%calcule d'erreur Entrainement
        passiveError(i-1,2)=n;
   
   vr{j,1}=passiveError;
   vr{j,2}=activeError;
   vr{j,3}=S;   
end

for i=1:size(vr,1)
    vv(i)=size(vr{i},1);
end
min_vr=min(vv);
for i=1:size(vr,1)
    ss=vr{i,1};
    vx1=[vx1 ss(1:min_vr,1)];
    vx2=[vx2 ss(1:min_vr,2)];
end
passiveErrorf(:,1)=sum(vx1,2)/size(vr,1);
passiveErrorf(:,2)=sum(vx2,2)/size(vr,1);

for i=1:size(vr,1)
    ss=vr{i,2};
    vx3=[vx3 ss(1:min_vr,1)];
    vx4=[vx4 ss(1:min_vr,2)];
end
activeErrorf(:,1)=sum(vx3,2)/size(vr,1);
activeErrorf(:,2)=sum(vx4,2)/size(vr,1);

ep=std(vx1');
ep1=std(vx3');

standard_diviation_active(alphai,betai)=ep1(size(activeErrorf(:,1),1));     
standard_diviation_passive(alphai,betai)=ep(size(activeErrorf(:,1),1));%
Pr{alphai,betai}=passiveErrorf;
Ar{alphai,betai}=activeErrorf;
Sr{alphai,betai}=S;
rejet(alphai,betai)=cpt;

   end   
end
fprintf('\n');
finalEr=passiveErrorf(i,1);
perf=100-(sum(passiveErrorf(:,1))/(i));
test_finalEr=activeErrorf(i,1);
test_perf=100-(sum(activeErrorf(:,1))/(i));
fprintf('Données d entrainement : \n');
fprintf('taux d erreur final : %f %%\n',finalEr);
fprintf('performance : %f %%\n', perf);
fprintf('Données de test : \n');
fprintf('taux d erreur final : %f %%\n',test_finalEr);
fprintf('performance : %f %%\n', test_perf);

plotactive;
plotpassive;
    


