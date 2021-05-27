
%This code is an implementation of the K-nn Active Learning Under Local Smoothness (KALLS). For theoritical details see the publication: "K-NN active learning under local smoothness assumption" by  B. N. Njike and X. Siebert

%Input:
%  epsilon: accuracy parameter
%  delta: Confidence parameter
%  beta,c: margin noise parameters
%  alpha,L: smoothness parameters
%  n:label budget   
%Output:
%  activeErrorf
%  passiveErrorf

close all
clear all
clc
Synth_data=1;
parameter=1;

passiveError=[];
activeError=[];
size_train_data=2*(10^6);%size of training data
size_test_data=10^4;%size of testing data


Q={};
I={};
S={};
T=1;
r=0;
i=1;

switch Synth_data
    case 1 %Uniform distribution
        rng('shuffle');
        data1=-1 + (0.5-(-1)).*rand(size_train_data,1); 
        data2=-1 + (0.5-(-1)).*rand(size_train_data,1);
        Train_X=[data1 data2];
        Train_X=round(Train_X,4);
        sl=regfcturner(data1,data2);
        sl=sl';
        Train_Y=binornd(1,sl);
        clear data1 data2 

        data1=-1 + (0.5-(-1)).*rand(size_test_data,1); 
        data2=-1 + (0.5-(-1)).*rand(size_test_data,1);
        Test_X=[data1 data2];
        Test_X=round(Test_X,4);
        sl1=regfcturner(data1,data2);
        sl1=sl1';
        Test_Y=binornd(1,sl1);
        clear data1 data2 sl1

    case 2 %Normal distribution
 
        d=2;
        mu=zeros(d,1);
        sigma=eye(d);
        Train_X = round(mvnrnd(mu,sigma,size_train_data),4);
        Test_X= round(mvnrnd(mu,sigma,size_test_data),4);
        a=[1,-1,1,-1,1,-1];
        a1=a(1:d);
        v=sqrt(sum((Train_X-a1).^2,2));
        v1=sqrt(sum((Test_X-a1).^2,2));
        train1=Train_X(v>2,:);
        sl1=(1/5).*ones(size(train1,1),1);
        test1=Test_X(v1>2,:);
        sl11=(1/5).*ones(size(test1,1),1);
        train2=Train_X(v<=2,:);
        sl2=1.-(2/5.*v(v<=2));
        test2=Test_X(v1<=2,:);
        sl12=1.-(2/5.*v1(v1<=2));
        
        Train_X=[train1;train2];
        sl=[sl1;sl2];
        Train_Y=binornd(1,sl)';

        Test_X=[test1;test2];
        sl3=[sl11;sl12];
        Test_Y=binornd(1,sl3)';

end
        
for x=1:10
perm=randperm(size(Train_X,1));
Train_X=Train_X(perm,:);
Train_Y=Train_Y(perm);

Train_tree = createns(Train_X,'nsmethod','kdtree');      



ii=1;
switch parameter
    case 1 %uniform distribution
        epsilon=0.09;%accuracy parameter
        delta=0.1;%Confidence parameter
        beta=1;%margin noise parameter
        c=15;%margin noise parameter
        L=4*pi;%smoothness parameter
        alpha=1;%smoothness parameter
        n=4*(10^5);%label budget
        d=size(Train_X,2);
        cpt=0;
    case 2 %Normal distribution
        beta=1;
        alpha=1;
        n=4*(10^5);
        epsilon=0.05;
        delta=0.15;
        for h=1:size(Train_X,1)
            l(h)=norm(Train_X(h,:));
        end
        f1=find(l==max(l));
        m1=(1/(2*pi*sqrt(det(sigma))))*exp(-1/2*(Train_X(f1,:)'-mu)'*sigma*(Train_X(f1,:)'-mu));
        clear l
        for h=1:size(Train_X,1)
            l(h)=norm(Train_X(h,:)-a1);
        end
        l=l(l<=2);
        coords=Train_X(l<=2);
        f2=find(l==max(l));
        m2=(1/(2*pi*sqrt(det(sigma))))*exp(-1/2*(Train_X(f2,:)'-mu)'*sigma*(Train_X(f2,:)'-mu));
        clear l coords
        r_a=1;
        r_2=2;
        if mod(d,2)==0
        v=(pi^(d/2)/factorial((d/2)))*(r_a^d);
        v1=(pi^(d/2)/factorial((d/2)))*(r_2^d);
        else
        v=(pi^((d-1)/2)/1*3)*(r_a^d)*(2^((d+1)/2));
        v1=(pi^((d-1)/2)/1*3)*(r_2^d)*(2^((d+1)/2));
        end    
        k=(pi^(d/2)/factorial((d/2)));
        f_i=2/5;
        c=2^(beta)*8;
        L=7*log(1/epsilon);
       
end
 
deltachap=maxi(epsilon,beta,c);
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
        S.Neighbors_active_set{i}=Q.coord_N{r};
        S.label_Neighbors_active_set{i}=Q.label_N{r};
        S.Neighbors_indices_active_set{i}=Q.indices_N{r};
        S.size_N(i)=I.size_N(r);      
        
        KNN1=Learn(S);
        
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
             [ind1,N1] = knnsearch(Train_tree,Train_X(ii,:),'K',size(Train_X,1));
             N1(N1==0)=[]; 
             ind(ind==Index)=[];
   
             T=reliable(Train_X(ii,:),alpha,L,S,d, N1,size(Train_X,1));
       
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
         
                     [rcoord,rlabel,Y]=confidentlabel(ii,ks,deltas,Train_X(ind,:),Train_Y(ind));
                     t1=size(rcoord,1);
                 else
                     [rcoord,rlabel,Y,t1,s1]=confidentlabelF1(ii,ks,deltas,Train_X(ind,:),Train_Y(ind),S,i_1,ind);
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
        [testPredp, prob] = predict(KNN,Test_X);
        passiveError(i-1,1)=error_calc(testPredp,Test_Y);
        passiveError(i-1,2)=n;
end

fprintf('\n');
test_Error=activeError(i,1);
test_perf=100-(sum(activeErrorf(:,1))/(i));
fprintf('Test data : \n');
fprintf('Error rate : %f %%\n',test_Error);
fprintf('performance : %f %%\n', test_perf);

h1=plot(activeErrorf(:,2),activeErrorf(:,1),'b');        
hold on        
h2=plot(passiveErrorf(:,2),passiveErrorf(:,1),'r');       
hold on       
legend([h1,h2],'KALLS','1-NN-standard');       
hold on       
ylabel('Error(%)') % x-axis label
xlabel('Number of learned examples')
hold on

    


