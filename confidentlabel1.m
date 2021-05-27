function [rcoord,rlabel,Y,t,r]=confidentlabel1(x,ks,deltas,valCoord,valLabel,S,idx,N)

i=1;v1=[];v2=[];r=0;s=0;l=[];
v1=S.Neighbors_indices_active_set{1,idx};
v2=S.Neighbors_active_set{1,idx};
l=S.label_Neighbors_active_set{1,idx};
rcoord=[];rlabel=[];

if isempty(find(v1==N(i)))==0
    rcoord(i,:)=v2(find(v1==N(i)),:);
    rlabel(i)=l(find(v1==N(i)));
    r=r+1;
else
    rcoord(i,:)=valCoord(i,:);
    rlabel(i)=valLabel(i);
    s=s+1;
end
i=i+1;
while i<=min(ks,100)
    Bs=sqrt((2/i)*(log(1/deltas)+log(log(1/deltas))+log(log(exp(1)*i))));
    if abs(((1/i)*sumlabel(rlabel))-1/2)>2*Bs 
        break;
    else
        if isempty(find(v1==N(i)))==0
            rcoord(i,:)=v2(find(v1==N(i)),:);
            rlabel(i)=l(find(v1==N(i)));
            r=r+1;
        else
            rcoord(i,:)=valCoord(i,:);
            rlabel(i)=valLabel(i);
            s=s+1;
        end
        i=i+1;
    end
end
v=var(rlabel);
v1=v;
while i<=ks && v1<=0.2                         
    
    Bs=sqrt((2/i)*(log(1/deltas)+log(log(1/deltas))+log(log(exp(1)*i))));            
    if abs(((1/i)*sumlabel(rlabel))-1/2)>2*Bs             
        break;            
    else
        if isempty(find(v1==N(i)))==0            
            rcoord(i,:)=v2(find(v1==N(i)),:);            
            rlabel(i)=l(find(v1==N(i)));            
            r=r+1;            
        else
            rcoord(i,:)=valCoord(i,:);           
            rlabel(i)=valLabel(i);            
            s=s+1;           
        end
        i=i+1;              
    end
    v1=var(rlabel);          
end

Y=calcullabelF(rlabel);
  
end