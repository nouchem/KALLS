
function [rcoord,rlabel,Y]=confidentlabel(Index,ks,deltas,valCoord,valLabel,kp)
i=1;
rcoord=[];rlabel=[];
rcoord(i,:)=valCoord(i,:);
rlabel(i)=valLabel(i);
i=i+1;
while i<=min(100,ks)
    Bs=sqrt((2/i)*(log(1/deltas)+log(log(1/deltas))+log(log(exp(1)*i))));
    if abs(((1/i)*sumlabel(rlabel))-1/2)>2*Bs 
        break;
    else        
        rcoord(i,:)=valCoord(i,:);
        rlabel(i)=valLabel(i);      
        i=i+1;         
    end    
end

v=var(rlabel);
v1=v;          

while i<=ks && v1<0.2 
    Bs=sqrt((2/i)*(log(1/deltas)+log(log(1/deltas))+log(log(exp(1)*i))));  
    if abs(((1/i)*sumlabel(rlabel))-1/2)>2*Bs    
        break;      
    else
        rcoord(i,:)=valCoord(i,:);
        rlabel(i)=valLabel(i);           
        i=i+1;    
    end
    v1=var(rlabel);   
end

Y=calcullabel(rlabel);
   
end
