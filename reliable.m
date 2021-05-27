function T=reliable(x,alpha,L,S,d,N,w)

c=4;
for i=1:size(S.coord_active_set,1)
    v(i)=size(N(N<=norm(S.coord_active_set(i,:)-x)),2)/w;
end
T=1;

for i=1:size(S.coord_active_set,1)
    if v(i)<=((c/L)*S.LB(i))^(d/alpha)
        T=0;
        break;
    end
end

end
  
