function ks=minimum(deltachap,deltas)
c=1000;
ks=(c/(deltachap.^2))*(log(1/deltas)+log(log(1/deltas))+log(log(512*sqrt(exp(1))/deltas)));
end