function y=calcullabel(rlabel)

nchap=(1/max(size(rlabel)))*sumlabel(rlabel);

if nchap>=1/2 
    y=1;
else
    y=0;
end

end