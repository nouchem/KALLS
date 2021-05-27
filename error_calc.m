function rate=error_calc(pred,label)

error=0;

for i=1:max(size(label))
    if(pred(i) ~= label(i))
        error=error+1;
    end
end

rate=error*100/max(size(label));

end