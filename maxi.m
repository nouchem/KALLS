function a=maxi(epsilon,betta,c)

if epsilon/2>=(epsilon/(2*c)).^(1/(betta+1))
    a=epsilon/2;
else
    a=(epsilon/(2*c)).^(1/(betta+1));
end
end