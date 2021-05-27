%Regression function  "1/2*[1-((sin(2*pi*x)sin(2*pi*y))]
        
function s=regfcturner(x,y)


k=sin(2.*pi.*x).*sin(2.*pi.*y);

s=(1/2).*(1.-k);
end