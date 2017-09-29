function v=gausfunc(a,b,c,x)

v = a.*exp( -( (x-b).^2 / (2.*c.^2) ) );
end