function am=Cij(x,wi,wj,b1,b,ti,tj,dtj,y)
am=-(wi./b1).*exp(-wi.*b.*(ti-x)).*((2*b.*b-1)*sin(b1.*wi.*(ti-x))-2*b.*b1.*cos(b1.*wi.*(ti-x))).*...
    cos(b1.*wj.*(x-tj+dtj)).*exp(-((x-tj+dtj)./y).^2);
end