function [f g] = funLS2(r, params)
f = norm(r,2);
%g = r./f;
g = r;

end