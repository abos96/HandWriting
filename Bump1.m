function y = Bump1(x,ux,nx)

y = exp(-1./(nx.^2-(x-ux).^2)).*(abs(x-ux) <= nx); 

