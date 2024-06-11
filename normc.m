function y = normc(x)
% this function normalizes entries in x by the magnitude of each column 
% y is the same size as x but contains entries normalized by the magnitude of
% each column

[r,c] = size(x);
for ii=1:c,
    y(:,ii) = x(:,ii)/norm(x(:,ii));
end