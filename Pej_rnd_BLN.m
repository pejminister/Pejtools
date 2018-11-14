function x = Pej_rnd_BLN(N, mu, v)
if size(N,2)>1
    error('N is assumed to be a nx1 column')
end


if numel(mu)==1
    mu = ones(size(N))*mu;
end
if numel(v)==1
    v = ones(size(N))*v;
end

t = mu + randn(size(N)).*sqrt(v);
r = Pej_Transform_logistic(t);
x = binornd(N,r);
end
