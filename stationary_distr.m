P = [0.977 0.023; 0.074 0.926];
[V , D] = eig(P')
V(: , 1)/sum(V(: , 1))