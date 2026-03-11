function [p1,p2,n1,n2,r] = match_points(SP,TP,SN,TN,Btree)

idx = knnsearch(Btree, SP');

p1 = SP;
n1 = SN;

p2 = TP(:,idx);
n2 = TN(:,idx);

d = p1 - p2;
dist = dot(d, n1 + n2)';

r = abs(dist);
end
