function plotroutine(n)
%
%  plotroutine(A)
%
%  Plots a two-dimensional world represented by the matrix A 
%  where nonzero entries represent living cells. 
%
[m,n]=size(n);
[x, y] = find(n);
hold on
for K = 1 : length(x)
  rectangle('Position', [y(K)-1, abs(x(K)-m) 1, 1], 'Linewidth', 2, 'facecolor', 'k','EdgeColor','k');
end

for j=0:n
  plot([j j],[0 m],'k-','Linewidth', 2),
end
for i=0:m
  plot([0 n],[i i],'k-','Linewidth', 2)
end
axis equal
