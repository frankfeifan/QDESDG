function plot_field(u, M, G, mesh)

if nargin == 3
  mesh = true;
end

if isempty(u)
  plot(G.x(G.EToV([1,2,4,3,1],:)), G.y(G.EToV([1,2,4,3,1],:)), 'k');
  return
end

trisurf(M.T, M.x, M.y, u, 'LineStyle', 'None')
view(2)
shading interp
axis image
if mesh
  hold on
  plot3(G.x(G.EToV([1,2,4,3,1],:)), G.y(G.EToV([1,2,4,3,1],:)), ...
        ones(size(G.EToV([1,2,4,3,1],:))) * max(u(:)), 'k');
  hold off
end
