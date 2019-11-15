function [] = station_data(pre, tspan)

if nargin < 2
  tspan = [];
end

V = SaveStreamData('Read', [pre,'Vel.dat']);
T = SaveStreamData('Read', [pre,'Time.dat']);
T = T / 31556926;


if ~isempty(tspan)
  I = find(T > tspan(1), 1, 'first');
  J = find(T > tspan(2), 1, 'first');
end

semilogy(T(I:J), max(abs(V(:,I:J))))

ind = find_indices(pre, V);
RI = diff(T(ind(2:2:end)));
disp([T(ind(2)), RI(end), RI(end-1)])
