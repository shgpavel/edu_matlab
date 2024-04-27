% Pavel Shago V2
%%  Task 2

eq_1 = @(x, y) 2.*x + 8.*y.*x.*sin(y) - 3;
eq_2 = @(x, y) 4 - 3.*x.^2 - 5.*y.^2;

[x, y] = meshgrid(linspace(-10, 10, 400), linspace(-10, 10, 400));

z1 = eq_1(x, y);
z2 = eq_2(x, y);

h = figure;
surf(x, y, z1, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;

surf(x, y, z2, 'FaceColor', 'blue', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
title('Поверхности функций eq_1 и eq_2 на одном графике');
xlabel('x');
ylabel('y');
legend({'eq_1', 'eq_2'});
a = [-4, 2, -100];
plot3(a(1), a(2), a(3), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'black');
view(45, 20);
light;
material shiny;

distance_func = @(w) sqrt((w(1) - a(1)).^2 + (w(2) - a(2)).^2 + (eq_2(w(1), w(2)) - a(3)).^2);
min_distance = fminsearch(@(w) distance_func(w), [5, 2.3]);
plot3([a(1), min_distance(1)], [a(2), min_distance(2)], 
  [a(3), eq_2(min_distance(1), min_distance(2))], 'g--', 'LineWidth', 1);



hold off;
waitfor(h);
