% Pavel Shago V2
%%  Task 2

%%  1

%{
hold on;
grid on;
xlabel('x');
ylabel('y');

f = @(x) (0.5 - x).*(x <= 0.5) + cos(pi*x).*(x > 0.5);
w = @(x) x.^2;
T = @(mn, x) sin(mn * x) ./ x;

a = 0;
b = 1;
i = 0;
n_values = [5, 10, 20, 50];

for n=n_values
  sum = @(x) 0;
  for i = 1:n
    mu_x=i;
    mu = @(x) tan(x);
    mu_guess = (mu_x-0.5) .*pi;
    mu_n = fsolve(mu, mu_guess);

    numerator = integral(@(x) T(mu_n, x) .* f(x) .* w(x), a, b, 'ArrayValued', true);
    denominator = integral(@(x) (T(mu_n, x).^2) .* w(x), a, b, 'ArrayValued', true);

    a_i = numerator/denominator;
    prev_sum = sum;
    sum = @(x) prev_sum(x)+a_i .*(sin(mu_n .* x)./x);
  end
  handle = ezplot(sum, [0,1]);
  set(handle,'LineWidth',3);
end

handle = ezplot(f, [0,1]);
set(handle,'LineWidth',4);
legend('F_5(x)', 'F_{10}(x)', 'F_{20}(x)', 'F_{50}(x)', 'f(x)');
hold off;
%}


%%  2
%{
eq_1 = @(x, y) 2.*x + 8.*y.*x.*sin(y) - 3;
eq_2 = @(x, y) 4 - 3.*x.^2 - 5.*y.^2;

[x, y] = meshgrid(linspace(-5, 5, 50));

z1 = eq_1(x, y);
z2 = eq_2(x, y);

h = figure;
surf(x, y, z1, 'FaceColor', 'yellow', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;

surf(x, y, z2, 'FaceColor', 'blue', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
title('Поверхности функций eq_1 и eq_2');
xlabel('x');
ylabel('y');
legend({'eq_1', 'eq_2'});
a = [-4, 2, -100];
plot3(a(1), a(2), a(3), 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'black');
view(-45, 20);
light;
material shiny;

distance_func = @(w) sqrt((a(1) - w(1)).^2 + (a(2) - w(2)).^2 + (a(3) - eq_2(w(1), w(2))).^2);
min_distance = fminunc(@(w) distance_func(w), [100, -100]);
b = [min_distance(1), min_distance(2), eq_2(min_distance(1), min_distance(2))];

plot3([a(1), b(1)], [a(2), b(2)], [a(3), b(3)], 'g--', 'LineWidth', 1);
plot3(b(1), b(2), b(3), 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'black');

difference = z1 - z2;
threshold = 0.01;
[C, ~] = contourc(x, y, difference, [threshold, -threshold]);

for i = 1:length(C)
    if C(1, i) == threshold || C(1, i) == -threshold
        j = i + 1;
        num_points = C(2, i);
        contour_x = C(1, j:j+num_points-1);
        contour_y = C(2, j:j+num_points-1);

        contour_z1 = interp2(x, y, z1, contour_x, contour_y);
        contour_z2 = interp2(x, y, z2, contour_x, contour_y);

        contour_z = (contour_z1 + contour_z2) / 2;

        plot3(contour_x, contour_y, contour_z, 'r', 'LineWidth', 2);
    end
end

[~, idx_x] = min(abs(x(1, :) - b(1)));
[~, idx_y] = min(abs(y(:, 1) - b(2)));

[dx, dy] = gradient(z2, x(1, 2) - x(1,1), y(2, 1) - y(1, 1));
dzdx = dx(idx_y, idx_x);
dzdy = dy(idx_y, idx_x);

normal_b = [-dzdx, -dzdy, 1];
normal_b = [normal_b(1), normal_b(2), normal_b(3)] / norm(normal_b);

quiver3(b(1), b(2), b(3), normal_b(1), normal_b(2), normal_b(3), 'r', 'LineWidth', 2);

ab = b - a;
length_normal_b = norm(normal_b);
length_ab = norm(ab);

cos_angle = dot(normal_b, ab) / (length_normal_b * length_ab);

angle_deg = rad2deg(acos(cos_angle));
annotation('textbox', [0.18, 0.78, 0.1, 0.1], 'String', ['Angle ', num2str(angle_deg)],
           'EdgeColor', 'none', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'red');


hold off;
waitfor(h);
%}

%% 3
%%{
f = @(y) sin(y) ./ log(y + 1) - 3 * cos(y) ./ (2 * sqrt(y));
F = @(x) arrayfun(@(x) integral(f, 0, x), x);

x = linspace(0, 10, 100);
F_values = arrayfun(F, x);

x_vals = linspace(0, 10, 100);
zero_points = [];

for i = 1:numel(x_vals) - 1
  x0 = x_vals(i);
  x1 = x_vals(i + 1);
  try
    zero = fzero(F, [x0, x1]);
    zero_points = [zero_points, zero];
  catch
  end
end

disp('Нули F(x)');
disp(num2str(zero_points));

dFdx = diff(F_values) ./ diff(x);
x_mid = (x(1:end-1) + x(2:end)) / 2;

extremums = [];
for i = 2:numel(dFdx)
  if dFdx(i - 1) * dFdx(i) <= 0
    extremum_x = interp1(dFdx([i-1, i]), x_mid([i-1, i]), 0);
    extremums = [extremums, extremum_x];
  end
end

disp('Экстремумы F(x)');
disp(num2str(extremums));

h = figure;
plot(x, F_values, 'b-', 'LineWidth', 4);
hold on;
grid on;

plot(zero_points, zeros(size(zero_points)), 'ro', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', 'Нули F(x)');

extremum_y = arrayfun(F, extremums);
plot(extremums, extremum_y, 'bo', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', 'Экстремумы F(x)');

legend('F(x)', 'Нули F(x)', 'Экстремумы F(x)');
xlabel('x');
ylabel('F(x)');
title('График функции F(x) с отмеченными нулями и экстремумами');
grid on;
hold off;
waitfor(h);
%%}

