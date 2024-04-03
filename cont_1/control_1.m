% Pavel Shago V1.7
%%  1
%{
eq_1 = @(x, y) y.^2-1.5.*x.*y + x.^2 - 3.*(y + cos(x)) + 3.*(x + sin(2.*y));
eq_2 = @(x, y) -y.^2 + 5.*x.^2 + 3.*x.*y.* cos(x.*y./2) - 7;
eq_comb = @(k) [eq_1(k(1), k(2)), eq_2(k(1), k(2))];
figure;
ezplot(eq_1, [-4 3 -2 4]);
hold on;
grid on;
ezplot(eq_2, [-4 3 -2 4]);
i_ps = [];

text(-3.5, 3.5, {'Выберите точки пересечения кривых',
 'кликая на графике.',
 'Для выхода из выбора нажмите Enter.',}, 'FontSize', 11);
[x, y, button] = ginput();

for i = 1:length(x)
  if button(i) == 13
    break;
  endif

  i_p = fsolve(eq_comb, [x(i), y(i)]);
  i_ps = [i_ps; i_p];
  plot(i_p(1), i_p(2), '*', 'MarkerSize', 10, 'Color', 'r');
  text(i_p(1), i_p(2), sprintf('%.10f, %.10f', i_p(1), i_p(2)),
  'FontSize', 13, 'VerticalAlignment', 'bottom');
end
%}

%%  2
%{
sum = 0;
x = [10 11 3*pi/2];
n = 100000;
format long;
for i = 1:3
  for k = 1:n
    sum += -(1)^(k+1) * sin(k * x(i))/(pi * k);
  end
  disp(sum);
  sum = 0;
end
%}

%%  3

eq_1 = @(theta) sqrt(2.*abs(cos(1.5.*theta)));
eq_2 = @(x, y) x.^2 + 1.1.*x.*y - y.^2 + 5.*y - x - 6;
a = [1, -2];

figure;
hold on;
ezpolar(eq_1, [0, 2*pi]);
grid on;
ezplot(eq_2);
plot(a(1), a(2), 'ro', 'MarkerSize', 5);

dist_y = @(x, y) ((((y-(5/2)).^2) / (1.1*x*y)) - a(2)).^2;
dist_x = @(x, y) ((((x - (1/2)).^2) / (1.1*x*y)) - a(1)).^2;
dist_comb = @(xy) sqrt(dist_x(xy(1), xy(2)) + dist_y(xy(1), xy(2)));
min_distance_point = fminunc(@(xy) dist_comb(xy), [1 1]);
disp(min_distance_point);
plot(min_distance_point(1), min_distance_point(2), 'ro', 'MarkerSize', 5);
plot([a(1), min_distance_point(1)], [a(2), min_distance_point(2)], 'g--', 'LineWidth', 1);

curve_points_x = @(theta) eq_1(theta) .* cos(theta);
curve_points_y = @(theta) eq_1(theta) .* sin(theta);
distances = @(theta) sqrt((curve_points_x(theta) - a(1)).^2 + (curve_points_y(theta) - a(2)).^2);
min_dist = fminunc(@(theta) distances(theta), 0);
min_dist_p = [curve_points_x(min_dist), curve_points_y(min_dist)];

plot(min_dist_p(1), min_dist_p(2), 'ro', 'MarkerSize', 5);
plot([a(1), min_dist_p(1)], [a(2), min_dist_p(2)], 'g--', 'LineWidth', 1);
text_p = [(a(1) + min_dist_p(1)) / 2; (a(2) + min_dist_p(2)) / 2];
text(text_p(1), text_p(2), ['Мин. расстояние: ', num2str(min_distance)], 'HorizontalAlignment', 'center');

