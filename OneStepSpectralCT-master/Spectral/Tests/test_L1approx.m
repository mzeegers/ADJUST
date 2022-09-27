x = linspace(-2,2, 10000);
delta = 0.1;
x_spacing = x(2) - x(1);

% % Test Long2014's hyperbola
% h_x = hyperbola(x, delta, 0);
% h_dot_x = hyperbola(x, delta, 1);
% h_dotdot_x = hyperbola(x, delta, 2);
% 
% figure(1); plot(x, h_x);
% figure(2); plot(x, h_dot_x, 'o', x, gradient(h_x) / x_spacing, '+');
% figure(3); plot(x, h_dotdot_x, 'o', x, gradient(gradient(h_x)) / x_spacing^2, '+');
% 
% % Test Huber function
% h_x = Huber(x, delta, 0);
% h_dot_x = Huber(x, delta, 1);
% h_dotdot_x = Huber(x, delta, 2);
% 
% figure(4); plot(x, h_x);
% figure(5); plot(x, h_dot_x, 'o', x, gradient(h_x) / x_spacing, '+');
% figure(6); plot(x, h_dotdot_x, 'o', x, gradient(gradient(h_x)) / x_spacing^2, '+');

% Test Green prior
h_x = GreenPrior(x, 0);
h_dot_x = GreenPrior(x, 1);
h_dotdot_x = GreenPrior(x, 2);

figure(7); plot(x, h_x);
figure(8); plot(x, h_dot_x, 'o', x, gradient(h_x) / x_spacing, '+');
figure(9); plot(x, h_dotdot_x, 'o', x, gradient(gradient(h_x)) / x_spacing^2, '+');
