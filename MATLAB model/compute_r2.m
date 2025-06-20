function R2 = compute_r2(x, y, plot_r2)
% Compute the coefficient of determination (R^2) between two variables
% with the fit forced to pass through the origin (0, 0).
%
% Inputs:
%   x - First variable (independent variable), as a vector
%   y - Second variable (dependent variable), as a vector
%
% Output:
%   R2 - The coefficient of determination (R^2)

% Ensure x and y are column vectors
x = x(:);
y = y(:);

% Fit the model with slope only (forced through 0,0)
slope = sum(x .* y) / sum(x .^ 2);
y_fit = slope * x; % Predicted values

% Compute residuals and total variance
SS_res = sum((y - y_fit).^2); % Residual sum of squares
SS_tot = sum((y - mean(y)).^2); % Total sum of squares

% Compute R^2
R2 = 1 - (SS_res / SS_tot);

if plot_r2 == 'y'
% Plot the fit line over the observed range of x
x_fit = linspace(min(x), max(x), 100); % Generate x values for the fit line
y_fit_line = slope * x_fit;            % Compute corresponding y values
plot(x_fit, y_fit_line, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Fit Line (0,0)');
end

% Display the R^2 value
fprintf('R^2 = %.4f\n', R2);


end
