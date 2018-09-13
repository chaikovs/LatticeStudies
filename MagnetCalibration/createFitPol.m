function [fitresult, gof] = createFitPol(I, current_mean)
%CREATEFIT(I,CURRENT_MEAN)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : I
%      Y Output: current_mean
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 21-Nov-2016 17:23:52


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( I, current_mean );

% Set up fittype and options.
ft = fittype( 'poly4' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Create a figure for the plots.
figure( 'Name', 'untitled fit 1' );

% Plot fit with data.
subplot( 2, 1, 1 );
h = plot( fitresult, xData, yData,'' );
legend( h, 'current_mean vs. I', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel I
ylabel current_mean
grid on

% Plot residuals.
subplot( 2, 1, 2 );
h = plot( fitresult, xData, yData, 'residuals' );
legend( h, 'untitled fit 1 - residuals', 'Zero Line', 'Location', 'NorthEast' );
% Label axes
xlabel I
ylabel current_mean
grid on


