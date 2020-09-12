function [ h,n,m ] = bivariate_Gaussian( dim,mu,sig2,theta )

% creates two-dimensional Gaussian filters

%%%%%%%%%%%%%%%%%   input   %%%%%%%%%%%%%%%%%%%
% dim = size of filter [N M]
% mu = mean of filter []
% sig2

%%%%%%%%%%%%%%%%   output   %%%%%%%%%%%%%%%%%%
% h = filter
% n = dimension 1 axis
% m = dimension 2 axis

% written by Brad Bazow 09/2018

M = dim( 1 ); % length of filter dimension 1
N = dim( 2 ); % length of filter dimension 1
mu_x = mu( 1 ); % mean of dimension 1
mu_y = mu( 2 ); % mean of dimension 2
sig2_x = sig2( 1 ); % variance of dimension 1
sig2_y = sig2( 2 ); % variance of dimension 2

% if dimension 1 > 1 -> form axis
if N > 1
    n = linspace( -( N-1 )/2,( N-1 )/2,N ).';
else
    n = 0;
end
% if dimension 2 > 1 -> form axis
if M > 1
    m = linspace( -( M-1 )/2,( M-1 )/2,M ).';
else
    m = 0;
end

[ y,x ] = meshgrid( n,m ); % create tow dimensional axes
xx = x*cosd( theta )-y*sind( theta );
yy = x*sind( theta )+y*cosd( theta );
% create Gaussian
h = exp( -( ( xx-mu_x ).^2/( 2*sig2_x ) + ( yy-mu_y ).^2/( 2*sig2_y ) ) );
h = h/sum( h( : ) ); % normalize for unity norm

end