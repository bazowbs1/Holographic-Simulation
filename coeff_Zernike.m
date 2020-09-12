function [coeff_Zer, ph] = coeff_Zernike(coeff,ph)
     % this code is to calculate the Zernike coefficients from 2D polynomial
     % matlab fitting with 5 order in x and y
     % coeff_Zer = array of 21 Zernike coefficients
     % ph is combination phase of 21 orders corresponding to 21 Zernike coefficients 
     %% Calculate Zernike coeff
   
Z = [1 0 0 -sqrt(3) 0 0 0 0 0 0 sqrt(5) 0 0 0 0 0 0 0 0 0 0;                             % P00
     0 2 0 0 0 0 0 -2*sqrt(8) 0 0 0 0 0 0 0 3*sqrt(12) 0 0 0 0 0;                        % P10
     0 0 2 0 0 0 -2*sqrt(8) 0 0 0 0 0 0 0 0 0 3*sqrt(12) 0 0 0 0;                        % P01
     0 0 0 2*sqrt(3) 0 sqrt(6) 0 0 0 0 -6*sqrt(5) -3*sqrt(10) 0 0 0 0 0 0 0 0 0;         % P20
     0 0 0 0 2*sqrt(6) 0 0 0 0 0 0 0 -6*sqrt(10) 0 0 0 0 0 0 0 0;                        % P11
     0 0 0 2*sqrt(3) 0 -sqrt(6) 0 0 0 0 -6*sqrt(5) 3*sqrt(10) 0 0 0 0 0 0 0 0 0;         % P02
     0 0 0 0 0 0 0 3*sqrt(8) 0 sqrt(8) 0 0 0 0 0 -12*sqrt(12) 0 -4*sqrt(12) 0 0 0;       % P30
     0 0 0 0 0 0 3*sqrt(8) 0 3*sqrt(8) 0 0 0 0 0 0 0 -12*sqrt(12) 0 -12*sqrt(12) 0 0;    % P21
     0 0 0 0 0 0 0 3*sqrt(8) 0 -3*sqrt(8) 0 0 0 0 0 -12*sqrt(12) 0 12*sqrt(12) 0 0 0;    % P12
     0 0 0 0 0 0 3*sqrt(8) 0 -sqrt(8) 0 0 0 0 0 0 0 -12*sqrt(12) 0 4*sqrt(12) 0 0;       % P03
     0 0 0 0 0 0 0 0 0 0 6*sqrt(5) 4*sqrt(10) 0 sqrt(10) 0 0 0 0 0 0 0;                  % P40
     0 0 0 0 0 0 0 0 0 0 0 0 8*sqrt(10) 0 4*sqrt(10) 0 0 0 0 0 0;                        % P31
     0 0 0 0 0 0 0 0 0 0 12*sqrt(5) 0 0 -6*sqrt(10) 0 0 0 0 0 0 0;                       % P22
     0 0 0 0 0 0 0 0 0 0 0 0 8*sqrt(10) 0 -4*sqrt(10) 0 0 0 0 0 0;                       % P13
     0 0 0 0 0 0 0 0 0 0 6*sqrt(5) -4*sqrt(10) 0 sqrt(10) 0 0 0 0 0 0 0;                 % P04
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10*sqrt(12) 0 5*sqrt(12) 0 sqrt(12) 0;                % P50
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10*sqrt(12) 0 15*sqrt(12) 0 5*sqrt(12);             % P41
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20*sqrt(12) 0 -10*sqrt(12) 0 -10*sqrt(12) 0;          % P32
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20*sqrt(12) 0 10*sqrt(12) 0 -10*sqrt(12);           % P23    
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10*sqrt(12) 0 -15*sqrt(12) 0 5*sqrt(12) 0;            % P14
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10*sqrt(12) 0 -5*sqrt(12) 0 sqrt(12);               % P05
     ];
     %% Cartesian From up to 5th Order
    [x,y] = size(ph);
    xx = (1:1:x); yy = (1:1:y);
    %xx = xx-xx(round( size( ph,2 ) )/2)-1;yy = yy-yy(round( size( ph,1 ) )/2)-1;
    [X,Y] = meshgrid(xx,yy);
    Z00 = ones(y,x);
    Z01 = 2*X;
    Z02 = 2*Y;
    Z03 = sqrt(3)*(2*X.^2+2*Y.^2-1);
    Z04 = sqrt(6)*(2*X.*Y);
    Z05 = sqrt(6)*(X.^2-Y.^2);
    Z06 = sqrt(8)*(3*X.^2.*Y+3*Y.^3-2*Y);
    Z07 = sqrt(8)*(3*X.^3+3*X.*(Y.^2)-2*X);
    Z08 = sqrt(8)*(3*X.^2.*Y-Y.^3);
    Z09 = sqrt(8)*(X.^3-3*X.*(Y.^2));
    Z10 = sqrt(5)*(6*X.^4+12*(X.^2).*(Y.^2)+6*Y.^4-6*X.^2-6*Y.^2+1);
    Z11 = sqrt(10)*(4*X.^4-3*X.^2+3*Y.^2-4*Y.^4);
    Z12 = sqrt(10)*(8*X.^3.*Y+8*X.*(Y.^3)-6*X.*Y);
    Z13 = sqrt(10)*(X.^4-6*(X.^2).*(Y.^2)+Y.^4);
    Z14 = sqrt(10)*(4*X.^3.*Y-4*X.*(Y.^3));
    Z15 = sqrt(12)*(10*X.^5+20*(X.^3).*(Y.^2)+10*X.*(Y.^4)-12*X.^3 - 12*X.*(Y.^2)+3*X);
    Z16 = sqrt(12)*(10*Y.^5+20*(Y.^3).*(X.^2)+10*(X.^4).*Y-12*Y.^3 - 12*(X.^2).*Y+3*Y);
    Z17 = sqrt(12)*(5*X.^5-10*(X.^3).*(Y.^2)-15*X.*(Y.^4)-4*(X.^3)+12*X.*(Y.^2));
    Z18 = sqrt(12)*(15*(X.^4).*Y+10*(X.^2).*(Y.^3)-5*Y.^5-12*(X.^2).*Y+4*(Y.^3));
    Z19 = sqrt(12)*(X.^5-10*(X.^3).*(Y.^2)+5.*X.*(Y.^4));
    Z20 = sqrt(12)*(5*(X.^4).*Y-10*(X.^2).*(Y.^3)+Y.^5);
    Zer = [];
    Zer(:,:,1)  = Z00; Zer(:,:,2)  = Z01; Zer(:,:,3)  = Z02; Zer(:,:,4)  = Z03;
    Zer(:,:,5)  = Z04; Zer(:,:,6)  = Z05; Zer(:,:,7)  = Z06; Zer(:,:,8)  = Z07;
    Zer(:,:,9)  = Z08; Zer(:,:,10) = Z09; Zer(:,:,11) = Z10; Zer(:,:,12) = Z11;
    Zer(:,:,13) = Z12; Zer(:,:,14) = Z13; Zer(:,:,15) = Z14; Zer(:,:,16) = Z15;
    Zer(:,:,17) = Z16; Zer(:,:,18) = Z17; Zer(:,:,19) = Z18; Zer(:,:,20) = Z19;
    Zer(:,:,21) = Z20;
    coeff_Zer = Z( 1:length( coeff ),1:length( coeff ) )\coeff; %% the Zernike coefficient

    %% generate Aberration from Zerniki
    ph = zeros(y,x);
    for i = 1:length( coeff )
        ph = ph+coeff_Zer(i)*Zer(:,:,i); % summing of aberration
    end
    ph = ph;
    %figure; surf(ph); axis square; colorbar; shading interp; title('Aberration - Zernike');
end