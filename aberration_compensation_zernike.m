function [ phase_zernike ] = aberration_compensation_zernike( phase0,X,Y,num )

Xdata = []; Ydata = []; ROI_Vec = [];
figure; imagesc( 180/pi*phase0 );colorbar;axis xy 
title(['Please choose ',num2str( num ),' regions of flat phase by double left click']);

for i=1:num
    h = imrect;
    vertices = wait( h );
    xmin = round( vertices( 1 ) );ymin = round( vertices( 2 ) );width = round( vertices( 3 ) );height = round( vertices( 4 ) );
    X = [ xmin xmin+width ];
    Y = [ ymin ymin+height ];
    
    xx = X( 1 ):X( 2 );
    yy = Y( 1 ):Y( 2 );
    [ XX,YY ] = meshgrid( xx,yy );
    Xdata = [ Xdata;XX( : ) ];Ydata = [ Ydata;YY( : ) ];
    
    ROI = phase0( Y( 1 ):Y( 2 ),X( 1 ):X( 2 ) );
    ROI_Vec = [ ROI_Vec;ROI( : ) ];
end

surface_reg = fit( [ Xdata,Ydata ],ROI_Vec,'poly55' );
coeff_reg = ( coeffvalues( surface_reg ) )'; %return the values of coefficients

[ x_tmp,y_tmp ] = size( phase0 );
xx = 1:x_tmp;yy = 1:y_tmp;
[ X,Y ] = meshgrid( xx,yy );

X1D = X( : ); Y1D = Y( : );
sf = zeros( length( X1D ),1 );
for ii = 1:length( X1D )
    if mod( ii,100000 ) == 0
        disp( ii )
    end
    sf( ii ) = surface_reg( Y1D( ii ),X1D( ii ) );
end
sf = ( reshape( sf,y_tmp,x_tmp ) )';
[ coeff_zernike,phase_zernike ] = coeff_Zernike( coeff_reg,sf );

end
