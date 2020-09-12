%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation of holographic recording and reconstruction using an off-axis
% telecentric configuration. Fresnel diffraction or the angular spectrum
% method is used for computational diffraction.
%
% written by Brad Bazow ( September, 2020 )
% The Johns Hopkins University Applied Physics Labratory ( JHU/APL )
% The Catholic University of America Electrical and Computer Science
% Department ( CUA/EECS )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; 

diffraction = 'angular spectrum'; % 'fresnel', 'angular spectrum'
lowpass_filter = 'gaussian'; % 'ideal', 'butterworth', 'gaussian'

%% parameters
lambda = 632.8e-9; % wavelength ( m )
k0 = 2*pi/lambda; % wavenumber ( rad/m )
pixel_pitch = 7.1250e-06; % pixel spacing ( m/pixel )
us_rate = 4; % image upsample rate
magnification = 3; % magnification of telecentric configuration
f1 = 16.5e-2; % focal length of ideal lens 1 ( m )
f2 = magnification*f1; % focal length of ideal lens 2 ( m )
I0 = 1; % maximum intensity of complex field object
hz0 = 150e-9; % maximum thickness of complex field object
k = 1.2126e-8; % extinction coefficient for BK7 optical glass
alpha = 4*pi*k/lambda; % absorption coefficient for BK7 optical glass ( 1/m )
n0 = 1; % refractive index of medium
n1 = 1.5151; % refractive index of object

%% generate USAF 1951 phase object

zlims = [ 0 hz0 ]; % min and max height profile

z = double( rgb2gray( imread( 'usaf_1951.png' ) ) );
z = z/max( z( : ) );
z = flipud( abs( z-1 ) );
z = ( zlims( 2 )-zlims( 1 ) )*( z-min( z( : ) ) )/( max( z( : ) )-min( z( : ) ) ) + zlims( 1 );

[ M,N ] = size( z ); % image dimensions

% zero pad to force square image
if M<N
    z = cat( 1,zlims( 1 )*ones( N-M,N ),z );
elseif N<M
    z = cat( 2,zlims( 1 )*ones( M,M-N ),z);
end
z = z( 1:364,1:364 ); % crop image
[ M,N ] = size( z ); % image dimensions

x = 0:N-1; y = 0:M-1;
[ X,Y ] = meshgrid( x,y );

% interpolate onto upsampled grid
upsample_rate = [ 4 4 ];
xq = linspace( 0,N-1,upsample_rate( 2 )*N );
yq = linspace( 0,M-1,upsample_rate( 1 )*M );
[ Xq,Yq ] = meshgrid( xq,yq );
z = interp2( X,Y,z,Xq,Yq );

opd = ( n1-n0 )*z; % optical path length
image_phase = 2*pi*opd/lambda;

T = exp( -1*alpha*opd ); % transmittance via Beer's law
image_intensity = T*I0;

im = sqrt( image_intensity ).*exp( 1i*image_phase );

%% Predict object field at the recording plane

[ M,N ] = size( im );

image_intensity = im.*conj( im );
image_phase = atan2( imag( im ),real( im ) );

K = us_rate*max( M,N ); % object dimensions

% zero pad image
pad1 = zeros( K,( K-N )/2 );
pad2 = zeros( ( K-M )/2,N );
psi0 = [ pad1,[ pad2;im;pad2 ],pad1 ]; % object

object_intensity = psi0.*conj( psi0 );
object_phase = atan2( imag( psi0 ),real( psi0 ) );

%% propagate object to lens 1

m = ( 0:K-1 )-K/2;
n = ( 0:K-1 )-K/2;
kx = 2*pi*m/( us_rate*M*pixel_pitch );
ky = 2*pi*n/( us_rate*N*pixel_pitch );
[ Kx,Ky ] = meshgrid( kx,ky );
% object field in Fourier domain
Psi0 = fftshift( fft2( fftshift( psi0 ) ) );

switch diffraction
    case 'fresnel'
        H0 = exp( -1i*k0*f1 )*exp( 1i*f1*( Kx.^2+Ky.^2 )/( 2*k0 ) );
        psi = fftshift( ifft2( fftshift( Psi0.*H0 ) ) );
    case 'angular spectrum'
        H0 = exp( -1i*k0*f1*sqrt( 1-( Kx/k0 ).^2-( Ky/k0 ).^2 ) );
        psi = fftshift( ifft2( fftshift( Psi0.*H0 ) ) );
end

% apply phase transfer function of lens 1
x = 0:K-1; x = pixel_pitch*( x-x( end/2 ) );
y = 0:K-1; y = pixel_pitch*( y-y( end/2 ) );
[ X,Y ] = meshgrid( x,y );
lens_function = exp( 1i*k0*( X.^2+Y.^2 )/( 2*f1 ) );
psi = psi.*lens_function;

%% propagate object to lens 2

Psi = fftshift( fft2( fftshift( psi ) ) );

switch diffraction
    case 'fresnel'
        H1 = exp( -1i*k0*( f1+f2 ) )*exp( 1i*( f1+f2 )*( Kx.^2+Ky.^2 )/( 2*k0 ) );
        psi = fftshift( ifft2( fftshift( Psi.*H1 ) ) );
    case 'angular spectrum'
        H1 = exp( -1i*k0*( f1+f2 )*sqrt( 1-( Kx/k0 ).^2-( Ky/k0 ).^2 ) );
        psi = fftshift( ifft2( fftshift( Psi.*H1 ) ) );
end

% apply phase transfer function of lens 2
lens_function = exp( 1i*k0*( X.^2+Y.^2 )/( 2*f2 ) );
psi = psi.*lens_function;

%% propage object to CCD

Psi = fftshift( fft2( fftshift( psi ) ) );

switch diffraction
    case 'fresnel'
        H1 = exp( -1i*k0*f2 )*exp( 1i*f2*( Kx.^2+Ky.^2 )/( 2*k0 ) );
        psi = fftshift( ifft2( fftshift( Psi.*H1 ) ) );
    case 'angular spectrum'
        H1 = exp( -1i*k0*f2*sqrt( 1-( Kx/k0 ).^2-( Ky/k0 ).^2 ) );
        psi = fftshift( ifft2( fftshift( Psi.*H1 ) ) );
end


% reduce diffraction plane
psi_crop = psi( K/2-M/2+1:K/2+M/2,K/2-M/2+1:K/2+M/2 ); 

% zero pad in frequency domain
psi = zeros( us_rate*M );
psi( us_rate*M/2-M/2+1:us_rate*M/2+M/2,us_rate*M/2-M/2+1:us_rate*M/2+M/2 )=...
    fftshift( ifft2( fftshift( psi_crop ) ) );
Psi = fftshift( fft2( fftshift( psi ) ) );
Psi = fliplr( flipud( Psi ) );

figure
subplot( 1,2,1 )
imagesc( 1e3*x,1e3*y,abs( Psi ).^2 )
axis xy;xlabel( 'x (mm)' ); ylabel( 'y (mm)' );title( 'Intensity' )
subplot( 1,2,2 )
imagesc( 1e3*x,1e3*y,180/pi*atan2( imag( Psi ),real( Psi ) ) )
axis xy;xlabel( 'x (mm)' ); ylabel( 'y (mm)' );title( 'Phase' )

%% interference at the hologram plane

% generate reference wave
A = ( min( abs( Psi( : ) ) )+max( abs( Psi( : ) ) ) )/2;
theta = 0.75; % reference wave angle ( deg )
m = 1:K; m = m-m( end/2 );
n = 1:K; n = n-n( end/2 );
dx = pixel_pitch/( us_rate );
dy = pixel_pitch/( us_rate );
[ X,Y ] = meshgrid( m*dx,n*dy );
% off-axis plane wave
Psir = A*exp( 1i*k0*sind( theta )*( X+Y ) );

hologram = abs( Psi + Psir ).^2;

Kr = sqrt( Kx.^2+Ky.^2 );

Hologram = ifft2( fftshift( hologram ) );

figure
subplot( 1,2,1 )
imagesc( 1e-3*kx,1e-3*ky,20*log10( fftshift( abs( Hologram ) ) ) )
xlabel( 'Spatial Frequency k_{x} ( mm^{-1} )' )
ylabel( 'Spatial Frequency k_{y} ( mm^{-1} )' )
title( 'Hologram Spectral Power ( dB )' )
caxis( [ -100 0 ] );colorbar;colormap( jet );xlim( [ -50 50 ] );ylim( [ -50 50 ] );
axis xy;set(gca,'FontSize',20,'FontWeight','b')
subplot( 1,2,2 )
semilogy( 1e-3*kx,diag( fftshift( abs( Hologram ) ) ),'b','LineWidth',1 )
xlabel( 'Spatial Frequency k ( mm^{-1} )' )
ylabel( 'Hologram Spectral Power' )
grid;set(gca,'FontSize',14,'FontWeight','b')
ylim( [ 1e-6 1e1 ] );xlim( [ -50 50 ] );

%% reconstruction
hologram = hologram - mean( hologram( : ) ); % remove DC term

psi_est = hologram.*Psir;
Psi_est = fftshift( fft2( fftshift( psi_est ) ) );

k_shift = k0*sind( theta )/( 2*pi*us_rate ); % m^-1
[ kx_shift,indkx ] = min( abs( kx+k_shift ) );
[ ky_shift,indky ] = min( abs( ky+k_shift ) );

kr = 40e3;

switch lowpass_filter
    case 'ideal'
        lpf = ( Kx.^2 + Ky.^2 ) <= kr^2;      
    case 'butterworth'
        D0 = kr; % passband width
        n = 8; % filter order
        D = sqrt( Kx.^2+Ky.^2 );
        D( D==0 ) = 1; % prevent divide by zero
        lpf = 1-1./( 1+( D0./D ).^( 2*n ) ); % frequency response
        lpf( length( kx )/2+1,length( ky )/2+1 ) = 1;      
    case 'gaussian'
        nstd = 3;
        [ lpf,~,~ ] = bivariate_Gaussian( [ K K ],[ 0 0 ],[ ( kr/( kx( 2 )-kx( 1 ) )/nstd )^2 ( kr/( kx( 2 )-kx( 1 ) )/nstd )^2 ],0 );
        lpf = lpf/max( lpf( : ) );
end

Psi_est = Psi_est.*lpf;
psi_est = fftshift( ifft2( fftshift( Psi_est ) ) );

%psi_est = psi_est.*conj( psi_est_mag );

psi0_est = conj( -f1/f2*exp( -2i*k0*( f1+f2 ) ) )*psi_est;
Psi0_est = fftshift( fft2( fftshift( psi0_est ) ) );

% downsample to original image size
ds_factor = size( psi0_est,1 )/size( im,1 );

% reconstruct intensity and phase
reconstructed_phase = atan2( imag( psi0_est ),real( psi0_est ) );

%% plots

xi = 0:M-1; xi = pixel_pitch*( xi-xi( end/2 ) );
yi = 0:N-1; yi = pixel_pitch*( yi-yi( end/2 ) );

xi = xi/3.5759; yi = yi/3.5759;
x = x/3.5759*5.785/4.32; y = y/3.5759*5.785/4.32;

figure
imagesc( 1e3*x,1e3*y,1e9*lambda*reconstructed_phase/( 2*pi )/(n1-n0 )  )
axis xy;xlabel( 'x (mm)' ); ylabel( 'y (mm)' );zlabel( 'z (nm)' );title( 'Reconstructed' )
caxis( [ 0 150 ] );cb=colorbar;set(gca,'FontSize',14,'FontWeight','b');shading flat
set(get(cb,'title'),'string','nm');
ylim( [ -5.85 5.85 ] );xlim( [ -5.85 5.85 ] )
figure
plot( 1e3*xi*us_rate^2,1e9*lambda*image_phase( N/2,: )/( 2*pi )/(n1-n0 ),'k','LineWidth',2 );hold on
plot( 1e3*x,1e9*lambda*reconstructed_phase( K/2,: )/( 2*pi )/(n1-n0 ),'r','LineWidth',2 )
xlabel( 'x (mm)' ); ylabel( 'z (nm)' )
legend( 'True','Reconstructed')
set(gca,'FontSize',14,'FontWeight','b')
ylim( [ 0 152 ] );xlim( [ -5.85 5.85 ] );grid
colormap( jet )


