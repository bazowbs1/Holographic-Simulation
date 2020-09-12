%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation of holographic recording and reconstruction using an off-axis
% digital holographic microscopy configuration. The optical field at the
% recording plane is computed using a generalized approach based on matrix
% optics and Fresnel diffraction. The optical field is digitally propagated
% to the reconstruction plane using Fresnel diffraction, the angular
% spectrum method, or the Rayleigh-Sommerfeld formulae. Aberrations are
% removed using principal componets, higher order aberrations are removed
% using Zernike polynomials
%
% written by Brad Bazow ( September, 2020 )
% The Johns Hopkins University Applied Physics Labratory ( JHU/APL )
% The Catholic University of America Electrical and Computer Science
% Department ( CUA/EECS )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all; 

diffraction = 'fresnel'; % 'rayleigh-sommerfeld', 'angular spectrum', 'fresnel'
lowpass_filter = 'butterworth'; % 'ideal', 'butterworth', 'gaussian'

%% parameters
lambda = 632.8e-9; % wavelength ( m )
k0 = 2*pi/lambda; % wavenumber ( rad/m )
pixel_pitch = 7.1250e-06;% pixel spacing ( m/pixel )
us_rate = 4; % image upsample rate
f = 16.5e-2; % focal length of ideal lens ( m )
d0 = 20e-2; % distance from object plane to lens ( m )
di = 100e-2; % distance from lens to image plane ( m )
d = 20e-2; % distance from FPA to image plane ( m )
d1 = di-d; % distance from lens to FPA ( m )
magnification = di/d0; % magnification
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

v = 1:K; v = v-v(end/2)-1; 
u = 1:K; u = u-u(end/2)-1;

%spatial frequencies paper: Simultaneous amplitude-contrast ...
[ Ky,Kx ] = meshgrid( v/( K*pixel_pitch),u/( K*pixel_pitch ) );
%reconstructed image resolution
[ YY,XX ] = meshgrid( v*( lambda*d )/( K*pixel_pitch*magnification ),u*( lambda*d )/( K*pixel_pitch*magnification ) );

% ----- ABCD matrix (A~=0) -----
A = 1-( di-d )/f; B = ( 1-( di-d )/f )*d0+di-d; C = -1/f; D = 1-d0/f; z = d0+di-d;

factor = exp( -1i*k0*C*( XX.^2+YY.^2 )/( 2*A ) );
U_Fourier = fftshift( fft2( fftshift( CustomizedZoom( psi0,abs( A ) )/A ) ) );
kernel = exp( -1i*k0*z )*exp( 1i*B*A*( Kx.^2+Ky.^2 )/( 2*k0 ) );
psi = factor.*fftshift( ifft2( fftshift( U_Fourier.*kernel ) ) );

% reduce diffraction plane
psi_crop = psi( K/2-M/2+1:K/2+M/2,K/2-M/2+1:K/2+M/2 ); 

% zero pad in frequency domain
psi = zeros( us_rate*M );
psi( us_rate*M/2-M/2+1:us_rate*M/2+M/2,us_rate*M/2-M/2+1:us_rate*M/2+M/2 )=...
    fftshift( ifft2( fftshift( psi_crop ) ) );
Psi = fftshift( fft2( fftshift( psi ) ) );
 
x = 0:K-1; x = pixel_pitch*( x-x( end/2 ) );
y = 0:K-1; y = pixel_pitch*( y-y( end/2 ) );

figure
subplot( 1,2,1 )
imagesc( 1e3*x,1e3*y,abs( Psi ).^2 )
axis xy;xlabel( 'x (mm)' ); ylabel( 'y (mm)' );title( 'Intensity' )
subplot( 1,2,2 )
imagesc( 1e3*x,1e3*y,180/pi*atan2( imag( Psi ),real( Psi ) ) )
axis xy;xlabel( 'x (mm)' ); ylabel( 'y (mm)' );title( 'Phase' )

%% Generate Hologram

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

kx = ( ( 0:K-1 )-K/2 )/K/pixel_pitch;
ky = ( ( 0:K-1 )-K/2 )/K/pixel_pitch;
[ Kx,Ky ] = meshgrid( kx,ky );
Kr = sqrt( Kx.^2+Ky.^2 );

Hologram = ifft2( fftshift( hologram ) );

figure
subplot( 1,2,1 )
imagesc( 1e-3*kx,1e-3*ky,20*log10( fftshift( abs( Hologram ) ) ) )
xlabel( 'Spatial Frequency k_{x} ( mm^{-1} )' )
ylabel( 'Spatial Frequency k_{y} ( mm^{-1} )' )
title( 'Hologram Spectral Power ( dB )' )
caxis( [ -100 0 ] );colorbar;colormap( jet ); xlim( [ -10 10 ] );ylim( [ -10 10 ] );
axis xy;set(gca,'FontSize',20,'FontWeight','b')
subplot( 1,2,2 )
semilogy( 1e-3*kx,diag( fftshift( abs( Hologram ) ) ),'b','LineWidth',1 )
xlabel( 'Spatial Frequency k ( mm^{-1} )' )
ylabel( 'Hologram Spectral Power' )
grid;set(gca,'FontSize',14,'FontWeight','b')
ylim( [ 1e-6 1e1 ] );

%% reconstruction

hologram = hologram - mean( hologram( : ) ); % remove DC term

psi_est = hologram.*Psir;
Psi_est = fftshift( fft2( fftshift( psi_est ) ) );

%% propagate pbject to reconstruction plane
x = 0:K-1; x = pixel_pitch*( x-x( end/2 ) );
y = 0:K-1; y = pixel_pitch*( y-y( end/2 ) );
[ X,Y ] = meshgrid( x,y );
m = ( 0:K-1 )-K/2;
n = ( 0:K-1 )-K/2;
dkx = 1/( us_rate*M*pixel_pitch );
dky = 1/( us_rate*N*pixel_pitch );
[ Kx,Ky ] = meshgrid( 2*pi*m*dkx,2*pi*n*dky );

switch diffraction
    case 'fresnel'
        H = exp( -1i*k0*d )*exp( 1i*d*( Kx.^2+Ky.^2 )/( 2*k0 ) );
        psi_est = fftshift( ifft2( fftshift( Psi_est.*H ) ) );
    case 'rayleigh-sommerfeld'
        h = 1i*k0/( 2*pi )*exp( -1i*k0*sqrt( d^2+X.^2+Y.^2 ) )./( sqrt( d^2+X.^2+Y.^2 ) )*d./( sqrt( d^2+X.^2+Y.^2 ) ).*( 1+ones( size( X ) )./( 1i*k0*sqrt( d^2+X.^2+Y.^2 ) ) );
        H = fftshift( fft2( fftshift( h ),K,K ) );
        psi = fftshift( ifft2( fftshift( Psi.*H ) ) );        
        psi_temp = zeros( size( psi ) );
        psi_temp( K/2-M/2:K/2+M/2,K/2-M/2:K/2+M/2 ) = psi( K/2-M/2:K/2+M/2,K/2-M/2:K/2+M/2 );
        psi = psi_temp;
    case 'angular spectrum'
        H = exp( -1i*k0*d*sqrt( 1-( Kx/k0 ).^2-( Ky/k0 ).^2 ) );
        psi_est = fftshift( ifft2( fftshift( Psi_est.*H ) ) );
end

Psi_est = fftshift( fft2( fftshift( psi_est ) ) );

k_shift = k0*sind( theta )/( 2*pi*us_rate ); % m^-1
[ kx_shift,indkx ] = min( abs( kx+k_shift ) );
[ ky_shift,indky ] = min( abs( ky+k_shift ) );

[ Kx,Ky ] = meshgrid( kx,ky );
kr = 5e3;
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


%% PCA aberration compensation
ind( 1,: ) = K/2-N/2+1:K/2+N/2;
ind( 2,: ) = K/2-N/2+1:K/2+N/2;
[ Psi_est_pca ] = aberration_compensation_pca( Psi_est,ind );
psi_est_pca = fftshift( ifft2( fftshift( Psi_est_pca ) ) );
reconstructed_phase_pca = atan2( imag( psi_est_pca ),real( psi_est_pca ) );

%% Zernike aberration compensation
ds_factor = size( psi_est,1 )/size( im,1 );

reconstructed_phase = atan2( imag( psi_est ),real( psi_est ) );
reconstructed_phase = downsample(downsample( reconstructed_phase,ds_factor ).',ds_factor ).';

x = 1:K; x = x-x(end/2)-1; %M x N
y = 1:K; y = y-y(end/2)-1;
[ X,Y ] = meshgrid( x*pixel_pitch,y*pixel_pitch );

reconstructed_phase_pca = downsample( downsample( reconstructed_phase_pca,ds_factor ).',ds_factor ).';

num = 4; % number of rectangles used to estimate phase aberration
[ phase_zernike ] = aberration_compensation_zernike( reconstructed_phase_pca,X,Y,num );
reconstructed_phase_zernike = reconstructed_phase_pca - phase_zernike;

%% plots

x = 0:M-1; x = pixel_pitch*( x-x( end/2 ) );
y = 0:N-1; y = pixel_pitch*( y-y( end/2 ) );
[ X,Y ] = meshgrid( x,y );

xi = x*1.2996; yi = y*1.2996;

figure
subplot( 1,3,1 )
imagesc( 1e3*xi,1e3*yi,1e9*lambda*reconstructed_phase/( 2*pi )/(n1-n0 ) );axis xy
subplot( 1,3,2 )
imagesc( 1e3*xi,1e3*yi,1e9*lambda*reconstructed_phase_pca/( 2*pi )/(n1-n0 ) );axis xy
subplot( 1,3,3 )
imagesc( 1e3*xi,1e3*yi,1e9*lambda*reconstructed_phase_zernike/( 2*pi )/(n1-n0 ) );axis xy

figure
subplot( 1,2,1 )
imagesc( 1e3*xi,1e3*yi,1e9*lambda*reconstructed_phase_zernike/( 2*pi )/(n1-n0 ) )
axis xy;xlabel( 'x (mm)' ); ylabel( 'y (mm)' );zlabel( 'z (nm)' );title( 'Reconstructed' )
caxis( [ 0 150 ] );cb=colorbar;set(gca,'FontSize',14,'FontWeight','b');shading flat
set(get(cb,'title'),'string','nm');
subplot( 1,2,2 )
plot( 1e3*x*magnification,1e9*lambda*image_phase( N/2,: )/( 2*pi )/(n1-n0 ),'k','LineWidth',2 );hold on
plot( 1e3*xi,1e9*lambda*reconstructed_phase_zernike( N/2,: )/( 2*pi )/(n1-n0 ),'r','LineWidth',2 )
xlabel( 'x (mm)' ); ylabel( 'z (nm)' )
legend( 'True','Reconstructed')
set(gca,'FontSize',18,'FontWeight','b')
ylim( [ 0 152 ] ); xlim( [ -6 6 ] )
