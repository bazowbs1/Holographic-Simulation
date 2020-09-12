function [ Psi_pca ] = aberration_compensation_pca( Psi,ind )

Psi_pca = Psi;
ROI = Psi( ind( 1,: ),ind( 2,: ) ); %Region Of Interest in frequency domain size M x N
[ N,M ] = size( ROI );
%IFFT the cropped spectrum to get a subsampled hologram
IFFT_ROI = fftshift( ifft2( fftshift( ROI ) ) );
%get the exponential term
ConjPhase = exp( -1i*angle( IFFT_ROI ) );
%singular value decomposition
[ U,S,V ] = svd( ConjPhase ); 

SS = zeros( N,M );
num = 1;%size( S,1 ); %number of principal components will be taken
for i=1:num %take first 'num' values
    SS( i,i ) = S( i,i );
end

%least-squares fitting
Unwrap_U = unwrap( angle( U( :,1:2 ) ) );
SF_U1 = polyfit( 1:N,Unwrap_U( :,1 )',2 ); %second degree
SF_U2 = polyfit( 1:N,Unwrap_U( :,2 )',2 ); %second degree
EstimatedSF_U1 = polyval( SF_U1,1:N );
EstimatedSF_U2 = polyval( SF_U2,1:N );
New_U1 = exp( 1i*( EstimatedSF_U1' ) );
New_U2 = exp( 1i*(EstimatedSF_U2' ) );
U = U*0;
U( :,1:2 ) = [ New_U1 New_U2 ];

Unwrap_V = unwrap( angle( V( :,1:2 ) ) );
SF_V1 = polyfit( 1:N,Unwrap_V( :,1 ).',2 ); %second degree
SF_V2 = polyfit( 1:N,Unwrap_V( :,2 ).',2 ); %second degree
EstimatedSF_V1 = polyval( SF_V1,1:N );
EstimatedSF_V2 = polyval( SF_V2,1:N );
New_V1 = exp( 1i*( EstimatedSF_V1.' ) );
New_V2 = exp( 1i*( EstimatedSF_V2.' ) );
V = V*0;
V( :,1:2 ) = [ New_V1 New_V2 ];

%get the aberration term
Z = U*SS*V';

%FFT and replace the corresponding original region of the spectrum
Psi_pca( ind( 1,: ),ind( 2,: ) ) = fftshift( fft2( fftshift( IFFT_ROI.*Z ) ) );

end