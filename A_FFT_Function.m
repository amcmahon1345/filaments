function resampledFP1 = A_FFT_Function(X,Y)
% Fast fourier transform given X and Y data. f is the x scale for the
% fourier transform and scaledP1 is the y variable for the fourier
% transform.

%Plot the input function to be Fourier transformed
%{
figure
plot(X,Y)
xlabel('Distance along filament/nm')
ylabel('Intensity')
%}

L = length(X); 
fftY = fft(Y);
P2 = abs(fftY/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = (0:(L/2))/L; %Needs frequency factor

%This gives the option for the first datapoint not to count in the scaling
%if you are interested in higher frequency modes.
scale = max(P1(2:end));
scaledP1 = P1./scale;

resampledFP1 = interp1(f, scaledP1, [0:0.0005:0.5], 'cubic');

%{
figure
plot([0:0.0005:0.5],resampledNormalisedFP1)
xlabel('Frequency/Length')
ylabel('Normalised intensity')
%}