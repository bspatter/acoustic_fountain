function [ timedata ] = bifft( freqdata )
    [rows,columns]=size(freqdata);
    
    if columns>rows %I want each column to be the complete fourier information
        freqdata=transpose(freqdata);
    end
    
    L=2*max(size(freqdata))-2;
    N=min(size(freqdata));
    timedata=zeros(L,N);
    
    for i=1:N
        tempfreqdata=freqdata(:,i); %Assign the ith column
        
        tempfreqdata(1)=real(tempfreqdata(1)); %Force the DC bin to be real, otherwise the time domain has a constant imaginary component. As below, a sign change here doesn't change the realness. Not sure what's "most" correct lol...

        L=2*length(tempfreqdata)-2; %Length of time output, such that Length(freq)=L/2+1
        preIFFT=zeros(L,1); %Prepare a full-length vector for IFFT
        preIFFT(1:L/2+1)=tempfreqdata; %Assign the freq data to first half
        temp=conj(flipud(tempfreqdata)); %Flip freq data, conjugate and store
        preIFFT(L/2+1:L)=temp(1:L/2); %Assign the flipped data, skipping the first data point
        preIFFT(L/2+1)=real(preIFFT(L/2+1)); %Force the folding frequency to be real (abs or real part?) - I'm choosing the real part, since having a negative number here doesn't destroy the real-ness of the time data
        timedata(:,i)=ifft(conj(preIFFT)); %Fixed for FFT Convention
    end
    
    if columns>rows %switch it back if necessary
        timedata=transpose(timedata);
    end
        
    
    %Psuedocode for the preIFFT setup:
    %Step 1: 0 0 0 0 0 0 // Create Empty Vector of Zeroes
    %Step 2: 1 2 3 4 0 0 // Assign first half of Freq Data from 1 to L/2+1
    %Step 3: store 4 3 2 1 // FlipLR, Conjugate, and store
    %Step 4: 1 2 3 4 3 2 // Assign this stored info to L/2+1 to L
    %Step 5: 1 2 3 4 3 2 // Force the "4" to become real
    %Step 6: IFFT
end

