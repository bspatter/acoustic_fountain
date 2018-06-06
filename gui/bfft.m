function [ freqdata ] = bfft( timedata )
    [rows,columns]=size(timedata);
    
    if columns>rows %I want each column to be the complete fourier information
        timedata=transpose(timedata);
    end
    
    L=round(max(size(timedata))/2)*2; %Account for possibility that L is not a power of 2
    N=min(size(timedata));
    freqdata=zeros(L/2+1,N);
    
    for i=1:N
        temptimedata=timedata(:,i); %Assign the ith column
        temp=conj(fft(temptimedata)); %FFT Convention enforced
        freqdata(:,i)=temp(1:L/2+1); %Grab only the information before the folding frequency
    end
    
    if columns>rows %switch it back if necessary
        freqdata=transpose(freqdata);
    end

end

