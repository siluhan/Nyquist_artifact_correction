% returnimage = qft(input)
function returnimage = qft(input)

returnimage = fftshift(fftshift(fft(fft(fftshift(fftshift(squeeze(input),1),2),[],1),[],2),1),2);

