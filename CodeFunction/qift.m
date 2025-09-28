% returnimage = qift(input)
function returnimage = qift(input)

returnimage = ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(squeeze(input),1),2),[],1),[],2),1),2);

