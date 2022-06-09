
%this is for removing the scanner lines
%please refer to the independent code sca_lin_trial_1.m
%the threshold value and the 
function filt_im = scan_line_rem(im)
    
    log_thresh = 16 ;
    %img = Tiff(file_name);
    %img = read(img);
    %imgdisp = imcrop(img,[2530 1044 159 135]);
    %subplot(2,3,1);
    %imshow(imgdisp);

    % Compute the 2D fft.
    frequencyImage = fftshift(fft2(im));
    % Take log magnitude so we can see it better in the display.
    amplitudeImage = log(abs(frequencyImage));
    minValue = min(min(amplitudeImage));
    maxValue = max(max(amplitudeImage));
    %subplot(2,3,2);
    %imshow(amplitudeImage, []);
    %title(caption, 'FontSize', fontSize);
    %axis on;

    amplitudeThreshold = log_thresh;
    brightSpikes = amplitudeImage > amplitudeThreshold; % Binary image.
    %subplot(2, 3, 3);
    %imshow(brightSpikes);
    %axis on;
    
    %please check the filtering values for every set of data
    %brightSpikes(1715:1785, :) = 0;
    %brightSpikes(1700:1800, :) = 0;
    brightSpikes(1650:1850, :) = 0;
    brightSpikes(:,1:1900) = 0;
    brightSpikes(:,2100:end) = 0;
    %subplot(2, 3, 4);
    %imshow(brightSpikes);
    %axis on;

    % Filter/mask the spectrum.
    frequencyImage(brightSpikes) = 0;
    %frequencyImage = frequencyImage - brightSpikes;
    % Take log magnitude so we can see it better in the display.
    %amplitudeImage2 = log(abs(frequencyImage));
    %subplot(2, 3, 5);
    %imshow(log(abs(frequencyImage)),[]);
    %imshow(amplitudeImage2);
    %bleh
    %axis on;

    filteredImage = ifft2(fftshift(frequencyImage));
    filteredImage = abs(filteredImage);
    %imgdisp2 = imcrop(filteredImage,[2530 1044 159 135]);
    %subplot(2,3,6);
    %imagesc(filteredImage);
    %imshow(filteredImage,[]);
    %bleh
    %imshow(imgdisp2,[])
    filt_im = filteredImage;
end
