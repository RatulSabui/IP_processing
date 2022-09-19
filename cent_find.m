%% this function finds the holes on the angular distribution plate
function [cent,select_cent] = cent_find(img)

    %tagNames = img_fil.getTagNames();
    %imshow(B);
    %T = graythresh(B)
    T = 0.0005;
    bin = imbinarize(img,T);
    %figure('name','binary');
    imshow(bin)
    %bleh
    bin = bwareaopen(bin,1000);
    com = imcomplement(bin);
    rem = bwareaopen(com,1000);
    rem2 = rem - bwareaopen(rem,100000);
    %imshow(rem2);
    %rem3 = imrotate(rem2,-50);
    %imshow(rem3);
    [B,L] = bwboundaries(rem2,'noholes');
    %bwim = rem2;
    
    %%this part is just for plotting the circles and the boundaries
    imshow(label2rgb(L,@jet,[.5 .5 .5]))
    h = gca;
    h.Visible = 'On';
    hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2),boundary(:,1),'w','LineWidth',2)
    end

    stats = regionprops(L,'Area','Centroid');
    
    struct2table(stats)
    cent = stats;
%    threshold = 0.94;

    
    area_array = [];
    % loop over the boundaries
    % selecting the 2 biggest regions
    for k = 1:length(B)

        % obtain (X,Y) boundary coordinates corresponding to label 'k'
        %boundary = B{k};

        % compute a simple estimate of the object's perimeter
        %delta_sq = diff(boundary).^2;
        %perimeter = sum(sqrt(sum(delta_sq,2)));

        % obtain the area calculation corresponding to label 'k'
        %area = stats(k).Area;

        % compute the roundness metric
        %metric = 4*pi*area/perimeter^2;

        % display the results
%         metric_string = sprintf('%2.2f',metric);

        % mark objects above the threshold with a black circle
%         if metric > threshold
            %centroid = stats(k).Centroid;
            %plot(centroid(1),centroid(2),'ko');
%         end
          area_temp = stats(k).Area;
          area_array(end+1) = area_temp;
%         text(boundary(1,2)-35,boundary(1,1)+13,metric_string,'Color','y',...
%             'FontSize',14,'FontWeight','bold')

    end
    
    %the next part selects out the biggest 
    [max1, ind1] = max(area_array);
    area_array(ind1)      = -Inf;
    [max2, ind2] = max(area_array);
    
    max_area(1) = stats(ind1);
    max_area(2) = stats(ind2);
    %this var cocntains  the data about the two largest holes detected
    select_cent = max_area;
    struct2table(max_area)
    %[centers,radii] = imfindcircles(rgb,[20 25],'ObjectPolarity','dark')
    
    for k = 1:length(max_area)
        centroid = max_area(k).Centroid;
        plot(centroid(1),centroid(2),'ko');
    end
    
    
    %pixel_vals = impixel
    %imhist(B);
    %imshow(his);
    %bin = imbinarize(B);
    %imshow(bin);
    %bin = imquantize(B,50);
    %imshow(bin);
end