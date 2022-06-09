%% Calculates the PSL value of each pixel in the image
function img_psl = psl_calc(img,res,sensit,lat,ran)
    %imshow(img,[])
    %imhist(img)
    %set(gca,'YScale','log')
    data = double(img);
    Res = res;
    L = lat;
    S = sensit;
    G = (2^ran) - 1;
    h = 3.5;
    row=size(data,1);
    col=size(data,2);
    %PSLtotal=0;
    psl = (h*(10^(L/2))*((data./G).^2)*(Res/100)^2);
    %PSLnew = psl;
    PSLnew = psl/(0.65498+ 0.18075*exp(-13/35)+0.18225*exp(-13/264));
%     for i=1:row
%         for j =1:col
%             %signal(i,j)=data(i,j);
%             %signal(i,j)=data(i,j)-noise(i,j);
%             %gray value to PSL conversion FLA 7000
%             %PSL(i,j)=double(((Res/100)^2)*(10^(L*((signal(i,j)/G)-0.5))));
%             %gray value to PSL conversion for GE
%             PSL(i,j)=double(h*(10^(L/2))*((data(i,j)/G)^2)*(Res/100)^2);
%             PSLnew(i,j)= PSL(i,j)/(0.65498+ 0.18075*exp(-13/35)+0.18225*exp(-13/264));
%             if data(i,j)>60000
%                PSL(i,j)
%                PSLnew(i,j) 
%                %fprintf("bleh")
%             end
%             %PSLtotal=PSLtotal+PSLnew(i,j);
%         end
%     end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    img_psl = PSLnew;
    %imshow(img_psl,[])
    %imhist(img_psl)
    %set(gca,'YScale','log')
    %fprintf("bleh")
end