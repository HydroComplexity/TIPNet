function [X] = PlotFunGUI(fh, Links,Weights, WeightRange, nodesize, varnames_s,titletext,A,uplim,bothways)
%Circular network plots

% fh: figure handle
% pn: position
% Links: matrix of link strengths (represented by line weight)
% Weights: time lags or other weight (will be color of line)
% WeightRange: range of weights [min max number_of_weight_values]
% H_y_s: node sizes
% varnames_s: names of nodes for labels
% titletext: text for title
% A = 1 for arrow
% uplim = upper limit for link plotting if don't want to plot all links
% bothways = 1 to draw bi-directional links, 0 to only draw larger link

%to match legend: max line width should =10
axes(fh)

%set(gcf, 'color', 'none','inverthardcopy', 'off');
nvars = size(Links,2);

nodesize(isnan(nodesize))=0;

w = 4;

set(gca,'Visible','off')
hold on
axis equal
%polar plotting - big circle
theta = linspace(0,2*pi,nvars+1);
theta = theta(1:end-1);
[x, y] = pol2cart(theta,.7);
[xlab, ylab] = pol2cart(theta,.9);
xlab = (abs(xlab) +0.05).*sign(xlab)-.2; %shift label a bit

Weight_Ints = Weights;

minW = WeightRange(1);
maxW = WeightRange(2);
nW = WeightRange(3);

wt_coords =linspace(minW,maxW,nW+1);
wts = 1:(nW+1);

for w=1:nW
    Weight_Ints(Weights>=wt_coords(w) & Weights<=wt_coords(w+1))=wts(w);
end

Weight_Ints(Weights>wt_coords(w+1))=wts(w);

colorvect = jet(nW);
Weights = Weight_Ints;

%scale nodesize ('markersize' is pixels in diameter of circle, 72 pixels in
%an inch)
nodesize = sqrt(200*4.*nodesize./pi());
Rnode = nodesize./2;
Rnode = Rnode./72; %back to pixels

%show top links in color in front, lesser links in gray in back
BigLinks=Links;
dd = reshape(Links,1,size(Links,1)*size(Links,2));
if uplim ==0
 uplim = prctile(dd,99);  
end

BigLinks(Links<uplim)=0;

for i = 1:nvars
    for j = 1:nvars
        
        if i~=j && BigLinks(i,j)>0 %want line from node i to node j
            
            if bothways ==1
                wid = Links(i,j)*10;
                
                theta = atan((y(j)-y(i))/(x(j)-x(i))); %angle of incoming arrow
                xdist = abs(Rnode(j).*cos(theta));
                ydist = abs(Rnode(j).*sin(theta));
                if y(j)>y(i)
                    y2 = y(j)-ydist;
                else
                    y2 = y(j)+ydist;
                end
                
                if x(j)>x(i)
                    x2 = x(j) - xdist;
                else
                    x2 = x(j) + xdist;
                end
                    
                    
                ln = line([x(i), x2],[y(i), y2],'Color',...
                    colorvect(Weights(i,j),:),'LineWidth',wid,'Parent',fh);
                
                if A==1
                   arrow(ln,...
                       'FaceColor','none',...
                       'TipAngle',25);
                end
            else %only draw line if larger than opposite direction
                if Links(i,j)>Links(j,i)
                    wid = Links(i,j)*10;
                    
                    theta = atan((y(j)-y(i))/(x(j)-x(i))); %angle of incoming arrow
                    xdist = abs(Rnode(j).*cos(theta));
                    ydist = abs(Rnode(j).*sin(theta));
                    if y(j)>y(i)
                        y2 = y(j)-ydist;
                    else
                        y2 = y(j)+ydist;
                    end
                    
                    if x(j)>x(i)
                        x2 = x(j) - xdist;
                    else
                        x2 = x(j) + xdist;
                    end
                    
                    ln = line([x(i), x2],[y(i), y2],'Color',...
                        colorvect(Weights(i,j),:),'LineWidth',wid,'Parent',fh);
                    
                    if A==1
                        arrow(ln,...
                       'FaceColor','none',...
                       'TipAngle',25);
                    end
                end                           
            end               
        end
    end
end

%put smaller links in front of larger links
if bothways ==1
for i = 1:nvars   
    for j = 1:nvars
        
    if i~=j && BigLinks(i,j)>0 && BigLinks(i,j) < BigLinks(j,i) %want line from node i to node j

    wid = Links(i,j)*10;
    
    theta = atan((y(j)-y(i))/(x(j)-x(i))); %angle of incoming arrow
    xdist = abs(Rnode(j).*cos(theta));
    ydist = abs(Rnode(j).*sin(theta));
    if y(j)>y(i) %subtract y
    y2 = y(j)-ydist;
    else
    y2 = y(j)+ydist;
    end
    
    if x(j)>x(i)
    x2 = x(j) - xdist;
    else
    x2 = x(j) + xdist;    
    end

    
    ln = line([x(i), x2],[y(i), y2],'Color',...
        colorvect(Weights(i,j),:),'LineWidth',wid,'Parent',fh);
    
    if A==1
                   arrow(ln,...
                       'FaceColor','none',...
                       'TipAngle',25);
    end
    end    
    end
end
end


%plot nodes (plot color the same as arrow values)
for j = 1:nvars  
    if ~isnan(Weights(j,j))
    m(j) = plot(x(j),y(j),'o','markersize',nodesize(j),'LineWidth',1.5,...
        'Color','k','MarkerFaceColor',colorvect(Weights(j,j),:),'Parent',fh);
    else
    m(j) = plot(x(j),y(j),'o','markersize',nodesize(j),'LineWidth',1.5,...
        'Color','k','MarkerFaceColor','k','Parent',fh);
    end
    text(xlab(j),ylab(j),sprintf('%s',varnames_s{j}),'FontSize',10,'Parent',fh,'FontWeight','bold');
end

set(gca, 'color', 'none');

X=2;

end

