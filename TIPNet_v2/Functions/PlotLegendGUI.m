function [X] = PlotLegendGUI(fh,Weights,colorlabel)
% arrow and color legend generation for circular network plots

% plot nothing except for colorbar with labels and range of arrow sizes

axes(fh)
caxis([0 max(Weights)])


nW = length(Weights);
b = jet(nW);
colormap(b);
for i=1:nW
   labels{i}=num2str(Weights(i)); 
end
hc = lcolorbar(labels);


set(hc,'position',[0.85    0.525    0.02    0.4])
set(fh,'Visible','Off')

texth = text(.45, .75,sprintf('%s',colorlabel),'FontSize',10,'Parent',fh);
%texth = text(.1, .65,sprintf('bits/bit'),'FontSize',10,'Parent',fh);


% %draw arrows of different sizes
%     ln = line([.1, .25],[.3, .3],'Color',...
%         'b','LineWidth',10,'Parent',fh);
%     texth = text(.28, .3,sprintf('1'),'FontSize',12,'Parent',fh);
%  
%     ln = line([.1, .25],[.4, .4],'Color',...
%         'b','LineWidth',5,'Parent',fh);
% 
%     texth = text(.28, .4,sprintf('0.5'),'FontSize',12,'Parent',fh);
%     
%     ln = line([.1, .25],[.5, .5],'Color',...
%         'b','LineWidth',1,'Parent',fh);
% 
%     texth = text(.28, .5,sprintf('0.1'),'FontSize',12,'Parent',fh);
%     
%     ln = line([.1, .25],[.6, .6],'Color',...
%         [.5 .5 .5],'LineWidth',.5,'LineStyle',':','Parent',fh);
% 
%     texth = text(.28, .6,sprintf('< 0.1'),'FontSize',12,'Parent',fh);
%     
    set(fh,'color','none')

end



