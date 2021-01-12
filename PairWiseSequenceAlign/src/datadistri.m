clear ;close all; clc
Y =[29860, 29859, 29896, 29782, 29782, 29853, 29832, 29846, 29841, 29853, 29882, 29882, 29882, 29882, 29882, 29786, 29788, 29805, 29791, 29812, 29409, 29409, 29409, 29409, 29894, 29864, 29903, 29896, 29903];
figure
barh(Y,'LineWidth',1.5,'BarWidth',0.6,'FaceColor',[128 128 105]/255,'EdgeColor',[0 0 0]/255)
%set(gca,'YLim',[0 1.02]);%X轴的数bai据显示du范围
%set(gca,'ytick',[0:0.1:1.02]);%设置要显示坐标刻度
% legend({'--Accuracy'},'Location','best','FontName','Times New Roman','FontWeight','Bold','FontSize',14,'Box','On','Color',[0.941 0.941 0.941],'EdgeColor',[0 0 0])

title('Data Count Distribution','FontName','Times New Roman','FontWeight','Bold','FontSize',16) %添加图形标题
ylabel('SARS cov-2 Strains','FontName','Times New Roman','FontWeight','Bold','FontSize',15) %给x轴标注
xlabel('Quantity','FontName','Times New Roman','FontWeight','Bold','FontSize',15)%给y轴标注
set(gca,'linewidth',1.5,'FontName','Times New Roman','FontWeight','Bold','FontSize',15,'Box','On','XGrid','on','YGrid','on');
set(gca, 'ytick', 1:length(Y), 'yticklabel',{'MW411947.1', 'MW411948.1', 'MW411949.1', 'MT232662.1', 'MT232664.1', 'MW403692.1', 'MW403693.1', 'MW403694.1', 'MW403695.1', 'MW403696.1', 'MW406485.1', 'MW406487.1', 'MW406488.1', 'MW406490.1', 'MW406491.1', 'MW320729.1', 'MW320730.1', 'MW320731.1', 'MW320733.1', 'MW320735.1', 'MW064259.1', 'MW064260.1', 'MW064263.1', 'MW064264.1', 'MW342706.1', 'MW404672.1', 'MW404673.1', 'MW365356.1', 'MW365357.1'})