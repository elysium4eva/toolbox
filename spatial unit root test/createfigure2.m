function createfigure2(X1, YMatrix1)
%CREATEFIGURE2(X1, YMatrix1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 10-Aug-2020 17:11:53

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',1,'Parent',axes1);
set(plot1(1),'DisplayName','W(2,2)','Marker','+');
set(plot1(2),'DisplayName','W(4,4)','Marker','*');
set(plot1(3),'DisplayName','W(6,6)','Marker','o');
set(plot1(4),'DisplayName','W(8,8)','Marker','square');
set(plot1(5),'DisplayName','W(10,10)','Marker','diamond');
set(plot1(6),'DisplayName','W(12,12)','Marker','x');
set(plot1(7),'DisplayName','W(14,14)','Marker','hexagram');
set(plot1(8),'DisplayName','W(16,16)','Marker','pentagram');
set(plot1(9),'DisplayName','W(18,18)','Marker','^');
set(plot1(10),'DisplayName','W(20,20)','Marker','v');

% Create ylabel
ylabel({'Theta'});

% Create xlabel
xlabel({'T'});

box(axes1,'on');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.691806836679329 0.151448672002702 0.163863266413456 0.283711217183771],...
    'EdgeColor',[1 1 1]);

