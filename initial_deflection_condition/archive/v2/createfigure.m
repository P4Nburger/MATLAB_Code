function createfigure(X1, Y1)
%CREATEFIGURE(X1, Y1)
%  X1:  plot x データのベクトル
%  Y1:  plot y データのベクトル

%  MATLAB からの自動生成日: 15-Oct-2025 12:17:01

% figure を作成
figure1 = figure;

% axes を作成
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% plot を作成
plot(X1,Y1,'DisplayName','Max Deflection \delta_0 = 2.7200 mm',...
    'LineWidth',2);

% ylabel を作成
ylabel('Deflection w(x) [mm]','FontSize',12);

% xlabel を作成
xlabel('Position x [mm] (Ends at L_{new} = 17.000mm)','FontSize',12);

% title を作成
title('Buckled Shape (Exact Solution, L = 18mm, Push-in = 1mm)',...
    'FontSize',14);

box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% 残りの座標軸プロパティの設定
set(axes1,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',...
    [2.36847946094221 1.5 1]);
% legend を作成
legend(axes1,'show');

