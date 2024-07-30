clc; % 清除命令行
clear; % 清除工作区
close all; % 关闭所有图形窗口

% 数据输入：N表示数据点数，err1和err2表示两组误差数据
data = [50 3.60E-06 5.14E-06
        60 2.08E-06  3.01E-06
        70 1.30E-06  1.91E-06
        80 8.69E-07  1.29E-06
        90 6.08E-07  9.08E-07
        100 4.42E-07 6.64E-07];

% 提取数据
N = data(:,1); % 数据点数
err1 = data(:,2); % 第一组误差数据
err2 = data(:,3); % 第二组误差数据

% 定义斜率为-3的直线的起始点和比例系数
y_start = 1e-5; % 起始点的y值
order = 3; % 斜率为-3
xx = N; % 使用相同的N值作为x坐标
k = y_start / (xx(1)^(-order)); % 计算比例系数
yy = k * xx.^(-order); % 计算斜率为-3的直线的y值

%% 绘图
figure; % 创建新的图形窗口
hold on; % 保持当前图形窗口，允许在同一窗口中绘制多条曲线

% 绘制第一组误差数据，红色线条，线宽2
plot(N, err1, '*-r', 'LineWidth', 2);

% 绘制第二组误差数据，蓝色线条，线宽2
plot(N, err2, '*-b', 'LineWidth', 2);

% 绘制斜率为-3的直线，黑色虚线加菱形标记，线宽2
plot(xx, yy, 'k--d', 'LineWidth', 2);

hold off; % 释放当前图形窗口

% 设置坐标轴为对数坐标轴
set(gca, 'XScale', 'log', 'YScale', 'log');

% 添加图形标签和标题
xlabel('N'); % x轴标签
ylabel('Error'); % y轴标签
legend('Error 1', 'Error 2', 'Slope -3 Line'); % 图例
