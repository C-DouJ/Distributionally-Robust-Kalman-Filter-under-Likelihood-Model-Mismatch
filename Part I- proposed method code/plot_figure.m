colorsTab = [
    0 0 0       1.0;   % 第1条 纯黑 (强制置顶)
    1 0.11 0.31   1.0;   % 第11条 勃艮第红 (ΔE=46)
];




figure(2);
j = 1:ts;
plot(j*T,rmse_kf_1(1,:),'color',colorsTab(1,:) ,'linewidth',2.5);hold on;
plot(j*T,rmse_DRKF_1(1,:),'color',colorsTab(2,:) ,'linewidth',2.5);hold on; 
xlabel('Time (s)');
ylabel('RMSE_{pos} (m)');
legend('KF','Proposed','Interpreter', 'latex', 'Location', 'best');

figure(3);
j = 1:ts;
plot(j*T,rmse_kf_2(1,:),'color',colorsTab(1,:) ,'linewidth',2.5);hold on;
plot(j*T,rmse_DRKF_2(1,:),'color',colorsTab(2,:) ,'linewidth',2.5);hold on; 
xlabel('Time (s)');
ylabel('RMSE_{vel} (m/s)');
legend('KF','Proposed','Interpreter', 'latex', 'Location', 'best');

