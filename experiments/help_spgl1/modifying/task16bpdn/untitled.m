load info_lasso
h = figure;
plot(1:66,log(info_spg.rNorm2),'k','LineWidth',1.2)
hold on;
plot(1:37,log([info_spg.rNorm2(1);info_pqn1.rNorm2]),'b');
legend('spg','pqn')
xlabel('iterations')
ylabel('log(residual)')
axis tight
title('spg vs pqn for a lasso problem')
set(h,'LineWidth',10)