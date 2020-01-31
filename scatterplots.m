opts = detectImportOptions('decol_simulations.csv','NumHeaderLines',1)
a= readtable('decol_simulations.csv',opts);
xv=a{:,4};
yv=a{:,5};
kv=a{:,1};
hv=a{:,2};
fv=a{:,3};


figure(1)
colormap jet
hold on
scatter(xv,yv,10,kv,'o','filled');
colorbar
xlabel('Application Effectiveness')
ylabel('Duration of Protection (Hours)')
legend

figure(2)
colormap jet
hold on
scatter(xv,yv,10,hv,'o','filled');
colorbar
xlabel('Application Effectiveness')
ylabel('Duration of Protection (Hours)')
legend

figure(3)
colormap jet
hold on
scatter(xv,yv,10,fv,'o','filled');
colorbar
xlabel('Application Effectiveness')
ylabel('Duration of Protection (Hours)')
legend
