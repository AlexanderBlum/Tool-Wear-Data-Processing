function z_fit = fit_plunge(dx, z, n)

x = linspace(0, length(z)*dx, length(z))';
p = polyfit(x, z, n);

z_fit = polyval(p, x);

% z_res = z - z_fit;
% 
% % figure
% % plot(x, z);
% % hold on
% % plot(x, z_fit);
% 
% figure
% plot(x, z_res*1000);
% ylabel('deviation from fit line (nm)');
% xlabel('position along tool edge (\mum)');
% axis([-100 100 -100 100]);
end