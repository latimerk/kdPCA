rng(12192018);
N = 20;

saveResults = false;
figDir = 'tex/figs_src';


rRanges = [0    1;
           1.5 2.5];
M = size(rRanges,1);
       
rs = zeros(N,M);
for ii = 1:M
    rs(:,ii) = rand(N,1)*(rRanges(ii,2)-rRanges(ii,1)) + rRanges(ii,1);
end

psis = rand(N,M)*2*pi;


xs = cos(psis).*rs;
ys = sin(psis).*rs;
zs = sqrt(xs.^2+ys.^2);


figure(1);
fontSize_axis = 10;
fontSize_label = 12;
fontSize_title = 12;
clf

colors = [0    0.4470    0.7410;
          0.8500    0.3250    0.0980];
ms = 6;
mm = {'o','x'};
      


subplot(1,2,1)
hold on
for ii = 1:M
    plot(xs(:,ii),ys(:,ii),mm{ii},'Color',colors(ii,:),'MarkerSize',ms);
end

title('original space: ','FontSize',fontSize_title);
xlabel('x','FontSize',fontSize_label);
ylabel('y','FontSize',fontSize_label);
set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);
axis equal
axis square
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
hold off


subplot(1,2,2)
hold on

for ii = 1:M
    plot3(xs(:,ii),ys(:,ii),zs(:,ii),mm{ii},'Color',colors(ii,:),'MarkerSize',ms);
%     plot(xs(:,ii),zs(:,ii),mm{ii},'Color',colors(ii,:),'MarkerSize',ms);
end

title('higher dimensional space: ','FontSize',fontSize_title);
xlabel('x','FontSize',fontSize_label);
ylabel('y','FontSize',fontSize_label);
zlabel('\sqrt{x^2 + y^2}','FontSize',fontSize_label);
set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);
set(gca,'CameraPosition',[-24.2590 -27.4869 7.6523],'CameraTarget',[-0.1774 0.0727 1.2500]);
axis equal
axis square
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
zlim([0 2.5]);
hold off

if(saveResults)
    saveas(gcf,sprintf('%s/kernelDemo_raw.eps',figDir),'epsc');
end