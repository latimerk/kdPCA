rng(12192018);
N = 40;

saveResults = false;
figDir = 'tex/figs_src';
cRange = [-2.0 2.0];
zRange = [0 4];

rRanges = [0    0.5;
           1.0 2.0];
       
M = size(rRanges,1);
       
rs = zeros(N,M);
for ii = 1:M
    rs(:,ii) = rand(N,1)*(rRanges(ii,2)-rRanges(ii,1)) + rRanges(ii,1);
end

psis = rand(N,M)*2*pi;


xs = cos(psis).*rs;
ys = sin(psis).*rs;
zs = (xs.^2+ys.^2);


figure(1);
fontSize_axis = 10;
fontSize_label = 12;
fontSize_title = 12;
clf

colors = [0    0.4470    0.7410;
          0.8500    0.3250    0.0980];
ms = 6;
mm = {'o','x'};
      


subplot(2,2,1)
hold on
for ii = 1:M
    plot(xs(:,ii),ys(:,ii),mm{ii},'Color',colors(ii,:),'MarkerSize',ms);
end

title('original space: ','FontSize',fontSize_title);
xlabel('x','FontSize',fontSize_label);
ylabel('y','FontSize',fontSize_label);
set(gca,'TickDir','out','box','off','FontSize',fontSize_axis,'XTick',[-2 0 2],'YTick',[-2 0 2]);
axis equal
axis square
xlim(cRange);
ylim(cRange);
hold off

subplot(2,2,3);
hold on

[~, SCORE] = pca([xs(:) ys(:)]);
pca_xs = reshape(SCORE(:,1),[],M);
pca_ys = reshape(SCORE(:,2),[],M);
for ii = 1:M
    plot(pca_xs(:,ii),pca_ys(:,ii),mm{ii},'Color',colors(ii,:),'MarkerSize',ms);
end

xlabel('PC 1','FontSize',fontSize_label);
ylabel('PC 2','FontSize',fontSize_label);

set(gca,'TickDir','out','box','off','FontSize',fontSize_axis,'XTick',[-2 0 2],'YTick',[-2 0 2]);
axis equal
axis square
xlim(cRange);
ylim(cRange);
hold off


subplot(2,2,2)
hold on

for ii = 1:M
    plot3(xs(:,ii),ys(:,ii),zs(:,ii),mm{ii},'Color',colors(ii,:),'MarkerSize',ms);
%     plot(xs(:,ii),zs(:,ii),mm{ii},'Color',colors(ii,:),'MarkerSize',ms);
end

title('higher dimensional space: ','FontSize',fontSize_title);
xlabel('x','FontSize',fontSize_label);
ylabel('y','FontSize',fontSize_label);
zlabel('x^2 + y^2','FontSize',fontSize_label);
set(gca,'TickDir','out','box','off','FontSize',fontSize_axis,'XTick',[-2 0 2],'YTick',[-2 0 2],'ZTick',[0 2 4]);
set(gca,'CameraPosition',[-22.9958  -24.1920   11.2705],'CameraTarget',[0 0 2]);
axis equal
axis square
xlim(cRange);
ylim(cRange);
zlim(zRange);
hold off


subplot(2,2,4);
hold on

[~, SCORE] = pca([xs(:) ys(:) zs(:)]);
kpca_xs = reshape(SCORE(:,1),[],M);
kpca_ys = reshape(SCORE(:,2),[],M);
for ii = 1:M
    plot(kpca_xs(:,ii),kpca_ys(:,ii),mm{ii},'Color',colors(ii,:),'MarkerSize',ms);
end

xlabel('PC 1','FontSize',fontSize_label);
ylabel('PC 2','FontSize',fontSize_label);

set(gca,'TickDir','out','box','off','FontSize',fontSize_axis,'XTick',[-2 0 2],'YTick',[-2 0 2]);
axis equal
axis square
xlim(cRange);
ylim(cRange);
hold off


if(saveResults)
    saveas(gcf,sprintf('%s/kernelDemo_raw.eps',figDir),'epsc');
end