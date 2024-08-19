OUTBASE = './';
name = 'Fr and Ft';
this = [10:10:30 40:20:80 100]*100;
h = figure('visible','off'); hold on
for i=1:length(this)
    plot(r*R(this(i)),Fr(:,this(i)),'linewidth',max(5-0.75*i, 0.75));
end
for i=1:length(this)  
    plot(r*R(this(i)),Ft(:,this(i)),'--','linewidth',max(5-0.75*i, 0.75));
end

ylabel(name);
set(gca,'fontsize',16)
xlabel('R/L');
name = strrep(strrep(name,'\',''),'/',' over ');
print(h,'-dpng',[OUTBASE name '.png'],'-r300'); close(h);

%make legends
    legends = {'Fr','Ft'};

    h = figure('visible','off'); hold on
    plot(1,1,'linewidth',max(5-0.75*i, 0.75));
    plot(1,1,'--','linewidth',max(5-0.75*i, 0.75));
    axis off
    legend(legends,'interpreter','latex','location','northwest');
    print(h,'-dpng','-r300',[OUTBASE 'legends']); close(h);
