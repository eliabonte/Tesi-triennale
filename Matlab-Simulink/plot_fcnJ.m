function plot_fcnJ(iterazioni,J_visited, J_min_visited, J_min_est, fig)
color_map = lines(3);
figure(fig)
sb(1) = subplot(2,1,1); hold on;
p_min_visited = plot(iterazioni,J_min_visited,'-o','Color',color_map(3,:),...
         'MarkerEdgeColor',color_map(3,:),'MarkerFaceColor',"w",'MarkerSize',8,'Linewidth',2);
p_visited = plot(iterazioni,J_visited,'o','Color',color_map(2,:),...
         'MarkerEdgeColor',color_map(2,:),'MarkerFaceColor',"w",'MarkerSize',8,'Linewidth',2);
for iter = 1:iterazioni(end)
    plot([iter iter],[J_min_visited(iter) J_visited(iter)],'-k');
end
legend([p_visited p_min_visited],{'current observation','current best'},'Location','northeast','Orientation','Vertical');
ylabel('Observations J');
xlabel('iterations [-]');

sb(2) = subplot(2,1,2); hold on;
p_min_visited = plot(iterazioni,J_min_visited,'-o','Color',color_map(3,:),...
         'MarkerEdgeColor',color_map(3,:),'MarkerFaceColor',"w",'MarkerSize',8,'Linewidth',2);
p_min_est = plot(iterazioni,J_min_est,'o','Color',color_map(1,:),...
         'MarkerEdgeColor',color_map(1,:),'MarkerFaceColor',"w",'MarkerSize',8,'Linewidth',2);      
legend([p_min_visited p_min_est],{'current best','current est. best'},'Location','northeast','Orientation','Vertical');
ylabel('Observations J'); 
xlabel('iterations [-]');

linkaxes(sb,'x');

fig = fig + 1;
end

