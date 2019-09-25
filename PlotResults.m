function PlotResults(targets,outputs,Name)

    errors=round((targets-outputs),2);

    RMSE=round(sqrt(mean(errors(:).^2)),2);
    
    error_mean=round(mean(errors(:)),2);
    error_std=round(std(errors(:)),2);

    subplot(2,2,[1 2]);
    plot(targets,'k*');
    xlabel('Pixel#')
    ylabel('LST(K)')
    grid on
    hold on;
    plot(outputs,'ro');
   
    grid on
    legend('Measured LST using IRT','Estimated LST by SW algorithm');
    title(Name);

    subplot(2,2,3);
    bar(errors,0.4,'k');
    axis tight
    grid on
    legend('Error');
    title(['RMSE = ' num2str(RMSE)]);
    ylabel('Error(K)')
    xlabel('Pixel#')

    subplot(2,2,4);
    h=histfit(errors);
    h(1).FaceColor = [0.5 0.5 0.5];
    grid on
    title(['Error mean = ' num2str(error_mean) ', Error std = ' num2str(error_std)]);
    ylabel('Frequency')
    xlabel('Error (K)')
end