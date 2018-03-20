%% time-averaged correlations

% S = something you need to load
clear S2

window = 1;

frameStart = 1;
counter = 1;

for i = 1:window:size(S,3)
    counter = counter + 1;
    S2(:, counter, :) = mean(S(:,frameStart:i,:),2);
    frameStart = i + 1;
end

for j = 1:size(S2,3)
        Corr_DTI2(:,:,j) = corrcoef(S2(:,:,j)');
end
%% Correlations

for i = 1:100
    
    figure(1)
    
    subplot(1,2,1)
    imagesc(Corr_DTI(:,:,i))
    title('Original')
    axis square
    
    subplot(1,2,2)
    imagesc(Corr_DTI2(:,:,i))
    title('Time-averaged')
    axis square
    
    pause()
    
end

%% Average correlations (only when all temps are equal!)

for window = 1:100 
   
    S2 = [];
    Corr_DTI2 = []; 
    
    frameStart = 1;
    counter = 1;
    
    for i = 1:window:size(S,3)
        counter = counter + 1;
        S2(:, counter, :) = mean(S(:,frameStart:i,:),2);
        frameStart = i + 1;
    end
    
    for j = 1:size(S2,3)
        Corr_DTI2(:,:,j) = corrcoef(S2(:,:,j)');
    end
    
    subplot(1,2,1)
    imagesc(mean(Corr_DTI,3))
    title('Original')
    axis square
    
    subplot(1,2,2)
    imagesc(mean(Corr_DTI2,3))
    title('Time-averaged')
    axis square
    
    suptitle(num2str(window))
    
    pause(0.1)
end