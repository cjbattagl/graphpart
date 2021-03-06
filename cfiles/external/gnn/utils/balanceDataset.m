function balanceDataset
global dataSet dynamicSystem
%% if we are in a multiclass problem we cannot balance the dataset (to be done)
if (size(dataSet.trainSet.targets,1)>1)
    warn(0,'Dataset balancing is not yet implemented for multiclass problems. It will be disabled')
    dynamicSystem.config.useBalancedDataset=0;
else

    % class with less examples has 1 in the maskMatrix, the other the suitable factor less than 1.
    % This should be done for trainingSet and validationSet, but not in the testSet

    % If we've tried without balancing the dataset the maskMatrix must be
    % converted back from logical to numerical
    dataSet.trainSet.maskMatrix=+dataSet.trainSet.maskMatrix;
    if isfield(dataSet,'validationSet')
        dataSet.validationSet.maskMatrix=+dataSet.validationSet.maskMatrix;
    end

    indexpos=find(dataSet.trainSet.targets==1);
    indexneg=find(dataSet.trainSet.targets==-1);
    sizepos=size(indexpos,2);
    sizeneg=size(indexneg,2);
    if sizepos<sizeneg
        if sizepos==0 err(0, 'I can''t find any positive example in the training set'); end
        fact=sizepos/sizeneg;
        for i=1:sizeneg
            dataSet.trainSet.maskMatrix(indexneg(i),indexneg(i))=fact;
        end
    else
        if sizeneg==0 err(0, 'I can''t find any negative example in the training set'); end
        fact=sizeneg/sizepos;
        for i=1:sizepos
            dataSet.trainSet.maskMatrix(indexpos(i),indexpos(i))=fact;
        end
    end
    % validationSet
    if isfield(dataSet,'validationSet')
        indexpos=find(dataSet.validationSet.targets==1);
        indexneg=find(dataSet.validationSet.targets==-1);
        sizepos=size(indexpos,2);
        sizeneg=size(indexneg,2);
        if sizepos<sizeneg
            dynamicSystem.config.validationBalancedClass='negative';
            if sizepos==0 err(0, 'I can''t find any positive example in the validation set'); end
            dynamicSystem.config.validationBalancingFact=sizepos/sizeneg;
            for i=1:sizeneg
                dataSet.validationSet.maskMatrix(indexneg(i),indexneg(i))=dynamicSystem.config.validationBalancingFact;
            end
        else
            dynamicSystem.config.validationBalancedClass='positive';
            if sizeneg==0 err(0, 'I can''t find any negative example in the validation set'); end
            dynamicSystem.config.validationBalancingFact=sizeneg/sizepos;
            for i=1:sizepos
                dataSet.validationSet.maskMatrix(indexpos(i),indexpos(i))=dynamicSystem.config.validationBalancingFact;
            end
        end
    end
end
