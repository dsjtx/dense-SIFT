clc;
 clear;

 %% =====================================================================%%

 setDir  = 'D:\Pictures or Images\Dataset\New Dataset (300dpi 24bit depth)';
 imds = imageDatastore(setDir,'IncludeSubfolders',true,'LabelSource',...
     'foldernames');
 [testSet, trainingSet] = splitEachLabel(imds,0.3,'randomize');
 clear i out setDir
%% =====================Calculate Sift Features=======================%%
tic
calculate(testSet,'test');
test = toc;
tic
calculate(trainingSet,'train');
train = toc;
 
function calculate(store,part)
    %
    Txt = fopen(sprintf('%sValue.txt',part),'w');
    nameTxt = fopen(sprintf('%sNameValue.txt',part),'w');
    
 %% ================================Set Parameters=======================%%
    gridSpacing = 16;
    
    
    s='A':'Z';
    [r] = randValue(size(store.Labels,1));
    fprintf('%s Number : ',part);
    for i=1:length(r)
        fprintf('%i ',r(i));
    end
    fprintf('\n');
    
    for i=1:length(r)
        
        Txt = fopen(sprintf('%sValue.txt',part),'a');
        nameTxt = fopen(sprintf('%sNameValue.txt',part),'a');
        
        [~,col] = find(s==char(store.Labels(r(i))));
        
        filename = store.Files{r(i)};
        desc = denseSIFT(filename,gridSpacing);
        
        clc;
        
        filename = split(filename,'\');
        name = char(filename(end));
        fprintf('%i/%i-%s-%s %i',i,length(r),part,name,col);
        
        writeValue(desc,nameTxt,Txt,col,name)
        
        fclose(Txt);
        fclose(nameTxt);    
    end
    
end

function [ r ] = randValue(numfile)
    
    %%using to get random number to writing in text
    rng ('shuffle','v5uniform')
    r = randperm(numfile,numfile);
    cont = 0;
    r = r.';
    tot = unique(r);
    X = length(tot);
    Ncount = histc(r, tot);
    i=1;
    while cont < X
        if Ncount(i) ~= 1
            r(i) = randperm(numfile,1);
            cont=cont-1;
        elseif Ncount(i) == 1
            cont=cont+1;
        end
        
        if cont < X && (i == X || i > X)
            i=1;
        end
        i=i+1;
    end
    
    if cont == numfile
        clc;
    end
end

function writeValue(descriptor,nametext,text,folderLabel,name)

    fprintf(text,'%i ',folderLabel);
    fprintf(nametext,'%i ',folderLabel);
    
    for i=1:length(descriptor)
        %    fprintf(clustDescName,'%i:%5.8f ',i,desript(1,i));
        fprintf(text,'%i:%f ',i,descriptor(i));
        fprintf(nametext,'%i:%f ',i,descriptor(i));

    end
    fprintf('done\n')
    fprintf(text,'\n');
    fprintf(nametext,' %s\n',name);
end