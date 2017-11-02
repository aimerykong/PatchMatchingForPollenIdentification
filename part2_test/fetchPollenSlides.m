clear
close all
clc

%% \
speciesName = 'mixedTestSet';
%sampleFolder = strcat( '../', speciesName, '\', speciesName, '\Picea Mariana B1501' ); % glauca - modern, glauca fossil, mariana - modern, mariana fossil
%sampleFolder = strcat( '../', speciesName ); % glauca - modern, glauca fossil, mariana - modern, mariana fossil
sampleFolder = strcat( './mixedSamples'); % glauca - modern, glauca fossil, mariana - modern, mariana fossil
fprintf('Loading data...');

sampleList = dir(sampleFolder);
for sampleID = 3:length(sampleList) % glauca fossile
    fprintf( '\n\t%d/%d', sampleID-2, length(sampleList)-2 );
    
    subfolder = dir( fullfile(sampleFolder, sampleList(sampleID).name) );
    for subFolderID = 3:length(subfolder) % glauca fossile/NE1010
        subsubFolderList = dir( fullfile(sampleFolder, sampleList(sampleID).name, subfolder(subFolderID).name) );
                
        for subsubFolderID = 3:length(subsubFolderList)
            imglist = dir( fullfile(sampleFolder, sampleList(sampleID).name, subfolder(subFolderID).name, subsubFolderList(subsubFolderID).name) );
            im = imread( fullfile(sampleFolder, sampleList(sampleID).name, subfolder(subFolderID).name, subsubFolderList(subsubFolderID).name,imglist(3).name ));
            mat3D = zeros(size(im,1), size(im,2), length(imglist));
            dataType = whos('im');
            
            fprintf(' (%d%d)', subsubFolderID-2, length(subsubFolderList)-2);
            
            for imId = 3:length(imglist)                
                im = imread( fullfile(sampleFolder, sampleList(sampleID).name, subfolder(subFolderID).name, subsubFolderList(subsubFolderID).name,imglist(imId).name ));
                
                %im = imread( fullfile(sampleFolder, sampleList(sampleID).name, subfolder(subFolderID).name, imglist(imId).name) );
                if strcmp(dataType.class, 'uint8')
                    im = double(im) / (2^8-1);
                else
                    im = double(im) / (2^16-1);
                end
                mat3D(:, :, imId) = im; % image format from uint16 to double
            end
            pollenMaxMap = max(mat3D, [], 3);
            if ~isdir( fullfile( './database', speciesName) )
                mkdir(fullfile( './database', speciesName));
            end
            imwrite(pollenMaxMap, ...
                fullfile( './database', speciesName, strcat(subsubFolderList(subsubFolderID).name, '.jpg' )) );
        end
    end
    fprintf('done\n');
end


%{
    imglist = dir( fullfile(sampleFolder, strcat(sampleList(sampleID).name, '\*.TIF')) );
    im = imread( fullfile(sampleFolder, sampleList(sampleID).name, imglist(1).name) );
    dataType = whos('im');
       
    mat3D = zeros(size(im,1), size(im,2), length(imglist));
    for imId = 1:length(imglist)
        fprintf('.');
        im = imread( fullfile(sampleFolder, sampleList(sampleID).name, imglist(imId).name) );
        if strcmp(dataType.class, 'uint8')
            im = double(im) / (2^8-1);
        else
            im = double(im) / (2^16-1);
        end
        mat3D(:, :, imId) = im; % image format from uint16 to double
    end
    
    pollenMaxMap = max(mat3D, [], 3);
    if ~isdir( fullfile( '.\database', speciesName) )
        mkdir(fullfile( '.\database', speciesName));
    end
    imwrite(pollenMaxMap, ...
        fullfile( '.\database', speciesName, strcat(sampleList(sampleID).name, '.jpg' )) );
end
%}



