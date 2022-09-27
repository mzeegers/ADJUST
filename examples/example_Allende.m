DataContainer = load('/export/scratch3/zeegers/AllendeDataset/Allende EDS tomography/test_spectral_dataset.mat');
Data = DataContainer.data;
%%
Data = rescale(Data);

%%


imshow(Data(:,:,853,1))
