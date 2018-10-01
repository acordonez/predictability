function ongrid2=grid1togrid2(ongrid1,ncfile);
% ongrid2=grid1togrid2(ongrid1,ncfile)
%
% Regrids CESM output from one grid to another. This should not
% be run interactively because it takes a long time.
%
% Parameters:
% ongrid2: data on new grid
% ongrid1: data on old grid
% ncfile: file with regridding weights
%
% Ana Ordonez 10/2018

N=nc_varget(ncfile,'dst_grid_dims');
Ntot=N(1)*N(2);

dst=squeeze(nc_varget(ncfile,'row')); % nasa stereo gridcell location
src=squeeze(nc_varget(ncfile,'col')); % gx1v6 pop ocn gridcell location
S=squeeze(nc_varget(ncfile,'S')); % map weight

totwts=zeros(Ntot,1);
totwts(dst)=totwts(dst)+S;

NN=size(ongrid1);
if (length(NN)<3) % just 2 D input 
  tmp=squeeze(ongrid1); 
  tmp=tmp'; 
  tmp=tmp(:);
  ongrid2=zeros(Ntot,1);
  ongrid2(dst)=S.*tmp(src); % multiply each element of the matrix (not matrix multiplication)
  ongrid2=ongrid2./totwts;  % divide each element of the matrix
  ongrid2=reshape(ongrid2,N(1),N(2)); 
  ongrid2=ongrid2'; 
else
  ongrid2=zeros(NN(1),Ntot);
  ongrid1=permute(ongrid1,[1 3 2]);
  for n=1:NN(1)
    size(S);
    size(ongrid1(n,src));
    size(dst);
    ongrid2(n,dst)=S'.*ongrid1(n,src);
    end
    size(ongrid2);
    size(totwts);
    ongrid2=ongrid2./(ones(NN(1),1)*totwts');
    ongrid2=reshape(ongrid2,NN(1),N(1),N(2)); ongrid2=permute(ongrid2,[1 3 2]); 
end


