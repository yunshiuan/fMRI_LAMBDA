% function [] = cluster_threshold_beta(x_matrix, y_matrix, slices, dim_xy, dim_z, FWHM, dim_resampled, ...
%                                      mask_name, mask_bytes, mask_plot, p_corrected, p_voxel, iterations, savename)
%
% Cluster_threshold_beta written by Scott D. Slotnick 
%
% Copyright 2003-2006 All rights reserved
%
% This program is intended for free use by the neuroimaging community - do not use for financial gain
%
% Modified 12/28/03, for the most recent version, visit http://www2.bc.edu/~slotnics/scripts.htm
%
% If results are used in a publication, please reference -
% Slotnick, S.D., Moo, L.R., Segal, J.B., Hart J. (2003) Distinct prefrontal cortex activity associated with 
% item memory and source memory for visual shapes. Cognitive Brain Research 17, 75-82.
% OR
% Slotnick, S.D. (YEAR). Cluster_threshold. Retrieved MONTH DAY, YEAR, from Web site: 
% http://www2.bc.edu/~slotnics/scripts.htm [MONTH, DAY, YEAR refer to download date]
%
% If you find a bug or have a suggestion, please e-mail sd.slotnick@bc.edu
% 
% This program determines the cluster extent threshold to correct for multiple comparisons in the analysis of neuroimaging
% data. One alternative method is Bonferroni correction for multiple comparisons, which is reasonable if an a priori region-
% of-interest is used, but usually too strict for whole brain analysis. The general idea behind the cluster extent threshold 
% method is to model the entire imaging volume, assume an individual voxel type I error, smooth the volume with a gaussian 
% kernel (if necessary), and then count the number of voxel clusters of each size. After a number of iterations are run, a 
% probability associated with each cluster extent (i.e. number of contiguous voxels) is calculated across runs, and the 
% cluster extent threshold yielding the desired correction for multiple comparisons can be enforced.
% 
% X_MATRIX and Y_MATRIX refer to the acquisition matrix. Number of SLICES is self explanatory
%
% DIM_XY and DIM_Z define the original voxel dimensions in the
% xy/acquisition-plane and z/slice-dimension (in mm) (Original Voxel size)
% 
% FWHM is the full width half maximum of the gassian smoothing kernel (in mm). Set FWHM to 0 to specify no spatial filtering
%
% DIM_RESAMPLED is the resampled voxel resolution (assuming isotropic
% resampling, also in mm) [The normalized voxel size]
%
% MASK_NAME specifies the user defined data file (in analyze '.img' format, e.g. a representative EPI volume) within which 
% the simulation is conducted, of dimensions (X_MATRIX, Y_MATRIX, SLICES). Only voxels with the highest 75% of intensities are 
% maintained as others are assumed to represent noise. For no mask file (recommended), enter 'none' and a rectangular volume 
% will be assumed. Only a small (if any) difference will be obtained between mask vs. no mask, as thresholds are based on voxel
% population statistics. Still, the option to enter a user defined mask is included to entertain the detail oriented user. 
% A MASK_BYTES setting of 1 swaps bytes in user specified mask, if needed. A MASK_PLOT setting of 1 displays image.
%
% P_CORRECTED is the desired correction for multiple comparisons (i.e. the overall type I error) [use .05 or smaller]
% This is the p_FDR reported.

% P_VOXEL is the assumed voxel type I error [use .05 or smaller, .01 or
% lower recommended] (uncorrected p value as specified in SPM).
%
% ITERATIONS refers to the number of monte carlo simulations [25 for preliminary analysis, 1,000 (or 10,000 max) for publication]
%
% SAVENAME is the name of the .mat file in which all the results will be saved
%
% function [] = cluster_threshold_beta(x_matrix, y_matrix, slices, dim_xy, dim_z, FWHM, dim_resampled, ...
%                                      mask_name, mask_bytes, mask_plot, p_corrected, p_voxel, iterations, savename)
%
% Example, cluster_threshold_beta(64,64,30,4.5,4.5,0,3,'none',0,0,.05,.01,1000,'ctb_ofm_pc05_pv01') 
function [] = cluster_threshold_beta(x_matrix, y_matrix, slices, dim_xy, dim_z, FWHM, dim_resampled, mask_name, mask_bytes, mask_plot, p_corrected, p_voxel, iterations, savename)
search_step = 0.1;
search_iterations = 20;
resample_resolution = 4;
iteration_display = 1;
th_step = 1000; 
mask_threshold = .25;
FWHM = FWHM/dim_xy;
FWHM_z = FWHM*dim_xy/dim_z;
randn('state',sum(100*clock))
if strcmp(mask_name,'none')
  F_mask = ones(x_matrix, y_matrix, slices); 
else  
  fid = fopen([mask_name '.img'],'r');
  if mask_bytes
    F = fread(fid, inf, 'uint8');
  else
    F = fread(fid, inf, 'uint16');
  end   
  fclose(fid); 
  if mask_bytes
    F_bit = F(1:end); 
    fhi = F_bit(1:2:end-1); 
    flo = F_bit(2:2:end); 
    Fswap = fhi*256 + flo; 
    F_mask = reshape(Fswap, [x_matrix y_matrix slices 1]);
  else
    F_mask = reshape(F(1:end), [x_matrix y_matrix slices 1]);
  end
end
if strcmp(mask_name,'none') == 0
  figure(1), clf, hold on, colormap('gray')
  imagesc(F_mask(:,:,round(slices/2)))
  title('Center slice of mask. If grainy, quit (<ctrl> c) and switch byte swapping. If ok, press any key. ')
  if mask_plot == 1
    pause
  else
    pause(2)
    drawnow
  end
close
end
if FWHM ~= 0
  sdxy = 2*FWHM*resample_resolution;
  sdz = 2*FWHM_z*resample_resolution;
  g_resampled = zeros(x_matrix*resample_resolution,y_matrix*resample_resolution,slices*resample_resolution);
  for g_search = 1:search_iterations
    guassian_percent_complete = g_search/search_iterations
    for k = 1:slices*resample_resolution
      for j = 1:y_matrix*resample_resolution
        g_resampled(:,j,k) = exp(-(([1:x_matrix*resample_resolution]-x_matrix*resample_resolution/2).^2/sdxy+(j-y_matrix*resample_resolution/2)^2/sdxy+(k-slices*resample_resolution/2)^2/sdz));
      end
    end
    g_resampled = g_resampled./max(max(max(g_resampled)));
    gxy = g_resampled(x_matrix*resample_resolution/2-5*resample_resolution:x_matrix*resample_resolution/2+5*resample_resolution,y_matrix*resample_resolution/2,round(slices*resample_resolution/2));
    mxgxy = find(gxy > .5);
    fwxy = mxgxy(end)-mxgxy(1);
    if fwxy > FWHM*resample_resolution, sdxy = sdxy - sdxy*search_step; end
    if fwxy < FWHM*resample_resolution, sdxy = sdxy + sdxy*search_step; end

    gz = g_resampled(x_matrix*resample_resolution/2, y_matrix*resample_resolution/2, round(slices*resample_resolution/2)-5*resample_resolution:round(slices*resample_resolution/2)+5*resample_resolution);  
    mxgz = find(gz > .5);
    fwz = mxgz(end)-mxgz(1);
    if fwz > FWHM_z*resample_resolution, sdz = sdz - sdz*search_step; end
    if fwz < FWHM_z*resample_resolution, sdz = sdz + sdz*search_step; end
  end
end
for iteration = 1:iterations
  if iteration == 1 & iteration_display > 1, monte_carlo_iteration = iteration, end
  if iteration/iteration_display == round(iteration/iteration_display), monte_carlo_iteration = iteration, end
  brain(:, :, :) = randn(x_matrix, y_matrix, slices);
  brain_resampled = zeros(x_matrix*resample_resolution,y_matrix*resample_resolution,slices*resample_resolution);
  for k = 1:slices
    for i = 1:x_matrix
      for j = 1:y_matrix
        brain_resampled((i-1)*resample_resolution+1:(i-1)*resample_resolution+resample_resolution, (j-1)*resample_resolution+1:(j-1)*resample_resolution+resample_resolution, (k-1)*resample_resolution+1:(k-1)*resample_resolution+resample_resolution) = brain(i, j, k);        
      end
    end
  end
  if FWHM ~= 0
    fb = fftn(brain_resampled);
    fg = fftn(g_resampled);
    fbg = fb.*fg;
    ifbg = real(ifftn(fbg));
  else
    ifbg = brain_resampled;  
  end  
  brain_reconstruct = zeros(x_matrix, y_matrix, slices);
  for k = 1:slices

    for i = 1:x_matrix
      for j = 1:y_matrix
        brain_reconstruct(i, j, k) = mean(mean(mean(ifbg((i-1)*resample_resolution+1:(i-1)*resample_resolution+resample_resolution, (j-1)*resample_resolution+1:(j-1)*resample_resolution+resample_resolution, (k-1)*resample_resolution+1:(k-1)*resample_resolution+resample_resolution))));        
      end
    end
  end
  for k = 1:slices
    F_mask_slice = F_mask(:,:,k);
    F_mask_thr = mask_threshold*max(max(F_mask_slice));
    for i = 1:x_matrix
      for j = 1:y_matrix
        if F_mask_slice(i,j) > F_mask_thr
          brain_mask(i,j,k) = brain_reconstruct(i,j,k);
          brain_mask_count(i,j,k) = 1;
        else
          brain_mask(i,j,k) = 0;
          brain_mask_count(i,j,k) = 0;  
        end
      end
    end
  end
  brain_mask_thr = zeros(size(brain_mask));
  th_init = 0;
  for i = min(min(min(brain_mask))):max(max(max(brain_mask)))
    if (th_init == 0) & (sum(sum(sum(brain_mask > i))))/(sum(sum(sum(brain_mask_count)))) < p_voxel
      th_init = 1;
      brain_mask_thr = brain_mask > i - 1;
    end
  end
  count = sum(sum(sum(brain_mask_thr)));
  p_voxel_computed = count/sum(sum(sum(brain_mask_count)));
  cluster = zeros(x_matrix, y_matrix, slices);
  cluster_count = 1; 
  for i = 1:x_matrix
    for j = 1:y_matrix
      for k = 1:slices
        if brain_mask_thr(i,j,k) == 1
          cluster(i,j,k) = cluster_count;
          if (i > 1 & cluster(i-1,j,k) ~= 0) 
            cluster(cluster == cluster(i-1,j,k)) = cluster_count; 
          end
          if (i < x_matrix & cluster(i + 1,j,k) ~= 0) 
            cluster(cluster == cluster(i+1,j,k)) = cluster_count; 
          end
          if (j > 1 & cluster(i,j-1,k) ~= 0) 
            cluster(cluster == cluster(i,j-1,k)) = cluster_count; 
          end
          if (j < y_matrix & cluster(i,j+1,k) ~= 0) 
            cluster(cluster == cluster(i,j+1,k)) = cluster_count; 
          end
          if (k > 1 & cluster(i,j,k-1) ~= 0) 
            cluster(cluster == cluster(i,j,k-1)) = cluster_count; 
          end

          if (k < slices & cluster(i,j,k+1) ~= 0)
            cluster(cluster == cluster(i,j,k+1)) = cluster_count; 
          end
          cluster_count = cluster_count + 1; 
        end
      end
    end
  end
  cluster_bins = zeros(1, round(p_voxel*x_matrix*y_matrix*slices));
  for i = 1:round(p_voxel*x_matrix*y_matrix*slices)
    cluster_sum = 0;
    for l = 1:slices
      cluster_sum = cluster_sum + sum(sum(cluster(:,:,l) == i));
    end
    cluster_bins(i) = cluster_sum;
  end
  for i = 1:20
    cluster_hist(iteration, i) = sum(cluster_bins == i);
  end
end
cluster_sum = sum(cluster_hist);
cluster_prob = cluster_sum/sum(cluster_sum);
cum_prob(1) = cluster_prob(1);
figure(1), clf, hold on
for i = 2:20
  cum_prob(i) = cluster_prob(i) + cum_prob(i-1);
end
plot(1 - cum_prob)
plot(1:20, 1 - cum_prob, 'bo')
plot([1 20],[p_corrected p_corrected],'r')
xlabel('Cluster size (in original voxels)')
ylabel('1 - p(this cluster extent or smaller)')
title('Cluster extent correction for multiple voxel comparisons')
cluster_indices = find((1 - cum_prob) < p_corrected);
cluster_threshold_original = cluster_indices(1);
cluster_threshold_resampled = round(cluster_threshold_original*dim_xy^2*dim_z/dim_resampled^3 + .49999);
cluster_threshold_volume = cluster_threshold_resampled*dim_resampled^3;
text(1.5, .95, ['Use a cluster extent of ' num2str(cluster_threshold_original) ' original voxels'])
text(1.5, .85, ['which is equivalent to ' num2str(cluster_threshold_resampled) ' resampled voxels']) 
text(1.5, .75, ['both of which define a volume of ' num2str(cluster_threshold_volume) ' cubic mm'])
text(1.5, .65, ['to correct for multiple comparisons at p < ' num2str((round(p_corrected*1000))/1000)])
text(1.5, .55, ['assuming an individual voxel type I error of p = ' num2str((round(p_voxel*1000))/1000)])
axis([1 20 0 1])
eval(['save ' savename ' cluster_threshold_original cluster_threshold_resampled cluster_threshold_volume x_matrix y_matrix slices dim_xy dim_z FWHM dim_resampled mask_name mask_bytes p_corrected p_voxel iterations'])