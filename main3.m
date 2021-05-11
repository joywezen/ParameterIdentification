clear all
clc
%main function
 %setting truth
 p_tr = 6; %true position
 r_tr = 0.25; %true radius
 dir=1;
 coor=20;
 Fx=50;
 Fy=50;
 %sample amount;
 N = 100;
 %############generate prior geometry guess distribution#############
 %setting prior geometry guess interval
 p_inv = [5,15]; %interval of position guess 
 r_inv = [0.1,0.3]; %interval of radius guess
 %geometry guess mean
 p = mean(p_inv);
 r = mean(r_inv);
 %geometry guess standard deviation
 sig_p = 1; %cannot be too small otherwise tends to be deterministic
 sig_r = 1;
 %initialize geometry samples as normal distribution
 p_d = makedist('Normal',p,sig_p); %mean = p, sigma = sig_r
 r_d = makedist('Normal',r,sig_r);
 %redefine sample distribution as truncate normal distribution
 p_td = truncate(p_d,p_inv(1),p_inv(2));
 r_td = truncate(r_d,r_inv(1),r_inv(2));
 %generate samples following truncated normal distribution
 p_sp = random(p_td,N,1);
 r_sp = random(r_td,N,1);

 %#############virtual prior displacement generating################
 pr_x_mat = zeros(40,N);%initializing 
 pr_y_mat = zeros(40,N);
 for i = 1:N
 [sensor_x,sensor_y] = Forward(p_sp(i),r_sp(i),i,dir,coor,Fx,Fy);
 pr_x_mat(:,i) = sensor_x;
 pr_y_mat(:,i) = sensor_y;
 end
 
 %##############virtual disp measurement generating#################
 [sensor_x_tr,sensor_y_tr] = Forward(p_tr,r_tr,N+1,dir,coor,Fx,Fy); 
  %measurement standard deviation
 sig_ms = abs(min(sensor_x_tr))*0.001; % 0.1% sensor error level
 er_x_mat=randn(40,N)*sig_ms;
 er_y_mat=randn(40,N)*sig_ms;
 ms_x_mat = bsxfun(@plus,er_x_mat,sensor_x_tr);
 ms_y_mat = bsxfun(@plus,er_y_mat,sensor_y_tr);
 
 %########################save Matlab data############################
 save('matlab_5_5_right_up.mat')
 
 %###################### assimilation type 1 #########################
 [tp1_pos_r]=enkf_mat(r_sp,pr_x_mat,er_x_mat,ms_x_mat);
 [tp1_pos_p]=enkf_mat(p_sp,pr_x_mat,er_x_mat,ms_x_mat);
 
 %###################### assimilation type 2 #########################
 [tp2_pos_r]=enkf_mat(r_sp,pr_y_mat,er_y_mat,ms_y_mat);
 [tp2_pos_p]=enkf_mat(p_sp,pr_y_mat,er_y_mat,ms_y_mat);
 
 %###################### assimilation type 3 #########################
 [tp3_pos_r]=enkf_mat(tp1_pos_r,pr_y_mat,er_y_mat,ms_y_mat);
 [tp3_pos_p]=enkf_mat(tp1_pos_p,pr_y_mat,er_y_mat,ms_y_mat);
 
 %##########################plot figure###############################
 
 h1=figure (1);
 
 [f_r_sp,xi_r_sp]=ksdensity(r_sp);
 plot(xi_r_sp,f_r_sp,'linewidth',2,'linestyle','-')
 hold on
 
 [f_tp1_pos_r,xi_tp1_pos_r]=ksdensity(tp1_pos_r);
 plot(xi_tp1_pos_r,f_tp1_pos_r,'linewidth',2,'linestyle','-.')
 
 [f_tp2_pos_r,xi_tp2_pos_r]=ksdensity(tp2_pos_r);
 plot(xi_tp2_pos_r,f_tp2_pos_r,'linewidth',2,'linestyle',':')
 
 [f_tp3_pos_r,xi_tp3_pos_r]=ksdensity(tp3_pos_r);
 plot(xi_tp3_pos_r,f_tp3_pos_r,'linewidth',2,'linestyle','--')
 
 line([r_tr r_tr], get(gca, 'ylim'));
 
 legend('R_{f}','No.1','No.2','No.3','R_{tr}')
 xlabel('R[mm]')
 set(gca,'fontsize',36)
 
 h2=figure (2);
 
 [f_p_sp,xi_p_sp]=ksdensity(p_sp);
 plot(xi_p_sp,f_p_sp,'linewidth',2,'linestyle','-')
 hold on
 
 [f_tp1_pos_p,xi_tp1_pos_p]=ksdensity(tp1_pos_p);
 plot(xi_tp1_pos_p,f_tp1_pos_p,'linewidth',2,'linestyle','-.')
  
 [f_tp2_pos_p,xi_tp2_pos_p]=ksdensity(tp2_pos_p);
 plot(xi_tp2_pos_p,f_tp2_pos_p,'linewidth',2,'linestyle',':')
 
 [f_tp3_pos_p,xi_tp3_pos_p]=ksdensity(tp3_pos_p);
 plot(xi_tp3_pos_p,f_tp3_pos_p,'linewidth',2,'linestyle','--')
 
 line([p_tr p_tr], get(gca, 'ylim'));
 
 legend('P_{f}','No.1','No.2','No.3','P_{tr}')
 xlabel('P(x,5)[mm]')
 set(gca,'fontsize',36)
 