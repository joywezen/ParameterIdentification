clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%     test A     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %setting truth
 p_tr = 12; %true position
 r_tr = 0.25; %true radius
 dir=1;
 coor=20;
 Fx=50;
 Fy=0;
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
 
 %###################### assimilation type A #########################
 [A_pos_r_x]=enkf_mat(r_sp,pr_x_mat,er_x_mat,ms_x_mat);
 [A_pos_r_xy]=enkf_mat(A_pos_r_x,pr_y_mat,er_y_mat,ms_y_mat);
 [A_pos_p_y]=enkf_mat(p_sp,pr_y_mat,er_y_mat,ms_y_mat);

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%     test B     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setting truth
 p_tr = 12; %true position
 r_tr = 0.25; %true radius
 dir=2;
 coor=10;
 Fx=0;
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
 

 
 %###################### assimilation type B #########################
 [B_pos_r_y]=enkf_mat(A_pos_r_xy,pr_y_mat,er_y_mat,ms_y_mat);
 [B_pos_p_y]=enkf_mat(A_pos_p_y,pr_y_mat,er_y_mat,ms_y_mat);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%     test C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setting truth
 p_tr = 12; %true position
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
 %###################### assimilation type B #########################
 [C_pos_r_x]=enkf_mat(B_pos_r_y,pr_x_mat,er_x_mat,ms_x_mat);
 [C_pos_p_x]=enkf_mat(B_pos_p_y,pr_x_mat,er_x_mat,ms_x_mat);
 [C_pos_p_xy]=enkf_mat(C_pos_p_x,pr_y_mat,er_y_mat,ms_y_mat);
 
 %##########################plot figure###############################
 h1=figure(1);
 
 [f_A_r,xi_A_r]=ksdensity(A_pos_r_xy); 
 [f_B_r,xi_B_r]=ksdensity(B_pos_r_y);
 [f_C_r,xi_C_r]=ksdensity(C_pos_r_x);
 hold on
 plot(xi_A_r,f_A_r,'linewidth',2,'linestyle','-');
 plot(xi_B_r,f_B_r,'linewidth',2,'linestyle','-.');
 plot(xi_C_r,f_C_r,'linewidth',2,'linestyle',':');
 legend('A','A&B','A&B&C')
 xlabel('R[mm]')
 set(gca,'fontsize',36)
 
 h2=figure(2);
 
 [f_A_p,xi_A_p]=ksdensity(A_pos_p_y); 
 [f_B_p,xi_B_p]=ksdensity(B_pos_p_y);
 [f_C_p,xi_C_p]=ksdensity(C_pos_p_xy);
 hold on
 plot(xi_A_p,f_A_p,'linewidth',2,'linestyle','-');
 plot(xi_B_p,f_B_p,'linewidth',2,'linestyle','-.');
 plot(xi_C_p,f_C_p,'linewidth',2,'linestyle',':');
 legend('A','A&B','A&B&C')
 xlabel('P[mm]')
 set(gca,'fontsize',36)
 
