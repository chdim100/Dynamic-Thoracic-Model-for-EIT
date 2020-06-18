function[imageC]=Gauss_Newton_Reconstruction(N,skipcurr,skipvolt,hyperparameter,H,meascur,geometry,recons,solver,prior,Vref)
%%%Performs Gauss-Newton EIT Reconstruction. Requires the EIDORS library
%%%tool
%%%%Inputs:
%N: Electrode Number
%skipcurr: Current skip-m protocol
%skipvolt: Voltage skip-n protocol
%hyperparameter: regularization hyperparameter (lambda)
%H: measurement frame, a N^2X1 collumn vector
%measucur: if set to 0: only tetrapolar measurements. If set to 1:
%Tetrapolar and bipolar measurements, including current source electrodes
%geometry: reconstruction model geometry. If 1: circular disk. If 2:
%thoracic boundary
%recons: Nodal/element reconstruction. 1 for one step nodal solver with Laplace prior
%2 for element-wise solvers.
% prior: reconstruction algorithm/ prior. Works only with elementwise
% reconstruction here. 1 for Gauss-Newton iterative approach. 2 for
% One-Step NOSER prior approach. 3 for Total Variation iterative approach
%solver: 1 for 'absolute', 2 for 'differential'
%Vref: reference measurement frame, a N^2X1 collumn vector


if meascur==1
    opt={'meas_current'};
else
    opt={'no_meas_current'};
end
if geometry==1
    imdl = mk_common_model ('d2d1c',N);
elseif geometry==2
    if N==16
        imdl=mk_common_model('d2T3',N);
    elseif N==32
        imdl=mk_common_model('d2T3',32);
    elseif N==64
        imdl=mk_common_model('h2T3',64);
    end
end
if N==32
    [st, els]= mk_stim_patterns (N,1,[0 skipcurr+1],[0 skipvolt+1;],opt,0.185);
else
    [st, els]= mk_stim_patterns (N,1,[0 skipcurr+1],[0 skipvolt+1;],opt,10);
end

imdl.fwd_model.stimulation= st;
imdl.fwd_model.meas_select= els;
if solver==1
    imdl.reconst_type = 'absolute';
elseif solver==2
    imdl.reconst_type = 'difference';
    img = mk_image(imdl, 1);
    if isempty(Vref)
        vh = fwd_solve(img);
        H1=vh.meas;
    else
        vh=Vref;
        H1=vh;
    end
end
if N~=16
    H=H(find(imdl.fwd_model.meas_select~=0));

end
imdl.hyperparameter.value = hyperparameter;
if recons==1
    imdl.solve=@nodal_solve;
    imdl.RtR_prior = @prior_laplace;
elseif recons==2
    if prior==1
        imdl.solve = @inv_solve_abs_GN;
        imdl.RtR_prior = @prior_laplace;
        imdl.inv_solve_gn.max_iterations = 120;
    elseif prior==2
        imdl.solve=@inv_solve_diff_GN_one_step;
        imdl.RtR_prior= @prior_noser;
        imdl.jacobian_bkgnd.value= 1;
    elseif prior==3
        imdl.solve=       @inv_solve_TV_pdipm;
        imdl.R_prior=     @prior_TV;
        imdl.parameters.term_tolerance= 1e-3;
        imdl.parameters.keep_iterations= 0;
    else
        while prior~=1&&prior~=2&&prior~=3
            fprintf('Press:\n')
            fprintf('1 for Iterative Gauss Newton\n')
            fprintf('2 for NOSER reconstruction\n')
            fprintf('3 for Total Variation\n')
            msg='Selection: ';
            prior=input(msg);
        end
    end
end
if solver==1
    imageC = inv_solve(imdl,H);
else
    imageC = inv_solve(imdl,H,H1);
end
show=0;
if show==1
    show_fem(imageC,1)
    caxis([-0.275 0.275])
    axis off
end
end