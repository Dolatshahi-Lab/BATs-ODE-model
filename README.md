# BATs-ODE-model
The files contained generate an ordinary differential equation model of tumor-directed cytotoxicity and predicts the rate of cytotoxicity for each BAT subtype. BATs are subtyped by expression of TIGIT and LAG3 as TIGIT+, LAG3+, TIGIT+/LAG3+, and double negative.

RB_Calibration_Driver_pp.m: Driver file for CaliPro. This is the main routine.
RB_eval_model_performance.m: function to compare the model runs against the experimental data and determine which of these runs "passed" (according to user specifications) and which runs "failed". Additionally returns the "passed" run parameters and "fail" run parameters in separate parameter structures.
RB_lhs_ode_predator_prey_ode_c.m: The ODEs for the predator prey model.
RB_lhs_ode_predator_prey_settings_new_c.m:This is a settings file for a simple predator/prey model with 2 equations and 4 parameters. can edit the number of samples (model runs), timepoints to analyze, initial parameter ranges, and specifies the model associated with those parameters.
RB_lhs_ode_run_new_c.m: function that samples parameter space using a Latin Hypercube Samplingstrategy, and runs the model using these parameters.
find_ads.m: function to perform alternative density subtraction to refine parameter space.
find_hdr.m
find_hdrcde.m: based on the hdr function within the hdrcde package written in the R programming language by Rob Hyndman, the statistician who discovered/created highest density region selection. Uses HDR to refine parameter space.
is_numeric_min.m
ks_density_twice.m
lhs_ode_default_output_labels_new.m
lhs_ode_define_run_matrices_new.m
lhs_ode_get_run_settings_new.m
lhs_ode_get_run_settings.m
lhs_ode_run_new_c.m
lhs_ode_unif_new.m
