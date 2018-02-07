"""
Run series of Grompy models simultaneously.

The number of simultaneous model runs is controlled by ModelOptions.max_proc in the model parameters file. Note that
each single run will still use the number of cores specified by run-escript (using the -t argument, see the escript
documentation), or the OMP_NUM_THREADS setting in your bashrc file. For example add the following line to your bashrc
file if you want grompy to use 10 threads for each model run: export OMP_NUM_THREADS=10

Example
-------
$ python grompy_parallel.py model_parameters.py

This will run all model experiments specified in the model_parameters.py file, with the number of simultaneous runs
specified by ModelOptions.max_proc

"""

import pdb
import sys
import itertools
import os
import inspect
import datetime
import random
import imp
import time

import numpy as np
import pandas as pd

import multiprocessing

import esys.escript as es
import esys.escript.linearPDEs
import esys.weipa

# local libraries
import lib.grompy_lib as grompy_salt_lib
import lib.mesh_functions as mesh_functions


def info(title):
    print title
    print 'module name:', __name__
    print 'parent process:', os.getppid()
    print 'process id:', os.getpid()


def run_model_scenario_and_analyze_results(output, Parameters, ModelOptions,
                                           mesh_function,
                                           run, model_scenario_name,
                                           scenario_parameters, scenario_param_names,
                                           df, model_output_folder,
                                           scriptdir, scenario_name,
                                           nscenarios, dfo=None):
    #pdb.set_trace()

    # override standard output:
    if not os.path.exists('log'):
        os.makedirs('log')
    sys.stdout = open('log/screen_output_%s_PID_%s.out' % (str(run), str(os.getpid())), "w")
    info(str(run))

    year = 365.25 * 24 * 60 * 60

    model_file_adj = ''

    run_id = 'S%i' % run
    #scenario_name = run_id

    print '-' * 30
    print 'model scenario id %s, run %i of %i' % (run_id, run + 1, nscenarios)
    print '-' * 30

    model_file_adj += 'run%s' % run_id

    # update default parameters in Parameter class
    for scenario_param_name, scenario_parameter in \
            zip(scenario_param_names, scenario_parameters):

        if scenario_parameter is not None:
            # find model parameter name to adjust
            model_param_name = scenario_param_name[:-2]

            print 'updating parameter %s from %s to %s' \
                  % (model_param_name,
                     str(getattr(Parameters, model_param_name)),
                     str(scenario_parameter))

            # update model parameter
            setattr(Parameters, model_param_name, scenario_parameter)

            # and add model param name to filename
            model_file_adj += '_%s_%s' % (str(model_param_name),
                                          str(scenario_parameter))

    # set filename for mesh

    mesh_fn = os.path.join(scriptdir,
                           'model_output',
                           '_%i_%s.msh' % (random.randint(0, 100), '%s.msh' % scenario_name))

    # get names and values of input parameters
    attributes = inspect.getmembers(
        Parameters, lambda attribute: not (inspect.isroutine(attribute)))
    attribute_dict = [attribute for attribute in attributes
                      if not (attribute[0].startswith('__') and
                              attribute[0].endswith('__'))]

    if df is None:
        # get attributes
        attribute_names = [attribute[0] for attribute in attributes
                           if not (attribute[0].startswith('__') and
                                   attribute[0].endswith('__'))]

        # set up pandas dataframe to store model input params
        ind = [0]
        columns = attribute_names
        df = pd.DataFrame(index=ind, columns=columns)

    # store input parameters in dataframe
    for a in attribute_dict:
        if a[0] in df.columns:
            if type(a[1]) is list:
                df.ix[run, a[0]] = str(a[1])
            else:
                df.ix[run, a[0]] = a[1]

    #
    start_time = datetime.datetime.now()
    start_date_str = '%i-%i-%i' % (start_time.day,
                                   start_time.month,
                                   start_time.year)
    start_time_str = '%i:%i:%i' % (start_time.hour,
                                   start_time.minute,
                                   start_time.second)
    df.ix[run, 'start_date'] = start_date_str
    df.ix[run, 'start_time'] = start_time_str

    # run a single model scenario
    model_results = grompy_salt_lib.run_model_scenario(
        scenario_name,
        model_output_folder,
        model_file_adj,
        ModelOptions,
        Parameters,
        mesh_function,
        mesh_fn)

    end_time = datetime.datetime.now()
    runtime = end_time - start_time
    df.ix[run, 'computing_runtime_sec'] = runtime.total_seconds()

    print 'processing model results'
    # get model results
    (mesh, surface, sea_surface, k_vector, P, Conc,
     rho_f, viscosity, h, q, q_abs, nodal_flux,
     Pdiff, Cdiff,
     pressure_differences_max, concentration_differences_max,
     pressure_differences_mean, concentration_differences_mean,
     dts, runtimes, nsteps, output_step,
     boundary_conditions, boundary_fluxes, boundary_flux_stats,
     reached_steady_state) = model_results

    dt = dts[-1]
    runtime = runtimes[-1]

    flux = q

    [specified_pressure_bnd,
     specified_pressure,
     specified_concentration_bnd,
     active_concentration_bnd,
     specified_concentration,
     specified_concentration_rho_f,
     rch_bnd_loc,
     active_seepage_bnd] = boundary_conditions

    [flux_surface_norm,
     land_flux_in,
     land_flux_out,
     submarine_flux,
     submarine_flux_in,
     submarine_flux_out] = boundary_fluxes

    [total_flux_over_surface_norm,
     total_rch_flux,
     total_seepage_flux,
     total_land_flux_in,
     total_land_flux_out,
     total_submarine_flux,
     total_submarine_flux_in,
     total_submarine_flux_out,
     ext_rch,
     ext_seepage,
     ext_inflow,
     ext_outflow,
     ext_inflow_land,
     ext_outflow_land,
     ext_inflow_sea,
     ext_outflow_sea,
     ext_outflow_land_threshold,
     ext_outflow_sea_threshold,
     min_land_flux, max_land_flux,
     min_seepage_flux, max_seepage_flux,
     min_submarine_flux, max_submarine_flux] = boundary_flux_stats

    laststep = output_step - 1

    newcols = ['model_scenario_id',
               'P_min', 'P_max',
               'C_min', 'C_max',
               'h_min', 'h_max',
               'vx_min', 'vx_max', 'vy_min', 'vy_max',
               'max_pressure_change', 'max_concentration_change',
               'runtime', 'nsteps', 'dt_final',
               'total_flux_over_surface',
               'total_rch_flux',
               'total_seepage_flux',
               'total_land_flux_in',
               'total_land_flux_out',
               'total_submarine_flux',
               'total_submarine_flux_in',
               'total_submarine_flux_out',
               'ext_inflow_land',
               'ext_outflow_land',
               'ext_inflow_sea',
               'ext_outflow_sea',
               'ext_outflow_land_exc_threshold',
               'ext_outflow_sea_exc_threshold',
               'min_land_flux', 'max_land_flux',
               'min_seepage_flux', 'max_seepage_flux',
               'min_submarine_flux', 'max_submarine_flux']

    if run == 0:
        df = grompy_salt_lib.add_cols_to_df(df, newcols)
        df['model_scenario_id'] = ''

    # store model results in dataframe
    df.ix[run, 'model_scenario_id'] = run_id
    if model_scenario_name is not None:
        df.ix[run, 'model_scenario_name'] = model_scenario_name
    df.ix[run, 'P_min'] = es.inf(P)
    df.ix[run, 'P_max'] = es.sup(P)
    df.ix[run, 'C_min'] = es.inf(Conc)
    df.ix[run, 'C_max'] = es.sup(Conc)
    df.ix[run, 'h_min'] = es.inf(h)
    df.ix[run, 'h_max'] = es.sup(h)
    df.ix[run, 'vx_min'] = es.inf(flux[0])
    df.ix[run, 'vx_max'] = es.sup(flux[0])
    df.ix[run, 'vy_min'] = es.inf(flux[1])
    df.ix[run, 'vy_max'] = es.sup(flux[1])
    df.ix[run, 'max_pressure_change'] = es.Lsup(Pdiff)
    df.ix[run, 'max_concentration_change'] = es.Lsup(Cdiff)
    df.ix[run, 'runtime'] = runtime
    df.ix[run, 'nsteps'] = nsteps
    df.ix[run, 'dt_final'] = dt
    df.ix[run, 'total_flux_over_surface'] = total_flux_over_surface_norm[1]
    df.ix[run, 'total_rch_flux'] = total_rch_flux[1]
    df.ix[run, 'total_seepage_flux'] = total_seepage_flux[1]
    df.ix[run, 'total_land_flux_in'] = total_land_flux_in[1]
    df.ix[run, 'total_land_flux_out'] = total_land_flux_out[1]
    df.ix[run, 'total_submarine_flux'] = total_submarine_flux[1]
    df.ix[run, 'total_submarine_flux_in'] = total_submarine_flux_in[1]
    df.ix[run, 'total_submarine_flux_out'] = total_submarine_flux_out[1]
    df.ix[run, 'ext_inflow_land'] = ext_inflow_land[1]
    df.ix[run, 'ext_outflow_land'] = ext_outflow_land[1]
    df.ix[run, 'ext_inflow_sea'] = ext_inflow_sea[1]
    df.ix[run, 'ext_outflow_sea'] = ext_outflow_sea[1]
    df.ix[run, 'ext_outflow_land_exc_threshold'] = ext_outflow_land_threshold[1]
    df.ix[run, 'ext_outflow_sea_exc_threshold'] = ext_outflow_sea_threshold[1]
    df.ix[run, 'min_land_flux'] = min_land_flux
    df.ix[run, 'max_land_flux'] = max_land_flux
    df.ix[run, 'min_seepage_flux'] = min_seepage_flux
    df.ix[run, 'max_seepage_flux'] = max_seepage_flux
    df.ix[run, 'min_submarine_flux'] = min_submarine_flux
    df.ix[run, 'max_submarine_flux'] = max_submarine_flux

    # check if model output folder exists and create if not
    if not os.path.exists(model_output_folder):
        os.makedirs(model_output_folder)

    # save model output to vtk file
    if ModelOptions.save_vtk_files is True:
        vtk_folder = os.path.join(model_output_folder, 'vtk_files')

        if not os.path.exists(vtk_folder):
            os.makedirs(vtk_folder)

        fn_VTK = os.path.join(vtk_folder,
                              '%s_final_output.vtu' % model_file_adj)
        print 'saving vtk file of model results: %s' % fn_VTK

        nodata = -99999
        flux_surface_plot = nodal_flux * surface + \
            nodata * es.whereZero(surface)
        if sea_surface is None:
            sea_surface_save = surface
        else:
            sea_surface_save = sea_surface

        esys.weipa.saveVTK(fn_VTK,
                           pressure=P,
                           concentration=Conc,
                           h=h,
                           flux=flux,
                           qx=flux[0],
                           qy=flux[1],
                           kx=k_vector[0],
                           ky=k_vector[1],
                           nodal_flux=nodal_flux,
                           surface=surface,
                           sea_surface=sea_surface_save,
                           nodal_flux_surface=flux_surface_plot,
                           specified_pressure_bnd=specified_pressure_bnd,
                           active_seepage_bnd=active_seepage_bnd,
                           recharge_bnd=rch_bnd_loc,
                           active_concentration_bnd=active_concentration_bnd,
                           flux_surface_norm=flux_surface_norm,
                           land_flux_in=land_flux_in,
                           land_flux_out=land_flux_out,
                           submarine_flux=submarine_flux,
                           submarine_flux_in=submarine_flux_in,
                           submarine_flux_out=submarine_flux_out)

        df.ix[run, 'vtk_filename'] = fn_VTK

    if ModelOptions.save_variables_to_csv is True:

        var_folder = os.path.join(model_output_folder, 'variables')

        if not os.path.exists(var_folder):
            os.makedirs(var_folder)

        q_abs = (flux[0]**2 + flux[1]**2) ** 0.5

        model_vars = [P, Conc, h, q_abs * year, flux[0] * year, flux[1] * year]
        varlabels = ['P', 'conc', 'h', 'v', 'vx', 'vy']

        for var, varlabel in zip(model_vars, varlabels):
            xya, va = grompy_salt_lib.convert_to_array(var)

            filename = os.path.join(var_folder,
                                    '%s_%s_final.csv' % (varlabel,
                                    model_file_adj))
            csv_str = 'x,y,%s\n' % varlabel
            for x, y, vai in zip(xya[:, 0], xya[:, 1], va):
                csv_str += '%0.3f,%0.3f,%0.3e\n' % (x, y, vai)

            print 'writing final variable %s to file %s' % (varlabel, filename)

            fout = open(filename, 'w')
            fout.write(csv_str)
            fout.close()

        # write node and element locations and connectivity to separate file
        mesh_fn = os.path.join(var_folder,
                               'mesh_%s.fly' % model_file_adj)
        mesh.write(mesh_fn)

    # write pressure and concentration differences
    diff_folder = os.path.join(model_output_folder, 'P_and_C_change')
    if not os.path.exists(diff_folder):
        os.makedirs(diff_folder)

    df_diff = pd.DataFrame(columns=['timestep', 'dt', 'runtime',
                                    'pressure_change_mean',
                                    'pressure_change_max',
                                    'concentration_change_mean',
                                    'concentration_change_max'],
                           index=np.arange(nsteps))

    #df_diff['timestep'] = np.arange(nsteps)
    df_diff['pressure_change_mean'] = pressure_differences_mean
    df_diff['pressure_change_max'] = pressure_differences_max
    df_diff['concentration_change_mean'] = concentration_differences_mean
    df_diff['concentration_change_max'] = concentration_differences_max
    df_diff['dt'] = dts
    df_diff['runtime'] = runtimes

    filename = os.path.join(diff_folder,
                            'P_and_C_changes_%s_final.csv' % model_file_adj)

    print 'saving P and C changes to %s' % filename
    df_diff.to_csv(filename, index_label='timestep')

    dfo = df.loc[[run]]

    # save model runs input params and results to .csv file

    ind_results_folder = os.path.join(model_output_folder, 'results_individual_runs')
    if not os.path.exists(ind_results_folder):
        os.makedirs(ind_results_folder)

    today = datetime.datetime.now()
    today_str = '%i-%i-%i' % (today.day, today.month, today.year)
    filename = os.path.join(ind_results_folder,
                            'individual_model_results_%s_%s_run_%s_of_%i.csv'
                            % (scenario_name, today_str, str(run), nscenarios))

    print 'saving model run input & output data to %s' % filename
    dfo.to_csv(filename, index_label='model_run')

    # store output in multiprocessing queue
    output.put(dfo)

    return df


def main():
    #year = 365.25 * 24 * 60.0 * 60.0

    #brackish_conc = 0.035 * 0.99

    # option to actively wait for closure of first open process in list if the max number of simultaneous processes
    # is reached. If False, the script simply waits, with increments specified below
    actively_close_processes = False

    # increments between two checks for the number of active processes (sec)
    wait_increments = 60.0

    #################################
    # import model scenario settings
    #################################

    scriptdir = os.path.dirname(os.path.realpath(__file__))

    print ''

    if len(sys.argv) > 1 and sys.argv[-1][-14:] != 'grompy.py':
        #scenario_name = sys.argv[-1]
        #model_input_subfolder = os.path.join(scriptdir, 'model_input',
        #                                     scenario_name)
        inp_file_loc = os.path.join(scriptdir, sys.argv[-1])

        print 'model input files: ', inp_file_loc

        try:
            model_parameters = imp.load_source('model_parameters', inp_file_loc)
        except IOError:
            msg = 'cannot find parameter file %s' % inp_file_loc
            raise IOError(msg)

        ModelParameters = model_parameters.ModelParameters
        ModelOptions = model_parameters.ModelOptions
        ParameterRanges = model_parameters.ParameterRanges

    else:

        print 'running model input data from file ' \
              'model_input/model_parameters.py'

        from model_input.model_parameters import ModelParameters, ModelOptions
        from model_input.model_parameters import ParameterRanges

    scenario_name = ModelOptions.scenario_name

    # select function to create mesh:
    if ModelParameters.mesh_type == 'coastal':
        #mesh_function = mesh_functions.setup_coastal_mesh
        mesh_function = mesh_functions.setup_coastal_mesh_glover1959
    #elif ModelParameters.mesh_type == 'coastal':
    #    mesh_function = mesh_functions.setup_coastal_mesh_new
    elif ModelParameters.mesh_type == 'rectangle':
        mesh_function = mesh_functions.setup_rectangular_mesh

    # run multiple model scenarios
    print '=' * 35
    print 'running model scenarios:'
    print '=' * 35

    # select folder to save results
    model_output_folder = ModelOptions.model_output_dir

    # check if there are previous runs in the model output folder:
    dfo = None
    firstrun = 0

    if os.path.exists(model_output_folder) and ModelOptions.overwrite_old_results is False:

        pt = model_output_folder

        output_files = [(os.path.getmtime(os.path.join(pt, fn)),
                        os.path.basename(os.path.join(pt, fn)))
                        for fn in os.listdir(model_output_folder)
                        if 'final_model_results' in fn and '.csv' in fn]
        output_files.sort()
        if len(output_files) > 1:
            latest_output_file = os.path.join(model_output_folder,
                                              output_files[-1][1])

            print 'found existing model runs file %s in folder %s' \
                  % (latest_output_file, model_output_folder)
            print 'will append new model results'

            dfo = pd.read_csv(latest_output_file)

            # find last number of model scenarios
            firstrun_lst = dfo['model_scenario_id'].tolist()
            firstrun_str = [f for f in firstrun_lst if str(f)[0] == 'S'][-1]
            firstrun = int(firstrun_str[1:]) + 1

            print 'numbering new model runs starting at %i' % firstrun

    if ModelOptions.model_scenario_list is not 'file':

        # get names of scenario parameters
        # import model parameter ranges

        scenario_param_names_raw = dir(ParameterRanges)
        scenario_param_names = [m for m in scenario_param_names_raw
                                if '__' not in m]

        # get all scenario parameter values
        scenario_parameter_list = [getattr(ParameterRanges, p)
                                   for p in scenario_param_names]

        # construct list with all parameter combinations
        if ModelOptions.model_scenario_list is 'combinations':
            scenario_parameter_combinations = \
                list(itertools.product(*scenario_parameter_list))
        else:
            nscens = np.sum(np.array([len(sp) for sp in scenario_parameter_list
                                      if sp is not None]))
            nparams = len(scenario_parameter_list)
            scenario_parameter_combinations = []

            if ModelOptions.initial_base_run is True:
                # add initial base run with unchanged model parameters
                scenario_parameter_combinations.append([None] * nparams)

            for j, sl in enumerate(scenario_parameter_list):
                if sl[0] is not None:
                    sc = [None] * nparams
                    for sli in sl:
                        sci = list(sc)
                        sci[j] = sli
                        scenario_parameter_combinations.append(sci)

        runs = np.arange(len(scenario_parameter_combinations))
        model_scenario_names = ['S%i' % run for run in runs]

    else:
        # read parameter combinations from file
        scen_file = os.path.join('model_input',
                                 ModelOptions.model_scenario_file)
        scen_df = pd.read_csv(scen_file)
        model_scenario_names = list(scen_df['model_scenario_name'])
        scenario_param_names = list(scen_df.columns)[1:]

        scenario_parameter_combinations = []
        for i in scen_df.index:
            scenario_parameter_combinations.append(list(scen_df.ix[i])[1:])

    #
    if ModelOptions.invert_model_run_order is True:
        print 'reversing order of parameter combinations'
        scenario_parameter_combinations = scenario_parameter_combinations[::-1]

    # read default model parameter file
    Parameters = ModelParameters()

    # get attributes
    attributes = inspect.getmembers(
        Parameters, lambda attribute: not (inspect.isroutine(attribute)))
    attribute_names = [attribute[0] for attribute in attributes
                       if not (attribute[0].startswith('__') and
                               attribute[0].endswith('__'))]

    # set up pandas dataframe to store model input params
    ind = np.arange(len(scenario_parameter_combinations))

    columns = attribute_names
    df = pd.DataFrame(index=ind, columns=columns)

    # set up empty list if only 1 param combination
    if scenario_parameter_combinations == []:
        scenario_parameter_combinations = [[None] * len(scenario_param_names)]

    nscenarios = len(scenario_parameter_combinations)

    # set up queue for multiprocessing
    #output = multiprocessing.Queue()
    output = multiprocessing.Manager().Queue()

    processes = []

    results = []

    print 'number of cores available: ', multiprocessing.cpu_count()

    max_proc = ModelOptions.max_proc

    output_timeout = 3.0

    start_time = time.time()

    # go through all model scenarios
    active_proc = 0
    for run, scenario_parameters in enumerate(scenario_parameter_combinations):


        # read default model parameter file
        Parameters = ModelParameters()
        model_scenario_name = model_scenario_names[run]

        # multproc approach example here:
        # https://stackoverflow.com/questions/2046603/is-it-possible-to-run-function-in-a-subprocess-without-threading-or-writing-a-se
        p = multiprocessing.Process(target=run_model_scenario_and_analyze_results,
                                    args=(output, Parameters,
                                          ModelOptions,
                                          mesh_function, run, model_scenario_name,
                                          scenario_parameters, scenario_param_names,
                                          df, model_output_folder,
                                          scriptdir, scenario_name, nscenarios))

        processes.append(p)

        if (len(multiprocessing.active_children()) -1) < max_proc:
            print 'starting a new process'
            p.start()
            active_proc += 1
            print 'number of active processes is now ', (len(multiprocessing.active_children()) -1), active_proc

        else:
            proc_close_count = 0
            while (len(multiprocessing.active_children()) -1) >= max_proc:

                if actively_close_processes is True:
                    print 'waiting for process to close before starting a new one'

                    if processes[proc_close_count].is_alive() is True:
                        processes[proc_close_count].join()
                        active_proc -= 1

                    proc_close_count += 1
                else:
                    pass
                    #wait

            print 'done, starting a new process'
            p.start()
            active_proc += 1

        output_txt = 'total number of processes to open = %i\n' % len(scenario_parameter_combinations)
        output_txt += 'total open processes = %i\n' % len(processes)
        output_txt += 'number of active processes acc. to multiproc module = %i\n' % (len(multiprocessing.active_children()) -1)
        #output_txt += 'number of active processes acc. to internal count of script = %i\n' % active_proc
        runtime = (time.time() - start_time) / 3600.0
        output_txt += 'time passed = %0.2f hrs\n' % runtime

        try:
            runtime_per_proc = runtime / (len(processes) - len(multiprocessing.active_children()) + 1)
            output_txt += 'time per process = %0.2f hrs\n or %0.0f sec' % (runtime_per_proc, runtime_per_proc * 3600.0)
        except:
            pass
        fout = open('parallel_runs_progress.txt', 'w')
        fout.write(output_txt)
        fout.close()

    print 'number of active processes is ', multiprocessing.active_children(), active_proc
    for i, pf in enumerate(processes):
        if pf.is_alive() is True:
            print 'process %i is alive, waiting to finalize' % i
            pf.join()

    # get the results
    print 'done running models, waiting for results'
    for i, p in enumerate(processes):
        try:
            results.append(output.get(timeout=output_timeout))
        except:
            print 'warning, could not get output for process %i' % i

    # add columns to df
    #dfo = results[0]
    #for output_run, result in enumerate(results[1:]):
    #   dfo = pd.merge()
    dfo = pd.concat(results)

    dfo = dfo[results[0].columns]

    # save model runs input params and results to .csv file
    today = datetime.datetime.now()
    today_str = '%i-%i-%i' % (today.day, today.month, today.year)
    runs_done = len(processes) - len(multiprocessing.active_children()) + 1
    filename = os.path.join(model_output_folder,
                            'final_model_results_%s_%s_%i_of_%i_runs.csv'
                            % (scenario_name, today_str, runs_done, nscenarios))

    # check if file exists already
    if os.path.exists(filename):
        backup_filename = filename + '_backup'
        os.rename(filename, backup_filename)
        print 'moved previous input & output data to %s' % backup_filename

    print 'saved final results as %s' % filename
    dfo.to_csv(filename)

    print 'done'

if __name__ == "__main__":
    main()