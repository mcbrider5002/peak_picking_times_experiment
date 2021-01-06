import os
import shutil
from pathlib import Path

from utils import XCMSParams, MZMineParams, MSDialParams
from experiments import Experiment, StaticTimes, UniformTimes, ChoiceTimes, RecursiveTimes, MarkovTimes

def init_dir(dirname):
    Path(dirname).mkdir(exist_ok=True)
    for f in os.listdir(dirname): 
        try: os.remove(os.path.join(dirname, f))
        except: shutil.rmtree(os.path.join(dirname, f))
    return dirname

def main():
    wrkdir = os.path.abspath(os.getcwd())
    datadir = init_dir(os.path.join(wrkdir, "data"))
    xcms_outs = init_dir(os.path.join(wrkdir, "xcms_results"))
    mzmine_outs = init_dir(os.path.join(wrkdir, "mzmine_results"))
    msdial_outs = init_dir(os.path.join(wrkdir, "msdial_results"))
    plotdir = os.path.join(wrkdir, "plots")
    
    chems, _, rt_path = Experiment.generate_chems(datadir, 500)
    with open(rt_path, 'r') as rts: min_rt, max_rt = map(float, rts.readline().split(','))
    
    def run_experiment(experiment_type, params, oratios, clean=False):
        exp = experiment_type(chems, min_rt, max_rt)
        exp_name = experiment_type.name.lower()
        
        def subdir(maindir): return init_dir(os.path.join(maindir, exp_name))
        xcms_params = XCMSParams(subdir(xcms_outs))
        mzmine_params = MZMineParams(subdir(mzmine_outs), os.path.join(wrkdir, "mzmine_template.xml"))
        msdial_params = MSDialParams("lcmsdda", subdir(datadir), subdir(msdial_outs), os.path.join(wrkdir, "msdialparam_dda.txt"))
        
        exp.generate_mzmls(subdir(datadir), params)
        xcms_filenames, mzmine_filenames, msdial_filenames = exp.run_peak_pickers(xcms_params=xcms_params, mzmine_params=mzmine_params, msdial_params=msdial_params)
        def results(): exp.match_peaks(plotdir, params, oratios, xcms_filenames=xcms_filenames, mzmine_filenames=mzmine_filenames, msdial_filenames=msdial_filenames, clean=clean)
        return results
    
    oratios = [0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
    
    static_params = [1.0, 2.0, 4.0, 6.0, 8.0, 16.0, 32.0]
    uniform_params = [(7.0, 7.0 / 20 * i) for i in range(20)]
    choice_params = [ls for ls in [[9] * i + [11] + [13] * i for i in range(1, 101, 10)]]
    recursive_params = [(v, p) for v in [1, 3, 5] for p in [0.1, 0.3, 0.5, 0.7, 0.9]]
    markov_params = [(min_rt, {"NOT_BUSY" : 0.6, "BUSY" : 10}, p) for p in [0.01, 0.03, 0.05, 0.1, 0.2]]
    params = [(StaticTimes, static_params), (UniformTimes, uniform_params), (ChoiceTimes, choice_params), (RecursiveTimes, recursive_params), (MarkovTimes, markov_params)]
    
    result_fns = [run_experiment(expt, expp, oratios, clean=(not i)) for i, (expt, expp) in enumerate(params)]
    for fn in result_fns: fn()
    
if __name__ == "__main__": main()