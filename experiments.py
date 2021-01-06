import random
import os
from pathlib import Path
import statistics as stats
import numpy as np
from decimal import Decimal
import matplotlib
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod

from vimms.ChemicalSamplers import DatabaseFormulaSampler
from vimms.Chemicals import ChemicalMixtureCreator
from vimms.MassSpec import IndependentMassSpectrometer
from vimms.Controller import SimpleMs1Controller
from vimms.Environment import Environment
from vimms.Common import *

from config import ConfigurationState
from utils import run_xcms, run_mzmine, run_msdial, run_peakonly

#TODO:
#add PeakOnly
#major question: can we break software by messing with scan times, and will interpolating fix that?
#maybe do a more thorough check MZMine time units are different by a factor of 60
#more parameter tuning?
#change saved stuff from generated mzmls to pickle all saved information into one object? - fix rts somehow
#do we want to print success of runs at the end?
#get mzmine centwave?
#can we train ML to get parameters automatically?
#can we train ML on *mistakes* or otherwise use these software (majority vote, include them as a feature, etc)?
#0.6 or 2.6 (10 MS2 * 0.2)
#check ADAP code for using time
#find some way to run MSDial on individual files

#sample unknown chemicals to turn off isotopes
#presenting how many scans per peak there are as a matter of intuition
#check that chemical overlap box in vimms actually matches chemical
#log scan retention time when chemical first appears in a scan as opposed to theoretical RT (vimms) - also log first fragmentation time + intensity - help check between picked and real picks
#disable noise/check it's equal across algorithms
#change wavelet parameters for xcms to check bump at 4.0 on 0.9 overlap
#get peakonly to work

#if no RPath, then we disable runXCMS (and print a warning that this is the reason XCMS was disabled)
#R also needs access to: library(xcms) and library(magrittr) - how do we ensure this? is there R package manager? do we just try/except and suggest that these packages need to be installed?
#Check peaks against chemical number to check parameters, same number of peaks from each one also good sanity check

class PPFile():

    class Peak():
        def __init__(self, min_rt, rt, max_rt, min_mass, max_mass, unit_conv=1):
            units = Decimal(unit_conv)
            self.min_rt = Decimal(min_rt) * units
            self.rt = Decimal(rt) * units
            self.max_rt = Decimal(max_rt) * units
            self.min_mass = None if min_mass is None else Decimal(min_mass)
            self.max_mass = None if max_mass is None else Decimal(max_mass)
            
        def __repr__(self): return "({}, {}, {})".format(self.min_rt, self.rt, self.max_rt)
            
        def compare_to_chem(self, chem, overlap_ratio):
            intersecting = max(0, min(chem["max_rt"], self.max_rt) - max(chem["min_rt"], self.min_rt))
            union = (chem["max_rt"] - chem["min_rt"] + self.max_rt - self.min_rt - intersecting)
            mass_bounds = self.min_mass is None or self.max_mass is None
            if(not mass_bounds):
                minm, maxm = self.min_mass * (1 - overlap_ratio * Decimal(0.1)), self.max_mass * (Decimal(1) + overlap_ratio * Decimal(0.1))
                mass_bounds = mass_bounds or (chem["m/z"] >= minm and chem["m/z"] <= maxm)
            return intersecting >= union * Decimal(overlap_ratio) and mass_bounds
            
    def __init__(self, peaks):
        self.peaks = peaks

    def __repr__(self): return str(self.peaks)
        
    @staticmethod
    def load_data(filepath, delim=','):
        with open(filepath, 'r') as f:
            headers = f.readline().strip().split(delim)
            peaks = [{headers[i] : x for i, x in enumerate(line.split(delim))} for line in f]
        return peaks
            
    @classmethod
    def from_xcms(cls, filepath):
        peaks = cls.load_data(filepath)
        def to_peak(p):
            left_header = next(k for k in p.keys() if k.endswith("Peak RT start\""))
            right_header = next(k for k in p.keys() if k.endswith("Peak RT end\""))
            min_mass = next(k for k in p.keys() if k.endswith("Peak m/z min\""))
            max_mass = next(k for k in p.keys() if k.endswith("Peak m/z max\""))
            return cls.Peak(p[left_header], p["\"row retention time\""], p[right_header], p[min_mass], p[max_mass])
        return PPFile([to_peak(p) for p in peaks])
    
    @classmethod
    def from_mzmine(cls, filepath):
        peaks = cls.load_data(filepath)
        def to_peak(p):
            left_header = next(k for k in p.keys() if k.endswith("Peak RT start"))
            right_header = next(k for k in p.keys() if k.endswith("Peak RT end"))
            min_mass = next(k for k in p.keys() if k.endswith("Peak m/z min"))
            max_mass = next(k for k in p.keys() if k.endswith("Peak m/z max"))
            return cls.Peak(p[left_header], p["row retention time"], p[right_header], p[min_mass], p[max_mass], unit_conv=60)
        return PPFile([to_peak(p) for p in peaks])
    
    @classmethod
    def from_msdial(cls, filepath):
        peaks = cls.load_data(filepath, delim='\t')
        def to_peak(p): return cls.Peak(p["RT left(min)"], p["RT (min)"], p["RT right (min)"], None, None, unit_conv=60)
        return PPFile([to_peak(p) for p in peaks])
    
    @staticmethod
    def canonise_rt(chem): 
        return {
                "min_rt" : Decimal(chem.rt) + Decimal(chem.chromatogram.min_rt), 
                "max_rt" : Decimal(chem.rt) + Decimal(chem.chromatogram.max_rt),
                "m/z" : Decimal(chem.mass)
               }
    
    def compare_to_chems(self, chems, overlap_ratio=0.9):
        TP, FP = 0, 0
        chems = chems.copy()
        for p in self.peaks:
            found = False
            for ch in chems:
                if p.compare_to_chem(ch, overlap_ratio):
                    chems.remove(ch)
                    TP += 1
                    found = True
                    break
            if(not found): FP += 1
        missed = len(chems)
        return TP, FP, missed

class Experiment():
    
    c = ConfigurationState()
    file_counter = 0
    name = None
    xlabel = None
    
    #choose a colour for each software package, and use a lighter shade for TP, and darker for FP
    xcms_tp_colour = (108.0 / 255, 211.0 / 255, 1.0)
    xcms_fp_colour = (0.0, 122.0 / 255, 174.0 / 255)
    mzmine_tp_colour = (1.0, 0.0, 0.0)
    mzmine_fp_colour = (174.0 / 255, 0.0, 0.0)
    msdial_tp_colour = (17.0 / 255, 1.0, 83.0 / 255)
    msdial_fp_colour = (0.0, 147.0 / 255, 40.0 / 255)

    @classmethod
    def generate_chems(cls, output_dir, no_chems):
        hmdb = load_obj(cls.c.HMDBPATH)
        df = DatabaseFormulaSampler(hmdb, min_mz=100, max_mz=1000)
        cm = ChemicalMixtureCreator(df, adduct_prior_dict={POSITIVE: {"M+H" : 1}})
        chemicals = cm.sample(no_chems, 1)
        min_rt, max_rt = min(chem.rt for chem in chemicals) * 0.9, max(chem.rt for chem in chemicals) * 1.1
        Path(output_dir).mkdir(exist_ok=True)
        with open(os.path.join(output_dir, "rts.txt"), 'w') as rts: rts.write("{},{}".format(min_rt, max_rt))
        save_obj(chemicals, os.path.join(output_dir, "chems.pkl"))
        return chemicals, os.path.join(output_dir, "chems.pkl"), os.path.join(output_dir, "rts.txt")

    def __init__(self, chems, min_rt, max_rt, filenames=None):
        self.chems = chems
        self.filenames = None if filenames is None else list(filenames)
        self.min_rt, self.max_rt = min_rt, max_rt

    @staticmethod
    @abstractmethod
    def time_gen(params): pass
    
    @staticmethod
    @abstractmethod
    def xtick_gen(params): pass
    
    def save_params(self, output_dir, params):
        save_obj(list(params), os.path.join(output_dir, "{}_params.pkl".format(cls.name)))

    def generate_mzmls(self, output_dir, params):
        scan_duration_dicts = self.time_gen(params)
        if(self.filenames is None): self.filenames = [os.path.join(output_dir, "time_exp_data_{:04d}.mzML".format(i)) for i, _ in enumerate(scan_duration_dicts)]
        self.file_counter += len(scan_duration_dicts)
        if(len(params) != len(self.filenames)): raise ValueError("Parameter and filename list not the same length!") 
        
        for f, d in zip(self.filenames, scan_duration_dicts):
            mass_spec = IndependentMassSpectrometer(POSITIVE, self.chems, None, scan_duration_dict=d)
            controller = SimpleMs1Controller()

            env = Environment(mass_spec, controller, self.min_rt, self.max_rt, progress_bar=True)
            set_log_level_warning()
            env.run()

            set_log_level_warning()
            env.write_mzML(output_dir, os.path.basename(f))

    def run_peak_pickers(self, xcms_params=None, mzmine_params=None, msdial_params=None, peakonly_params=None, clean=False):
        xcms_filenames, mzmine_filenames, msdial_filenames = None, None, None
        if(not xcms_params is None): xcms_filenames = run_xcms(self.c, self.filenames, xcms_params, clean=clean)
        if(not mzmine_params is None): mzmine_filenames = run_mzmine(self.c, self.filenames, mzmine_params, self.min_rt, self.max_rt, clean=clean)
        if(not msdial_params is None): msdial_filenames = run_msdial(self.c, msdial_params, self.min_rt, self.max_rt, clean=clean)
        return xcms_filenames, mzmine_filenames, msdial_filenames
        
    def plot(self, plotdir, xcms_results, mzmine_results, msdial_results, xticks, oratio):
        xcms_tp, xcms_fp, xcms_missed = list(zip(*xcms_results))
        mzmine_tp, mzmine_fp, mzmine_missed = list(zip(*mzmine_results))
        msdial_tp, msdial_fp, msdial_missed = list(zip(*msdial_results))
        xs = range(max(len(xcms_results), len(mzmine_results), len(msdial_results)))
        
        matplotlib.use('TKAgg') #qt is broken on one of my systems - this is workaround
        fig, ax = plt.subplots()
        ax.scatter(xs, xcms_tp, color=self.xcms_tp_colour, label="XCMS TP")
        ax.scatter(xs, mzmine_tp, color=self.mzmine_tp_colour, label="MZMine TP")
        ax.scatter(xs, msdial_tp, color=self.msdial_tp_colour, label="MSDial TP")
        ax.scatter(xs, xcms_fp, color=self.xcms_fp_colour, label="XCMS FP", marker='x')
        ax.scatter(xs, mzmine_fp, color=self.mzmine_fp_colour, label="MZMine FP", marker='x')
        ax.scatter(xs, msdial_fp, color=self.msdial_fp_colour, label="MSDial FP", marker='x')
        
        ax.set(xlabel=self.xlabel, ylabel="True/False Positives", title="{} True/False Positives {} Overlap".format(self.name, oratio))
        ax.set_ylim(ymin=0)
        ax.legend()
        
        plt.xticks(xs, xticks, rotation=25)
        plt.tight_layout()
        plt.savefig(os.path.join(plotdir, "{}_positives_{}.jpg".format(oratio, self.name)))
        plt.show()
        
        fig, ax = plt.subplots()
        ax.scatter(xs, xcms_missed, color=self.xcms_tp_colour, label="XCMS Missed")
        ax.scatter(xs, mzmine_missed, color=self.mzmine_tp_colour, label="MZMine Missed")
        ax.scatter(xs, msdial_missed, color=self.msdial_tp_colour, label="MSDial Missed")
        
        ax.set(xlabel=self.xlabel, ylabel="Chemicals Missed", title="{} Chemicals Missed {} Overlap".format(self.name, oratio))
        ax.set_ylim(ymin=0)
        ax.legend()
        
        plt.xticks(xs, xticks, rotation=25)
        plt.tight_layout()
        plt.savefig(os.path.join(plotdir, "{}_missed_{}.jpg".format(oratio, self.name)))
        plt.show()

    def match_peaks(self, plotdir, params, oratios, xcms_filenames=None, mzmine_filenames=None, msdial_filenames=None, clean=False):
        chems = [PPFile.canonise_rt(ch) for ch in self.chems]
    
        xcms_files = [PPFile.from_xcms(f) for f in xcms_filenames] if not xcms_filenames is None else []
        mzmine_files = [PPFile.from_mzmine(f) for f in mzmine_filenames] if not mzmine_filenames is None else []
        msdial_files = [PPFile.from_msdial(f) for f in [f for f in msdial_filenames if not os.path.basename(f).startswith("AlignResult")]] if not msdial_filenames is None else []
        
        Path(plotdir).mkdir(exist_ok=True)
        if(clean): 
            for f in os.listdir(plotdir): os.remove(os.path.join(plotdir, f))
        for oratio in oratios:
            xcms_results = [f.compare_to_chems(chems, overlap_ratio=Decimal(oratio)) for f in xcms_files]
            mzmine_results = [f.compare_to_chems(chems, overlap_ratio=Decimal(oratio)) for f in mzmine_files]
            msdial_results = [f.compare_to_chems(chems, overlap_ratio=Decimal(oratio)) for f in msdial_files]
        
            print("\n---\n")
            print(self.name, oratio)
            print(xcms_results)
            print(mzmine_results)
            print(msdial_results)
            
            self.plot(plotdir, xcms_results, mzmine_results, msdial_results, self.xtick_gen(params), oratio)

class StaticTimes(Experiment):
    name = "Fixed"
    xlabel = "Time Steps"
    
    @staticmethod
    def time_gen(params): return [{1 : v} for v in params]
    
    @staticmethod
    def xtick_gen(params): return [str(v) for v in params]

class UniformTimes(Experiment):
    name = "Uniform"
    xlabel = "Spread"
    
    @staticmethod
    def time_gen(params):
        def ugen(umean, uspread): return lambda: random.uniform(umean-uspread, umean+uspread) if umean <= uspread else random.uniform(0, 2*umean)
        return [{1 : ugen(*p)} for p in params]
        
    @staticmethod
    def xtick_gen(params): return ["%.2f" % uspread for _, uspread in params]

class ChoiceTimes(Experiment):
    name = "Choice"
    xlabel = "(Mean, Variance)"

    @staticmethod
    def time_gen(params):
        def chgen(ls): return lambda: random.choice(ls)
        return [{1 : chgen(p)} for p in params]
        
    @staticmethod
    def xtick_gen(params): return ["({:.2f}, {:.2f})".format(stats.mean(ls), stats.pvariance(ls)) for ls in params]

class RecursiveTimes(Experiment):
    name = "Recursive"
    xlabel = "(Timestep, p)"

    @staticmethod
    def time_gen(params):
        def itgen(v, p): 
            def recurse(v, p): return v + recurse(v, p) if p >= random.uniform(0, 1) else v
            return lambda: recurse(v, p)
        return [{1 : itgen(*p)} for p in params]
        
    @staticmethod
    def xtick_gen(params): return ["({:.2f}, {:.2f})".format(v, p) for v, p in params]

class BinaryStateTimes(Experiment):
    class StateGen():

        NOT_BUSY = 0
        BUSY = 1
        
        def __init__(self, min_rt, scan_times):
            self.state = self.BUSY
            self.rt = min_rt
            self.scan_times = (scan_times["NOT_BUSY"], scan_times["BUSY"])
            
        def scan_time(self):
            return self.scan_times[self.state]
        
        @abstractmethod
        def new_state(self): pass
        
        def __call__(self):
            self.state = self.new_state()
            stime = self.scan_time()
            self.rt += stime
            return stime

class MarkovTimes(BinaryStateTimes):
    name = "Markov"
    xlabel = "(Not Busy Scan Time, Busy Scan Time, p)"

    class StateRandom(BinaryStateTimes.StateGen):

        def __init__(self, min_rt, scan_times, p):
            super().__init__(min_rt, scan_times)
            self.p = p

        def new_state(self): return int(not self.state) if self.p >= random.uniform(0, 1) else self.state

    @classmethod
    def time_gen(cls, params):
        return [{1 : cls.StateRandom(*p)} for p in params]
        
    @staticmethod
    def xtick_gen(params): return ["({:.2f}, {:.2f}, {:.2f})".format(d["NOT_BUSY"], d["BUSY"], p) for _, d, p in params]

class HalfTimes(BinaryStateTimes):
    name = "Half"

    class StateHalves(BinaryStateTimes.StateGen):
        
        def __init__(self, min_rt, max_rt, scan_times):
            super().__init__(min_rt, scan_times)
            self.min_rt, self.max_rt = min_rt, max_rt
            
        def new_state(self): return self.NOT_BUSY if self.rt >= (self.min_rt + self.max_rt) / 2 else self.BUSY
        
    @classmethod
    def time_gen(cls, params):
        return [{1 : cls.StateHalves(*p)} for p in params]
        
    @staticmethod
    def xtick_gen(params): return ["({:.2f}, {:.2f})".format(d["NOT_BUSY"], d["BUSY"]) for _, _, d in params]