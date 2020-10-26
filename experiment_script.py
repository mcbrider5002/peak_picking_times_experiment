import random
import sys
import os
import platform
import glob
import itertools
from collections import OrderedDict
from pathlib import Path
import subprocess
import xml.etree.ElementTree
from decimal import Decimal
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from vimms.ChemicalSamplers import DatabaseFormulaSampler
from vimms.Chemicals import ChemicalMixtureCreator
from vimms.MassSpec import IndependentMassSpectrometer
from vimms.Controller import SimpleMs1Controller
from vimms.Environment import Environment
from vimms.Common import *

#TODO:
#add PeakOnly
#maybe do a more thorough check MZMine time units are different by a factor of 60
#what to do about MSDial DIA file?
#parameterise: MSDial running mode, alternative callables for times, plotdir and partitions
#fix that some of the time generating functions aren't joined up properly to config
#do something with removing things in the output dir?
#more parameter tuning?
#change saved stuff from generated mzmls to pickle all saved information into one object?
#do we want to print success of runs at the end?

#if no RPath, then we disable runXCMS (and print a warning that this is the reason XCMS was disabled)
#R also needs access to: library(xcms) and library(magrittr) - how do we ensure this? is there R package manager? do we just try/except and suggest that these packages need to be installed?
#Check peaks against chemical number to check parameters, same number of peaks from each one also good sanity check

'''To run the experiment, we first use ChemicalMixtureCreator to create some random chemicals and pass it a UniformMZFormulaSampler (looks like the simplest option)
I pass this to an IndependentMassSpectrometer along with a scan_duration_dict specifying a callable that generates some uniform times for the MS1 scan (and also test this with a constant time, the mean of the initial range of the random ones)
I specify a SimpleMS1Controller to schedule only MS1 scans, no MS2 scans (Controller is a class for the data acquisition strategy, i.e. how MS1 and MS2 scans are scheduled?)
I then run these in an environment
Repeat, fiddling with the mass spec's callable to generate more or less random scan durations (i.e. making the range of the distribution wider so long as I'm sampling uniformly), seeing if the output gets increasingly different as randomness increases
The specifics of how I do this can be inferred from the demo notebooks
Later, we interpolate this data (probably by running back through the simulator with even time intervals) and see if that causes the differences in peak-picking output to disappear?
'''

class ConfigurationState():
    '''A class to hold configuration state for running peak-picking software.
    Class variables contain default values that will be written to config.txt on initialisation, if no such file already exists.
    If a value for a variable can be found in config.txt, that value will be stored in the class instance, otherwise the default in class variables will be used.'''

    WRKDIR = os.path.abspath(os.getcwd())
    
    vimms_params = OrderedDict([
        ("DATADIR", os.path.join(WRKDIR, "data")),
        ("HMDBPATH", os.path.join(WRKDIR, "vimms", "tests", "fixtures", "hmdb_compounds.p")),
        ("CHEMICALS", "100")
    ])

    #times to use for generating .mzMLs
    times = OrderedDict([
        #("FIXEDTIMES", "0.1,0.25,0.55,0.75,1,1.5,2"),
        ("FIXEDTIMES", "1.0,2.0,4.0,8.0,16.0,32.0"),
        ("UNIFORMMEANS", "7.0"),
        ("UNIFORMSPREADS", ",".join(str(7.0 / 20 * i) for i in range(20))) #distance of boundaries of uniform distribution from the mean
    ])

    #params for xcms
    xcms_params = OrderedDict([
        ("RPATH", "Rscript"),
        ("XCMS", os.path.join(WRKDIR, "run_xcms.r")),
        ("XCMSPPM", 15),
        ("XCMSPWLOWER", 15),
        ("XCMSPWUPPER", 80),
        ("XCMSSNTHRESH", 5),
        ("XCMSNOISE", 1000),
        ("XCMSPREFILTERLOWER", 3),
        ("XCMSPREFILTERUPPER", 500),
        ("RUNXCMS", 1)
    ])

    #params for mzmine
    mzmine_params = OrderedDict([
        ("MZMINEPATH", WRKDIR),
        ("MZMINEPARAM", os.path.join(WRKDIR, "mzmine_template.xml")),
    ])
    _globs = glob.glob(os.path.join(mzmine_params["MZMINEPATH"], "MZmine-*-*"))
    _MZMINEDIST = _globs[0] if _globs else None
    if(not _MZMINEDIST is None):
        mzmine_params.update([
            ("MZMINEVER", os.path.basename(_MZMINEDIST).split('-')[1]),
            ("RUNMZMINE", 1)
        ])
    else:
        mzmine_params.update([
            ("MZMINEVER", "0.00"),
            ("RUNMZMINE", 0)
    ])

    #params for msdial
    msdial_params = OrderedDict([
        ("MSDIALPATH", WRKDIR),
        ("MSDIALPARAM", os.path.join(WRKDIR, "msdialparam_dda.txt")),
    ])
    _globs = glob.glob(os.path.join(msdial_params["MSDIALPATH"], "MSDIAL ver *"))
    _MSDIALDIST = _globs[0] if _globs else None
    if(not _MSDIALDIST is None):
        msdial_params.update([
            ("MSDIALVER", os.path.basename(_MSDIALDIST).split()[-1]),
            ("RUNMSDIAL", 1)
        ])
    else:
        msdial_params.update([
            ("MSDIALVER", "0.00"),
            ("RUNMSDIAL", 0)
    ])

    #params for peakonly
    peakonly_params = OrderedDict([
        ("RUNPEAKONLY", 1),
    ])
    
    @classmethod
    def generate_configs(cls, config_path):
        with open(config_path, 'w') as configs:
            def write_variable(name, value): configs.write("{} = {}\n".format(name, value))
            
            for k, v in cls.vimms_params.items(): write_variable(k, v)
            configs.write("\n")
            for k, v in cls.times.items(): write_variable(k, v)
            configs.write("\n")
            for k, v in cls.xcms_params.items(): write_variable(k, v)
            configs.write("\n")
            for k, v in cls.mzmine_params.items(): write_variable(k, v)
            configs.write("\n")
            for k, v in cls.msdial_params.items(): write_variable(k, v)
            configs.write("\n")
            for k, v in cls.peakonly_params.items(): write_variable(k, v)
    
    #try/except may be useful here
    @staticmethod
    def load_configs(config_path):
        with open(config_path, 'r') as configs:
            return {ln.split('=')[0].strip() : ln.split('=')[1].strip() for ln in configs if ln.strip()}
    
    def __init__(self):
    
        config_path = os.path.join(self.WRKDIR, "config.txt")
        #if(not os.path.isfile(config_path)): self.generate_configs()
        self.generate_configs(config_path)
        configs = self.load_configs(config_path)

        system_type = platform.system()
        
        def get_val(params): return lambda k: configs.get(k, params[k])
        get_vimms = get_val(self.vimms_params)
        get_times = get_val(self.times)
        get_xcms = get_val(self.xcms_params)
        get_mzmine = get_val(self.mzmine_params)
        get_msdial = get_val(self.msdial_params)
        
        self.DATADIR = get_vimms("DATADIR")
        self.HMDBPATH = get_vimms("HMDBPATH")
        self.CHEMICALS = int(get_vimms("CHEMICALS"))
        
        self.FIXEDTIMES = [float(t) for t in get_times("FIXEDTIMES").split(",") if t.strip()]
        self.UNIFORMMEANS = [float(m) for m in get_times("UNIFORMMEANS").split(",") if m.strip()]
        self.UNIFORMSPREADS = [float(s) for s in get_times("UNIFORMSPREADS").split(",") if s.strip()]

        self.RPATH = get_xcms("RPATH")
        self.XCMSRESULTS = os.path.join(self.WRKDIR, "xcms_results")
        self.XCMSPARAMS = [
            ("--ppm", get_xcms("XCMSPPM")),
            ("--pwlower", get_xcms("XCMSPWLOWER")),
            ("--pwupper", get_xcms("XCMSPWUPPER")),
            ("--snthresh", get_xcms("XCMSSNTHRESH")),
            ("--noise", get_xcms("XCMSNOISE")),
            ("--prefilterlower", get_xcms("XCMSPREFILTERLOWER")),
            ("--prefilterupper", get_xcms("XCMSPREFILTERUPPER"))
        ]
        self.XCMS = get_xcms("XCMS")
        self.RUNXCMS = get_xcms("RUNXCMS")

        self.MZMINEPATH = get_mzmine("MZMINEPATH")
        self.MZMINERESULTS = os.path.join(self.WRKDIR, "mzmine_results")
        self.MZMINEVER = configs.get("MZMINEVER", None)
        if(not self.MZMINEVER is None):
            self.MZMINEPARAM = get_mzmine("MZMINEPARAM")
            self.MZMINEDIST = os.path.join(self.MZMINEPATH, "MZMINE-{}-{}".format(self.MZMINEVER, system_type))
            self.RUNMZMINE = get_mzmine("RUNMZMINE")
            self.MZMINE = glob.glob(os.path.join(self.MZMINEDIST, "startMZmine*"))[0]
        else:
            self.RUNMZMINE = 0

        self.MSDIALPATH = get_msdial("MSDIALPATH")
        self.MSDIALRESULTS = os.path.join(self.WRKDIR, "msdial_results")
        self.MSDIALVER = configs.get("MSDIALVER", None)
        if(not self.MSDIALVER is None):
            self.MSDIALPARAM = get_msdial("MSDIALPARAM")
            self.MSDIALDIST = os.path.join(self.MSDIALPATH, "MSDIAL ver {}".format(self.MSDIALVER))
            self.RUNMSDIAL = get_msdial("RUNMSDIAL")
            self.MSDIAL = os.path.join(self.MSDIALDIST, "MsdialConsoleApp.exe")
        else:
            self.RUNMSDIAL = 0

        self.RUNPEAKONLY = configs.get("RUNPEAKONLY", 0)
        
        self.PLOTDIR = os.path.join(self.WRKDIR, "plots")
        self.PARTITION_NAMES = None
        self.PARTITIONS = None
        self.XLABELS = None
        self.XTICKS = None

def generate_mzmls(c):
    hmdb = load_obj(c.HMDBPATH)
    df = DatabaseFormulaSampler(hmdb, min_mz=100, max_mz=1000)
    cm = ChemicalMixtureCreator(df, adduct_prior_dict={POSITIVE: {"M+H" : 1}})
    chemicals = cm.sample(c.CHEMICALS,1)
    min_rt, max_rt = min(chem.rt for chem in chemicals) * 0.9, max(chem.rt for chem in chemicals) * 1.1
    for f in os.listdir(c.DATADIR): os.remove(os.path.join(c.DATADIR, f))
    Path(c.DATADIR).mkdir(exist_ok=True)
    with open(os.path.join(c.DATADIR, "rts.txt"), 'w') as rts: rts.write("{},{}".format(min_rt, max_rt))
    save_obj(chemicals, os.path.join(c.DATADIR, "chems.pkl"))

    def ugen(umean, uspread): return lambda: random.uniform(umean-uspread, umean+uspread)
    uniforms = [ugen(umean, uspread) for umean in c.UNIFORMMEANS for uspread in c.UNIFORMSPREADS if umean >= uspread]
    def chgen(ls): return lambda: random.choice(ls)
    choices = [chgen(ls) for ls in [[1] * i + [3] + [5] * i for i in range(1, 6)]]
    def itgen(v, p): return v + itgen(v, p) if p >= random.uniform(0, 1) else v
    iteratives = [itgen(v, p) for v in [1, 3, 5] for p in [0.1, 0.3, 0.5, 0.7, 0.9]]
    
    c.PARTITION_NAMES = ["Fixed", "Uniform", "Choice", "Recursive"]
    c.PARTITIONS = [0] + [len(c.FIXEDTIMES), len(uniforms), len(choices), len(iteratives)]
    c.XLABELS = ["Time Steps", "Spread", "Padding Width", "(Timestep, p)"]
    c.XTICKS = [[str(f) for f in c.FIXEDTIMES], ["%.2f" % u for u in c.UNIFORMSPREADS], list(range(1, 6)), [(v, p) for v in [1, 3, 5] for p in [0.1, 0.3, 0.5, 0.7, 0.9]]]
    scan_duration_dicts = [{1 : v} for v in (c.FIXEDTIMES + uniforms + choices + iteratives)]
    
    for i, d in enumerate(scan_duration_dicts):
        mass_spec = IndependentMassSpectrometer(POSITIVE, chemicals, None, scan_duration_dict=d)
        controller = SimpleMs1Controller()

        env = Environment(mass_spec, controller, min_rt, max_rt, progress_bar=True)
        set_log_level_warning()
        env.run()

        set_log_level_debug()
        env.write_mzML(c.DATADIR, "time_exp_data_%04d.mzML" % i)
        
def run_xcms(c, min_rt, max_rt):
    if(not c.RUNXCMS): return
    Path(c.XCMSRESULTS).mkdir(exist_ok=True)
    for f in os.listdir(c.XCMSRESULTS): os.remove(os.path.join(c.XCMSRESULTS, f))
    data = [f for f in glob.glob(os.path.join(c.DATADIR, "*.mzML"))]
    p = subprocess.run([c.RPATH, c.XCMS] + list(itertools.chain(*c.XCMSPARAMS)) + [c.XCMSRESULTS] + data)
    
def run_mzmine(c, min_rt, max_rt):
    if(not c.RUNMZMINE): return
    Path(c.MZMINERESULTS).mkdir(exist_ok=True)
    for f in os.listdir(c.MZMINERESULTS): os.remove(os.path.join(c.MZMINERESULTS, f))
                  
    et = xml.etree.ElementTree.parse(c.MZMINEPARAM)
    for filename in (f for f in glob.glob(os.path.join(c.DATADIR, "*.mzML"))):
        root = et.getroot()
        for child in root:
        
            if child.attrib["method"].endswith("RawDataImportModule"):
                for e in child:
                    for g in e: g.text = filename
                    
            if child.attrib["method"].endswith("CropFilterModule"):
                for e in child:
                    if e.attrib["name"] == "Scans":
                        ie = e.find("retention_time")
                        if(not ie is None):
                            if(not ie.find("min") is None): ie.find("min").text = str(float(min_rt) / 60)
                            if(not ie.find("max") is None): ie.find("max").text = str(float(max_rt) / 60)
                    
            if child.attrib["method"].endswith("CSVExportModule"):
                for e in child:
                    out_files = itertools.chain(*[e.findall(s) for s in ["current_file", "last_file"]])
                    for f in out_files: f.text = os.path.join(c.MZMINERESULTS, "{}_pp.csv".format(os.path.basename(filename).split('.')[0]))
                            
        new_xml = os.path.join(c.MZMINERESULTS, "{}.xml".format(os.path.basename(filename).split('.')[0]))
        et.write(new_xml)
    
    for f in glob.glob(os.path.join(c.MZMINERESULTS, "*.xml")): subprocess.run([c.MZMINE, f])
    
def run_msdial(c, min_rt, max_rt):
    if(not c.RUNMSDIAL): return
    Path(c.MSDIALRESULTS).mkdir(exist_ok=True)
    for f in os.listdir(c.MSDIALRESULTS): os.remove(os.path.join(c.MSDIALRESULTS, f))
    
    d = {
        "Retention time begin" : str(float(min_rt) / 60),
        "Retention time end" : str(float(max_rt) / 60)
    }
    def replace(ln):
        delim = ln.find(':')
        if(delim < 0 or delim > ln.find('#') and ln.find('#') > -1): return ln
        k, v = ln.split(':')[0], ":".join(ln.split(':')[1:])
        return "{}: {}\n".format(k, d.get(k, v.strip()))
    
    with open(c.MSDIALPARAM, 'r') as f:
        lines = [replace(ln) for ln in f]
    with open(c.MSDIALPARAM, 'w') as f:
        for ln in lines: f.write(ln)
    
    demos = os.path.join(c.MSDIALDIST, "MsdialConsoleApp demo files")
    mode = "lcmsdda"
    #input = os.path.join(demos, "LCMS_DIA")
    input = c.DATADIR
    output = c.MSDIALRESULTS
    #method = os.path.join(demos, "LCMS_DIA", "Msdial-lcms-dia-Param.txt")
    method = c.MSDIALPARAM
    subprocess.run([c.MSDIAL, mode,  "-i", input,  "-o", output, "-m", method])
    
    for ext in ["dcl", "aef", "pai2"]:
        format = "*.{}".format(ext)
        for f in glob.glob(os.path.join(c.DATADIR, format)): os.remove(f)
    
def run_peakonly(c, min_rt, max_rt):
    if(not c.RUNPEAKONLY): return
    
    #from peakonly-master.processing_utils.runner import FilesRunner
    #import PeakOnly.processing_utils.runner as peakonly
    #from PeakOnly.processing_utils.runner import FilesRunner
    '''import importlib.util
    #spec = importlib.util.spec_from_file_location("processing_utils.runner", os.path.join(c.WRKDIR, "PeakOnly", "processing_utils", "runner.py"))
    spec = importlib.util.spec_from_file_location("processing_utils.runner", os.path.join(c.WRKDIR, "PeakOnly"))
    peakonly = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(peakonly)'''
    '''import PeakOnly
    
    runner = PeakOnly.processing_utils.runner.FilesRunner("sequential", [classifier, segmentator], 0.0, 0, 0, 3, torch.device("cuda:0" if torch.cuda.is_available() else "cpu"))
    features = runner(os.listdir(DATADIR))
    print(features)'''

def match_peaks(c, chems): 
    class PPFile():
    
        class Peak():
            def __init__(self, min_rt, rt, max_rt, unit_conv=1):
                units = Decimal(unit_conv)
                self.min_rt = Decimal(min_rt) * units
                self.rt = Decimal(rt) * units
                self.max_rt = Decimal(max_rt) * units
                
            def __repr__(self): return "({}, {}, {})".format(self.min_rt, self.rt, self.max_rt)
                
            def compare_to_chem(self, chem, overlap_ratio):
                intersecting = max(0, min(chem["max_rt"], self.max_rt) - max(chem["min_rt"], self.min_rt))
                union = (chem["max_rt"] - chem["min_rt"] + self.max_rt - self.min_rt - intersecting)
                return intersecting >= union * Decimal(overlap_ratio)    
                
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
                return cls.Peak(p[left_header], p["\"row retention time\""], p[right_header])
            return PPFile([to_peak(p) for p in peaks])
        
        @classmethod
        def from_mzmine(cls, filepath):
            peaks = cls.load_data(filepath)
            def to_peak(p):
                left_header = next(k for k in p.keys() if k.endswith("Peak RT start"))
                right_header = next(k for k in p.keys() if k.endswith("Peak RT end"))
                return cls.Peak(p[left_header], p["row retention time"], p[right_header], unit_conv=60)
            return PPFile([to_peak(p) for p in peaks])
        
        @classmethod
        def from_msdial(cls, filepath):
            peaks = cls.load_data(filepath, delim='\t')
            def to_peak(p): return cls.Peak(p["RT left(min)"], p["RT (min)"], p["RT right (min)"], unit_conv=60)
            return PPFile([to_peak(p) for p in peaks])
        
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
    
    def canonise_rt(chem): return {"min_rt" : Decimal(chem.rt) + Decimal(chem.chromatogram.min_rt), "max_rt" : Decimal(chem.rt) + Decimal(chem.chromatogram.max_rt)}
    chems = [canonise_rt(ch) for ch in chems]
    
    oratio = 0.7
    #load all XCMS results, compare pp files to chemical
    files = [PPFile.from_xcms(f) for f in glob.glob(os.path.join(c.XCMSRESULTS, "*.csv"))]
    xcms_results = [f.compare_to_chems(chems, overlap_ratio=oratio) for f in files]
    #load all MZMine results, compare pp files to chemical
    files = [PPFile.from_mzmine(f) for f in glob.glob(os.path.join(c.MZMINERESULTS, "*.csv"))]
    mzmine_results = [f.compare_to_chems(chems, overlap_ratio=oratio) for f in files]
    #load all MSDial results, compare pp files to chemical
    files = [PPFile.from_msdial(f) for f in [f for f in glob.glob(os.path.join(c.MSDIALRESULTS, "*.msdial")) if not os.path.basename(f).startswith("AlignResult")]]
    msdial_results = [f.compare_to_chems(chems, overlap_ratio=oratio) for f in files]
    
    print(xcms_results)
    print(mzmine_results)
    print(msdial_results)
    
    #choose a colour for each software package, and use a lighter shade for TP, and darker for FP
    xcms_tp_colour = (108.0 / 255, 211.0 / 255, 1.0)
    xcms_fp_colour = (0.0, 122.0 / 255, 174.0 / 255)
    mzmine_tp_colour = (1.0, 0.0, 0.0)
    mzmine_fp_colour = (174.0 / 255, 0.0, 0.0)
    msdial_tp_colour = (17.0 / 255, 1.0, 83.0 / 255)
    msdial_fp_colour = (0.0, 147.0 / 255, 40.0 / 255)
        
    def plot(name, i, i2, xlabel, xticks):
        xs = list(range(i2 - i))
        xcms_tp, xcms_fp, xcms_missed = list(zip(*xcms_results[i:i2]))
        mzmine_tp, mzmine_fp, mzmine_missed = list(zip(*mzmine_results[i:i2]))
        msdial_tp, msdial_fp, msdial_missed = list(zip(*msdial_results[i:i2]))
        
        fig, ax = plt.subplots()
        ax.scatter(xs, xcms_tp, color=xcms_tp_colour, label="XCMS TP")
        ax.scatter(xs, mzmine_tp, color=mzmine_tp_colour, label="MZMine TP")
        ax.scatter(xs, msdial_tp, color=msdial_tp_colour, label="MSDial TP")
        ax.scatter(xs, xcms_fp, color=xcms_fp_colour, label="XCMS FP", marker='x')
        ax.scatter(xs, mzmine_fp, color=mzmine_fp_colour, label="MZMine FP", marker='x')
        ax.scatter(xs, msdial_fp, color=msdial_fp_colour, label="MSDial FP", marker='x')
        
        ax.set(xlabel=xlabel, ylabel="True/False Positives", title="{} True/False Positives".format(name))
        ax.set_ylim(ymin=0)
        ax.legend()
        
        plt.xticks(xs, xticks, rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(c.PLOTDIR, "{}_positives.jpg".format(name)))
        plt.show()
        
        fig, ax = plt.subplots()
        ax.scatter(xs, xcms_missed, color=xcms_tp_colour, label="XCMS Missed")
        ax.scatter(xs, mzmine_missed, color=mzmine_tp_colour, label="MZMine Missed")
        ax.scatter(xs, msdial_missed, color=msdial_tp_colour, label="MSDial Missed")
        
        ax.set(xlabel=xlabel, ylabel="Chemicals Missed", title="{} Chemicals Missed".format(name))
        ax.set_ylim(ymin=0)
        ax.legend()
        
        plt.xticks(xs, xticks, rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(c.PLOTDIR, "{}_missed.jpg".format(name)))
        plt.show()
    
    Path(c.PLOTDIR).mkdir(exist_ok=True)
    partitions = list(itertools.accumulate(c.PARTITIONS))
    for name, i, i2, xlbl, xt in zip(c.PARTITION_NAMES, partitions, partitions[1:], c.XLABELS, c.XTICKS):
        plot(name, i, i2, xlbl, xt)

def main():
    config = ConfigurationState()
    #if(not os.path.isdir(config.DATADIR) or not os.listdir(config.DATADIR)): generate_mzmls(config)
    generate_mzmls(config)
    with open(os.path.join(config.DATADIR, "rts.txt"), 'r') as rts: min_rt, max_rt = rts.readline().split(',')
    run_xcms(config, min_rt, max_rt)
    run_mzmine(config, min_rt, max_rt)
    run_msdial(config, min_rt, max_rt)
    #run_peakonly(config, min_rt, max_rt)
    chemicals = load_obj(os.path.join(config.DATADIR, "chems.pkl"))
    match_peaks(config, chemicals)
    
main()