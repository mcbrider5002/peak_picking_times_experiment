import itertools
import os
from pathlib import Path
import glob
import subprocess
import xml.etree.ElementTree

def clean_dir(directory, ext):
    for f in glob.glob(os.path.join(directory, "*.{}".format(ext))): os.remove(f)

class XCMSParams():
    def __init__(self, results_dir, ppm=None, pwlower=None, pwupper=None, snthresh=None, noise=None, prefilterlower=None, prefilterupper=None):
        self.results_dir = results_dir
        self.params = [
                        ("ppm", ppm),
                        ("pwlower", pwlower),
                        ("pwupper", pwupper),
                        ("snthresh", snthresh),
                        ("noise", noise),
                        ("prefilterlower", prefilterlower),
                        ("prefilterupper", prefilterupper),
                      ]
        
    def get_args(self): return ["--{} {}".format(opt, v) for opt, v in self.params if not v is None]

def run_xcms(c, filenames, xcms_params, clean=False):
    if(not c.RUNXCMS): return
    Path(xcms_params.results_dir).mkdir(exist_ok=True)
    if(clean): clean_dir(xcms_params.results_dir, "csv")
    p = subprocess.run([c.RPATH, c.XCMS] + xcms_params.get_args() + [xcms_params.results_dir] + filenames)
    return [os.path.join(xcms_params.results_dir, "{}_xcms_box.csv".format(os.path.basename(f).split('.')[0])) for f in filenames]
    
class MZMineParams():
    def __init__(self, results_dir, param_filename):
        self.results_dir = results_dir
        self.param_filename = param_filename
    
def run_mzmine(c, filenames, mzmine_params, min_rt, max_rt, clean=False):
    if(not c.RUNMZMINE): return
    Path(mzmine_params.results_dir).mkdir(exist_ok=True)
    if(clean):
        clean_dir(mzmine_params.results_dir, "csv")
        clean_dir(mzmine_params.results_dir, "xml")
                  
    et = xml.etree.ElementTree.parse(mzmine_params.param_filename)
    for filename in filenames:
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
                    for f in out_files: f.text = os.path.join(mzmine_params.results_dir, "{}_pp.csv".format(os.path.basename(filename).split('.')[0]))
                            
        new_xml = os.path.join(mzmine_params.results_dir, "{}.xml".format(os.path.basename(filename).split('.')[0]))
        et.write(new_xml)
    
    for f in glob.glob(os.path.join(mzmine_params.results_dir, "*.xml")): subprocess.run([c.MZMINE, f])
    return [os.path.join(mzmine_params.results_dir, "{}_pp.csv".format(os.path.basename(f).split('.')[0])) for f in filenames]
  
class MSDialParams():
    def __init__(self, mode, input_dir, output_dir, param_filename):
        self.mode = "lcmsdda" if mode.lower() == "dda" or mode.lower() == "lcmsdda" else "lcmsdia"
        self.input_dir, self.output_dir, self.param_filename = input_dir, output_dir, param_filename
        self.params = [("-i", input_dir), ("-o", output_dir), ("-m", param_filename)]
    
    def get_args(self): return [self.mode] + list(itertools.chain(*self.params))
  
def run_msdial(c, msdial_params, min_rt, max_rt, clean=False):
    if(not c.RUNMSDIAL): return
    Path(msdial_params.output_dir).mkdir(exist_ok=True)
    if(clean): clean_dir(msdial_params.output_dir, "msdial")
    
    d = {
        "Retention time begin" : str(float(min_rt) / 60),
        "Retention time end" : str(float(max_rt) / 60)
    }
    def replace(ln):
        delim = ln.find(':')
        if(delim < 0 or delim > ln.find('#') and ln.find('#') > -1): return ln
        k, v = ln.split(':')[0], ":".join(ln.split(':')[1:])
        return "{}: {}\n".format(k, d.get(k, v.strip()))
    
    with open(msdial_params.param_filename, 'r') as f:
        lines = [replace(ln) for ln in f]
    with open(msdial_params.param_filename, 'w') as f:
        for ln in lines: f.write(ln)
    
    subprocess.run([c.MSDIAL] + msdial_params.get_args())
    for ext in ["dcl", "aef", "pai2"]: clean_dir(msdial_params.input_dir, ext)
    return [os.path.join(msdial_params.output_dir, filename.split(".")[0] + ".msdial") for filename in os.listdir(msdial_params.input_dir)]
    
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