import os
import platform
import glob
from collections import OrderedDict

class ConfigurationState():
    '''A class to hold configuration state for running peak-picking software.
    Class variables contain default values that will be written to config.txt on initialisation, if no such file already exists.
    If a value for a variable can be found in config.txt, that value will be stored in the class instance, otherwise the default in class variables will be used.'''

    WRKDIR = os.path.abspath(os.getcwd())
    
    vimms_params = OrderedDict([
        ("HMDBPATH", os.path.join(WRKDIR, "vimms", "tests", "fixtures", "hmdb_compounds.p")),
    ])

    #params for xcms
    xcms_params = OrderedDict([
        ("RPATH", "Rscript"),
        ("XCMS", os.path.join(WRKDIR, "run_xcms.R")),
        ("RUNXCMS", 1)
    ])

    #params for mzmine
    mzmine_params = OrderedDict([
        ("MZMINEPATH", WRKDIR),
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
        if(not os.path.isfile(config_path)): self.generate_configs()
        configs = self.load_configs(config_path)

        system_type = platform.system()
        
        def get_val(params): return lambda k: configs.get(k, params[k])
        get_vimms = get_val(self.vimms_params)
        get_xcms = get_val(self.xcms_params)
        get_mzmine = get_val(self.mzmine_params)
        get_msdial = get_val(self.msdial_params)
        
        self.HMDBPATH = get_vimms("HMDBPATH")

        self.RPATH = get_xcms("RPATH")
        self.XCMS = get_xcms("XCMS")
        self.RUNXCMS = get_xcms("RUNXCMS")

        self.MZMINEPATH = get_mzmine("MZMINEPATH")
        self.MZMINEVER = configs.get("MZMINEVER", None)
        if(not self.MZMINEVER is None):
            self.MZMINEDIST = os.path.join(self.MZMINEPATH, "MZMINE-{}-{}".format(self.MZMINEVER, system_type))
            self.RUNMZMINE = get_mzmine("RUNMZMINE")
            self.MZMINE = glob.glob(os.path.join(self.MZMINEDIST, "startMZmine*"))[0]
        else:
            self.RUNMZMINE = 0

        self.MSDIALPATH = get_msdial("MSDIALPATH")
        self.MSDIALVER = configs.get("MSDIALVER", None)
        if(not self.MSDIALVER is None):
            self.MSDIALDIST = os.path.join(self.MSDIALPATH, "MSDIAL ver {}".format(self.MSDIALVER))
            self.RUNMSDIAL = get_msdial("RUNMSDIAL")
            self.MSDIAL = os.path.join(self.MSDIALDIST, "MsdialConsoleApp.exe")
        else:
            self.RUNMSDIAL = 0

        self.RUNPEAKONLY = configs.get("RUNPEAKONLY", 0)