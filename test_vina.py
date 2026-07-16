from docking_engine import download_vina_executable
import subprocess

def test():
    import urllib.request
    urllib.request.urlretrieve("https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/master/build/linux/release/vina", "vina")
    pass
