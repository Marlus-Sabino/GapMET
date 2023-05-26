import importlib
import subprocess

def is_module_installed(module_name):
    try:
        importlib.import_module(module_name)
        return True
    except ImportError:
        return False

def install_package(package_name):
    subprocess.check_call(['pip', 'install', package_name])

def check_libraries():
    packages = ['numpy', 'wheel', 'pandas', 'scikit-learn', 'pymap3d', 'statsmodels', 'scipy']
    for package_name in packages:
        if not is_module_installed(package_name):
            print(f"Installing {package_name}...")
            install_package(package_name)
            print(f"{package_name} is now installed.")
        else:
            print(f"{package_name} is already installed.")
